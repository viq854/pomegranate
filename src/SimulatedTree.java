/*
 * Program POMEGRANATE for cell lineage tree simulation and sampling
 * by Victoria Popic (viq@stanford.edu) 2014-2015
 *
 * MIT License
 *
 * Copyright (c) 2014 Victoria Popic.
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the "Software"),
 * to deal in the Software without restriction, including
 * without limitation the rights to use, copy, modify, merge, publish,
 * distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject to
 * the following conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
 * BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
 * ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
 * CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
*/

import java.awt.Color;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Random;

/**
 * Simulated cell lineage tree
 * Each node represents a cell population with a newly acquired mutation;  
 * the path from the GL root to the node contains all the mutations present in this cell population.
 * The tree can be sampled using randomized or localized sampling schemes.
 */
public class SimulatedTree {
	
	private ArrayList<CellPopulation> nodes; 
	private HashMap<CellPopulation, ArrayList<CellPopulation>> edges;
	private int numDeadNodes;
	private Random randGen = new Random();
	
	/**
	 * Creates an initial tree with the GL root node
	 */
	public SimulatedTree() {
		nodes = new ArrayList<CellPopulation>();
		edges = new HashMap<CellPopulation, ArrayList<CellPopulation>>();
		CellPopulation germlineRoot = new CellPopulation(); 
		germlineRoot.setGermline();
		nodes.add(germlineRoot);
		numDeadNodes = 0;
	}
	
	/**
	 * For each undead cell population node in the tree
	 * create a descendant population or undergo death (except root)
	 * with some probability
	 */
	public void grow() {
		ArrayList<CellPopulation> children = new ArrayList<CellPopulation>();
		for(CellPopulation node : nodes) {
			if(node.isDead()) continue;
		
			// population death
			float death_roll = randGen.nextFloat();
			if(death_roll < Parameters.PROB_DEATH && !node.isGermline()) {
				node.setDead();
				numDeadNodes++;
				continue;
			}
			
			// division
			Mutation childMut = null;
			float roll = randGen.nextFloat();
			if(roll < Parameters.PROB_SNV) {
				if(node.isCNV() && Parameters.UP_CNV_EFFECT){
					childMut = new Mutation.SNV((Mutation.CNV) node.getLastMutation());
				} else {
					childMut = new Mutation.SNV();
				}
			} else if(roll < (Parameters.PROB_SNV + Parameters.PROB_CNV)) {
				if(node.isCNV() || node.isGermline || !Parameters.UP_CNV_EFFECT) {
					childMut = new Mutation.CNV();
				} else {
					childMut = new Mutation.CNV((Mutation.SNV) node.getLastMutation());
				}
			}
			if(childMut == null) continue;
			
			CellPopulation child = new CellPopulation(); 
			child.setSize(randGen.nextInt(Parameters.MAX_POPULATION_SIZE));
			child.setMutations(node.getMutations()); // all the parent mutations 
			child.addMutation(childMut); // + new mutation
			children.add(child);
			ArrayList<CellPopulation> nbrs = edges.get(node);
			if(nbrs == null) {
				edges.put(node, new ArrayList<CellPopulation>());
			}
			edges.get(node).add(child);
		}
		nodes.addAll(children);
	}
	
	public int getNumNodes() {
		return nodes.size();
	}
	
	public int getNumDeadNodes() {
		return numDeadNodes;
	}
	
	/**
	 * Randomly pick a subset of undead nodes from the node list
	 */
	public ArrayList<CellPopulation> selectSubclones(ArrayList<CellPopulation> nodeList, int maxSubclones) {
		Collections.shuffle(nodeList);
		int numSubclonesToSample = 1 + (maxSubclones > 1 ? randGen.nextInt(maxSubclones-1) : 0);
		ArrayList<CellPopulation> subclones = new ArrayList<CellPopulation>();
		for(int i = 0; i < nodeList.size(); i++) {
			if(subclones.size() >= numSubclonesToSample) break;
			if(nodeList.get(i).isDead) continue;
			if(nodeList.get(i).isGermline) continue;
			subclones.add(nodeList.get(i));
		}
		return subclones;
	}
	
	/**
	 * Extract a sample using randomized sampling
	 */
	public TumorSample getSample() {
		ArrayList<CellPopulation> subclones = selectSubclones(nodes, Parameters.MAX_NUM_SUBCLONES);
		return createSample(subclones, Parameters.NUM_CELLS_PER_SAMPLE, getNormalContamination());
	}
	
	private TumorSample createSample(ArrayList<CellPopulation> subclones, int numCellsInSample, int numNormalCells) {
		TumorSample sample = new TumorSample();
		sample.setNumNormalCells(numNormalCells);
		int totalCellCount = 0;
		for(CellPopulation subclone : subclones) {
			subclone.sampleColors.add(sample.color);
			totalCellCount += subclone.size;
		}
		// simulate a multinomial distribution
		int numSubclones = subclones.size();
		float[] subcloneProbabilities = new float[numSubclones];
		float[] intervalLimits = new float[numSubclones];
		subcloneProbabilities[0] = (float) subclones.get(0).size/totalCellCount;
		for(int i = 0; i < numSubclones-1; i++) {
			subcloneProbabilities[i] = (float) subclones.get(i).size/totalCellCount;
			intervalLimits[i] = ((i > 0) ? intervalLimits[i-1] : 0) + subcloneProbabilities[i];
		}
		intervalLimits[numSubclones-1] = 1;

		for(int i = 0; i < numCellsInSample - numNormalCells; i++) {
			float trialCell = randGen.nextFloat();
			// find the interval to classify the cell
			for(int j = 0; j < numSubclones; j++) {
				if(trialCell < intervalLimits[j]) {
					sample.addCell(subclones.get(j));
					break;
				}
			}
		}
		return sample;
	}
	
	/**
	 * Returns the number of cells in the subtree anchored at root,
	 * including the root node
	 */
	private int getSubtreeSize(CellPopulation root) {
		int size = 0;
		ArrayList<CellPopulation> q = new ArrayList<CellPopulation>();
		q.add(root);
		while(q.size() > 0) {
			CellPopulation p = q.remove(0);
			size += (p.isDead ? 0 : p.size);
			if(edges.get(p) != null) {
				q.addAll(edges.get(p));
			}
		}
		root.setSubtreeSize(size);
		return size;
	}
	
	/**
	 * Returns all the nodes in the subtree anchored at root,
	 * including the root node
	 */
	private ArrayList<CellPopulation> getSubtreeNodes(CellPopulation root) {
		ArrayList<CellPopulation> subtree = new ArrayList<CellPopulation>();
		ArrayList<CellPopulation> q = new ArrayList<CellPopulation>();
		q.add(root);
		while(q.size() > 0) {
			CellPopulation p = q.remove(0);
			subtree.add(p);
			if(edges.get(p) != null) {
				q.addAll(edges.get(p));
			}
		}
		return subtree;
	}
	
	/**
	 * Returns k localized samples from the tree
	 * by sampling from k disjoint subtrees.
	 * If the tree does not contain k disjoint subtrees,
	 * sampling will be done from the maximum number of disjoint 
	 * subtrees, with the minimum number of samples overlapping
	 */
	public ArrayList<TumorSample> getKLocalizedSamples(int k) {
		if(edges.get(nodes.get(0)) == null) {
			System.err.println("Cannot collect samples from the tree, only the root node is present");
			System.exit(-1);
		}
		
		// stores the collected samples
		ArrayList<TumorSample> samples = new ArrayList<TumorSample>();
		// queue of disjoint subtree roots (initially undead subtree children of GL root)
		ArrayList<CellPopulation> subtreeRoots = new ArrayList<CellPopulation>();
		for(CellPopulation cp : edges.get(nodes.get(0))) {
			if(getSubtreeSize(cp) > 0) {
				subtreeRoots.add(cp);
			}
		}
		while(subtreeRoots.size() < k) { 
			ArrayList<CellPopulation> children = new ArrayList<CellPopulation>();
			int parentIdx = 0;
			for(int i = 0; i < subtreeRoots.size(); i++) {
				CellPopulation p = subtreeRoots.get(i);
				if(edges.get(p) != null) {
					for(CellPopulation cp : edges.get(p)) {
						if(getSubtreeSize(cp) > 0) {
							children.add(cp);
						}
					}
					if(children.size() > 0) {
						parentIdx = i; 
						break;
					}
				}
			}
			if(children.size() == 0) { // no nodes in the queue have children
				break;
			}
			subtreeRoots.remove(parentIdx);
			subtreeRoots.addAll(0, children);
		} 
		
		Collections.sort(subtreeRoots, new Comparator<CellPopulation>() {
		       public int compare(CellPopulation o1, CellPopulation o2) {
		    	   return o1.getSubtreeSize() < o2.getSubtreeSize() ? 1 : o1.getSubtreeSize() > o2.getSubtreeSize() ? -1 : 0;
		       }
		   });
		
		// define the subtrees
		HashMap<Integer, ArrayList<CellPopulation>> subtrees = new HashMap<Integer, ArrayList<CellPopulation>>();
		for(int i = 0; i < k; i++) {
			CellPopulation root = subtreeRoots.get(i % subtreeRoots.size());
			subtrees.put(i, getSubtreeNodes(root));
		}
		
		// select samples from each subtree
		for(int i = 0; i < k; i++) {
			ArrayList<CellPopulation> subclones = selectSubclones(subtrees.get(i), Parameters.MAX_NUM_SUBCLONES);
			
			// add a subclone from a neighboring subtree
			if(Parameters.MIX_NBR_SUBTREE_SUBCLONE) {
				if(i > 0) {
					subclones.addAll(selectSubclones(subtrees.get(i-1), 1));
				} else {
					subclones.addAll(selectSubclones(subtrees.get(k-1), 1));
				}
			}
			samples.add(createSample(subclones, Parameters.NUM_CELLS_PER_SAMPLE, getNormalContamination()));
		}
		return samples;
	}
	
	/**
	 * Generates a random contamination percentage
	 * based on the provided contamination thresholds
	 */
	private int getNormalContamination() {
		double percentNormal = Parameters.MIN_PERCENT_NORMAL_CONTAMINATION;
		if(Parameters.MAX_PERCENT_NORMAL_CONTAMINATION > Parameters.MIN_PERCENT_NORMAL_CONTAMINATION) {
			percentNormal += new Random().nextDouble()*(Parameters.MAX_PERCENT_NORMAL_CONTAMINATION - Parameters.MIN_PERCENT_NORMAL_CONTAMINATION);
		} 
		int numNormalCells = (int) ((double)(percentNormal*Parameters.NUM_CELLS_PER_SAMPLE)/100.0);
		return numNormalCells;
	}
	
	public void resetColors() {
		for(CellPopulation p : nodes) {
			p.sampleColors = new ArrayList<Color>();
		}
	}
	
	public String toString() {
		String t = "";
		for(CellPopulation n : edges.keySet()) {
			ArrayList<CellPopulation> nbrs = edges.get(n);
			for(CellPopulation n2 : nbrs) {
				if(n.isGermline) {
					t += "GL";
				} else {
					t += n.getName();
				}
				t += "\t" + n2.getName() + "\n";
			}
		}
		
		for(CellPopulation n : nodes) {
			t +=  n.getLastMutation() + "\n";
		}
		return t;
	}
	
	public String toDOT() {
		DecimalFormat df = new DecimalFormat("#.##");
		String t = "";
		t += "digraph G { \n";
		for(CellPopulation n : edges.keySet()) {
			ArrayList<CellPopulation> nbrs = edges.get(n);
			for(CellPopulation n2 : nbrs) {
				t += n.id + " -> " + n2.id + ";\n";
			}
		}
		for(CellPopulation n : nodes) {
			if(!n.isGermline()) {
				String color = "white";
				if(n.isDead) {
					color = "grey";
				}
				if(n.isCNV()) {
					t += n.id + " [shape=star style=filled fillcolor=" + color + " fontname=\"helvetica-bold\" fontsize=42 label=\"" + n.getName() + "\"];\n";
				} else {
					t += n.id + " [shape=circle style=filled fillcolor=" + color + " fontname=\"helvetica-bold\" fontsize=56 label=\"" + n.getName() + "\"" + " width=" + df.format(5*((double)n.size/Parameters.MAX_POPULATION_SIZE)) + " height=2 ];\n";
				}
			} else {
				t += n.id + " [label=\"GL\" fontname=\"arial-bold\" fontsize=56 width=5 height=5];\n";
			}
		}
		t += "}";
		return t;
	}
	
	public String toColoredDOT(ArrayList<TumorSample> samples) {
		DecimalFormat df = new DecimalFormat("#.##");
		String t = "";
		t += "digraph G { \n";
		t += "rankdir=TB;\n";
	    for(CellPopulation n : edges.keySet()) {
			ArrayList<CellPopulation> nbrs = edges.get(n);
			for(CellPopulation n2 : nbrs) {
				t += n.id + " -> " + n2.id + ";\n";
			}
		}
		for(CellPopulation n : nodes) {
			if(!n.isGermline()) {
				String color = "white";
				if(n.isDead) {
					color = "grey";
				} 
				if(n.sampleColors.size() > 0) {
					color = "\"";
					for(int i = 0; i < n.sampleColors.size(); i++) {
						Color c = n.sampleColors.get(i);
						if(i != 0) {
							color += ":";
						}
						color += String.format("#%02x%02x%02x", c.getRed(), c.getGreen(), c.getBlue());
					}
					color += "\"";
				}
				
				if(n.isCNV()) {
					t += n.id + " [shape=star style=filled fillcolor=" + color + " fontname=\"helvetica-bold\" fontsize=42 label=\"" + n.getName() + "\"];\n";
				} else {
					if(n.sampleColors.size() > 1) {
						t += n.id + " [shape=circle style=wedged color=" + color + " fontname=\"helvetica-bold\" fontsize=56 label=\"" + n.getName() + "\"" + " width=" + df.format(5*((double)n.size/Parameters.MAX_POPULATION_SIZE)) +" height=2 ];\n";
					} else {
						t += n.id + " [shape=circle style=filled fillcolor=" + color + " fontname=\"helvetica-bold\" fontsize=56 label=\"" + n.getName() + "\"" + " width=" + df.format(5*((double)n.size/Parameters.MAX_POPULATION_SIZE)) + " height=2 ];\n";
					}
				}
			} else {
				t += n.id + " [label=\"GL\" fontname=\"arial-bold\" fontsize=56 width=5 height=5];\n"; 
			}
		}
		
		t += "{";
		t += "rank=sink;\n";
		t += "Legend[shape=none, margin=0, label=";
		t += "<<TABLE border=\"0\" cellborder=\"0\" cellspacing=\"0\"> \n";
		t += "<TR>";
		for(int i = 1; i <= samples.size(); i++) {
			Color c = samples.get(i-1).color;
			String color = "\""+String.format("#%02x%02x%02x", c.getRed(), c.getGreen(), c.getBlue())+"\"";
			t += "<TD width=\"200\" height=\"200\" colspan=\"1\"><FONT POINT-SIZE=\"36.0\"><B>Sample " + i + "</B></FONT></TD><TD width=\"200\" height=\"200\" colspan=\"1\" BGCOLOR="+ color +"></TD>\n";
		}
		t += "</TR>";
		t += "</TABLE>>];\n";
	    t += "} \n";
	    
	    t += "}";
		return t;
	}

}
