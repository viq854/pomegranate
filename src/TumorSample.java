/*
 * Program POMEGRANATE for cell lineage tree simulation and sampling
 * by Victoria Popic (viq@stanford.edu) 2014
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
import java.util.HashMap;
import java.util.Random;

/**
 * Represents a tumor sample consisting of several cell populations
 */
public class TumorSample {

	protected HashMap<CellPopulation, Integer> cellPopulationCounts;
	protected int numNormalCells; // normal contamination 
	protected int numCNVAffectedSNVs;
	
	protected Color color;
	
	public TumorSample() {
		numNormalCells = 0;
		cellPopulationCounts = new HashMap<CellPopulation, Integer>();
		Random r = new Random();
		color = new Color(r.nextFloat(), r.nextFloat(), r.nextFloat());
	}
	
	public void addCell(CellPopulation cell) {
		if(cellPopulationCounts.containsKey(cell)) {
			int count = cellPopulationCounts.get(cell).intValue();
			cellPopulationCounts.put(cell, count + 1);
		} else {
			cellPopulationCounts.put(cell, 1);
		}
	}
	
	public void setNumNormalCells(int numCells) {
		numNormalCells = numCells;
	}
	
	public int getNumSubclones() {
		return cellPopulationCounts.keySet().size();
	}
	
	public HashMap<Mutation.SNV, Double> getSNVFrequencies() {
		// count how many cells contain each mutation
		if(Parameters.PROB_CNV == 0) {
			HashMap<Mutation, Integer> snvCounts = new HashMap<Mutation, Integer>();
			int totalNumCells = 0;
			for(CellPopulation c : cellPopulationCounts.keySet()) {
				for(Mutation snv : c.getMutations()) {
					if(snvCounts.containsKey(snv)) {
						int count = snvCounts.get(snv).intValue();
						snvCounts.put(snv, count + cellPopulationCounts.get(c));
					} else {
						snvCounts.put(snv, cellPopulationCounts.get(c));
					}
				}
				totalNumCells += cellPopulationCounts.get(c);
			}
			totalNumCells += numNormalCells;
					
			HashMap<Mutation.SNV, Double> freq = new HashMap<Mutation.SNV, Double>();
			for(Mutation snv : snvCounts.keySet()) {
				freq.put((Mutation.SNV) snv, (double) snvCounts.get(snv)/(2*totalNumCells));
			}
			return freq;
		}
		
		HashMap<Mutation.SNV, Integer> ref_haplotype_counts = new HashMap<Mutation.SNV, Integer>();
		HashMap<Mutation.SNV, Integer> var_haplotype_counts = new HashMap<Mutation.SNV, Integer>();
		HashMap<Mutation.SNV, Boolean> affectedSNVs = new HashMap<Mutation.SNV, Boolean>();
		for(CellPopulation c : cellPopulationCounts.keySet()) {
			int[] haplotype_ref_counts = new int[c.getMutations().size()]; // how many copies of the reference exit
			int[] haplotype_var_counts = new int[c.getMutations().size()]; // how many copies of the SNV exist
			for(int i = 0; i < c.getMutations().size(); i++) {
				Mutation m = c.getMutations().get(i);
				if(m instanceof Mutation.CNV) {
					continue;
				}
				Mutation.SNV snv = (Mutation.SNV) m;
				// find all the CNVs affecting this locus
				for(int j = 0; j < c.getMutations().size(); j++) {
					Mutation m2 = c.getMutations().get(j); 
					if(m2 instanceof Mutation.SNV) continue;
					Mutation.CNV cnv = (Mutation.CNV) m2;
					if(cnv.chr == snv.chr) {
						if((snv.position <= Mutation.CHROMOSOME_LENGTHS[snv.chr]/2 && cnv.arm == 0) || 
								(snv.position > Mutation.CHROMOSOME_LENGTHS[snv.chr]/2 && cnv.arm == 1)) {
							// matched the arm
							affectedSNVs.put(snv, true);
							if(j < i) {
								// CNV happened before the SNV occurred 
								haplotype_ref_counts[i]++;
							} else {
								if(cnv.haplotype == snv.haplotype) {
									haplotype_var_counts[i]++;
								} else {
									haplotype_ref_counts[i]++;
								}
							}
						}
					}
				}
				
				if(ref_haplotype_counts.containsKey(snv)) {
					int count_ref = ref_haplotype_counts.get(snv).intValue();
					ref_haplotype_counts.put(snv, count_ref + cellPopulationCounts.get(c)*(haplotype_var_counts[i] + haplotype_ref_counts[i] + 2));
					int count_var = var_haplotype_counts.get(snv).intValue();
					var_haplotype_counts.put(snv, count_var + cellPopulationCounts.get(c)*(haplotype_var_counts[i] + 1));
				} else {
					ref_haplotype_counts.put(snv, cellPopulationCounts.get(c)*(haplotype_var_counts[i] + haplotype_ref_counts[i] + 2));
					var_haplotype_counts.put(snv, cellPopulationCounts.get(c)*(haplotype_var_counts[i] + 1));
				}
			}
		}
		for(Mutation.SNV s : ref_haplotype_counts.keySet()) {
			for(CellPopulation c : cellPopulationCounts.keySet()) {
				boolean contains = false;
				int numCNVs = 0;
				for(int i = 0; i < c.getMutations().size(); i++) {
					Mutation m = c.getMutations().get(i);
					if(m instanceof Mutation.CNV) {
						Mutation.CNV cnv = (Mutation.CNV) m;
						if(cnv.chr == s.chr) {
							if((s.position <= Mutation.CHROMOSOME_LENGTHS[s.chr]/2 && cnv.arm == 0) || 
									(s.position > Mutation.CHROMOSOME_LENGTHS[s.chr]/2 && cnv.arm == 1)) {
								numCNVs++;
								affectedSNVs.put(s, true);							}
						}
					}
					if(s.name == m.name) {
						contains = true;
						break;
					}
				}
				if(!contains) {
					int count_ref = ref_haplotype_counts.get(s).intValue();
					ref_haplotype_counts.put(s, count_ref + cellPopulationCounts.get(c)*(2 + numCNVs));
				}
			}
			// add normal contribution
			int count_ref = ref_haplotype_counts.get(s).intValue();
			ref_haplotype_counts.put(s, count_ref + numNormalCells*(2));
		}
		HashMap<Mutation.SNV, Double> freq = new HashMap<Mutation.SNV, Double>();
		for(Mutation.SNV snv : ref_haplotype_counts.keySet()) {
			freq.put(snv, (double)var_haplotype_counts.get(snv)/ref_haplotype_counts.get(snv));
		}
		numCNVAffectedSNVs = affectedSNVs.keySet().size();
		return freq;
	}
	
	public String getCompositionString() {
		String s = "";
		for(CellPopulation c : cellPopulationCounts.keySet()) {
			s += c.getName() + ": " + cellPopulationCounts.get(c) + " (";
			for(Mutation m : c.getMutations()) {
				s += m.name + " ";
			}
			s += ")\n";
		}
		return s;
	}
}
