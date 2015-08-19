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
import java.util.ArrayList;

public class CellPopulation {
	
	/** Number of cells in the population */
	protected int size;
	/** List of mutations in the cells */
	protected ArrayList<Mutation> mutations;
	/** Number of cells in the population, including all the undead descendants */
	protected int subtreeSize;

	protected boolean isDead;
	protected boolean isGermline;
	protected int id;	
	private static int counter = 0;
	/** For visualization */
	protected ArrayList<Color> sampleColors;
	
	public CellPopulation() {
		size = 0;
		mutations = new ArrayList<Mutation>();
		isDead = false;
		isGermline = false;
		id = counter;
		counter++;
		sampleColors = new ArrayList<Color>();
	}
	
	public void setSize(int populationSize) {
		size = populationSize;
	}
	
	public int getSize() {
		return size;
	}

	public void addMutation(Mutation snv) {
		mutations.add(snv);
	}
	
	public void setMutations(ArrayList<Mutation> populationMutations) {
		mutations = new ArrayList<Mutation>(populationMutations);
	}
	
	public ArrayList<Mutation> getMutations() {
		return mutations;
	}
	
	public boolean isDead() {
		return isDead;
	}
	
	public void setDead() {
		isDead = true;
	}
	
	public boolean isGermline() {
		return isGermline;
	}
	
	public void setGermline() {
		isGermline = true;
	}
	
	public Mutation getLastMutation() {
		if(mutations.size() > 0) {
			return mutations.get(mutations.size()-1);
		} else {
			return null;
		}
	}
	
	public void setSubtreeSize(int size) {
		subtreeSize = size;
	}
	
	public int getSubtreeSize() {
		return subtreeSize;
	}
	
	public boolean isCNV() {
		Mutation last = getLastMutation();
		if(last != null && (last instanceof Mutation.CNV)) {
			return true;
		}
		return false;
	}
	
	public String getName() {
		return getLastMutation().name;// + "," + size;
	}
	
	@Override
	public boolean equals(Object o) {
		if(o != null && (o instanceof CellPopulation)) {
			return ((CellPopulation) o).id == id;
		}
		return false;
	}
	
	@Override 
	public int hashCode() {
		return ((Integer)id).hashCode();
	}
}
