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

import java.util.Random;

public class Mutation {
	public static final int NUM_CHROMOSOMES = 23; 
	public static final int[] CHROMOSOME_LENGTHS = {249250621, 243199373, 198022430, 191154276, 
													180915260, 171115067, 159138663, 155270560, 
													146364022, 141213431, 135534747, 135006516, 
													133851895, 115169878, 107349540, 102531392, 
													90354753, 81195210, 78077248, 63025520,
													59373566, 59128983, 51304566, 48129895};
	protected String name; // unique 
	protected int chr;
	protected int haplotype;
	private static int counter = 0;
	protected Random r = new Random();
	
	public Mutation() {
		name = "M" + counter;
		chr = r.nextInt(NUM_CHROMOSOMES);
		haplotype = r.nextInt(2);
		counter++;
	}
	
	public static class SNV extends Mutation {
		protected int position;
		public SNV() {
			super();
			position = r.nextInt(CHROMOSOME_LENGTHS[chr]);
		}
		
		public SNV(CNV parent) {
			super();
			chr = parent.chr;
			position = r.nextInt(CHROMOSOME_LENGTHS[chr]/2);
			if(parent.arm == 1) {
				position += CHROMOSOME_LENGTHS[chr]/2;
			}
		}
		
		public String toString() {
			return name + ": chr=" + (chr + 1) + ", pos="  + position + ", haplotype=" + haplotype;
		}
	}
	
	public static class CNV extends Mutation {
		protected int arm;
		public CNV() {
			super();
			arm = r.nextInt(2);
			name = "CNV_" + name;
		}
		
		public CNV(SNV parent) {
			super();
			chr = parent.chr;
			if(parent.position <= CHROMOSOME_LENGTHS[chr]/2) {
				arm = 0;
			} else {
				arm = 1;
			}
			name = "CNV_" + name;
		}
		
		public String toString() {
			return name + ": chr=" + (chr + 1) + ", arm="  + arm + ", haplotype=" + haplotype;
		}
	}
}