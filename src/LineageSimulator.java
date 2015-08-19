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

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Random;
import java.util.logging.ConsoleHandler;
import java.util.logging.Formatter;
import java.util.logging.Level;
import java.util.logging.LogRecord;
import java.util.logging.Logger;

import org.apache.commons.cli.BasicParser;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.uncommons.maths.random.BinomialGenerator;

/**
 * Lineage tree simulation and sampling pipeline
 */

public class LineageSimulator {
	protected static final Logger logger = Logger.getLogger("simulation.engine");
	
	public static void simulateLineageTrees(Args args) {
		int totalNumNodes = 0;
		// --- grow lineage trees --- //
		for(int t = 0; t < Parameters.NUM_TREES; t++) {
			// create the directory to store the results for each generated tree 
			File treeDir = new File(args.simPath + "/tree" + "_" + t);
			treeDir.mkdirs();
			// initial tree (only contains the root)
			SimulatedTree lineageTree = new SimulatedTree();
			// -- expand the tree --
			int iter = 0;
			while(iter < Parameters.NUM_ITERATIONS || /* there must be a min number of undead nodes */
					lineageTree.getNumNodes() < lineageTree.getNumDeadNodes() + Parameters.MIN_NUM_NODES + 1) { 
				if(lineageTree.getNumNodes() >= lineageTree.getNumDeadNodes() + Parameters.MAX_NUM_NODES + 1) {
					break;
				}
				lineageTree.grow();
				iter++;
			}
			writeOutputFile(treeDir.getAbsolutePath() + "/TREE_plain.txt", lineageTree.toString());
			if(args.generateDOT) {
				writeOutputFile(treeDir.getAbsolutePath() + "/TREE.dot", lineageTree.toDOT());
			}
			logger.fine("Generated tree " + t + " with " + lineageTree.getNumNodes() + " nodes.");
			totalNumNodes += lineageTree.getNumNodes();
			
			// --- sampling --- //
			for(int s = 0; s < Parameters.NUM_SAMPLES_ARRAY.length; s++) { 
				int numSamples = Parameters.NUM_SAMPLES_ARRAY[s];		
				ArrayList<TumorSample> samples = new ArrayList<TumorSample>();
				HashSet<CellPopulation> subclones = new HashSet<CellPopulation>();
				HashMap<Mutation.SNV, double[]> multiSampleFrequencies = new HashMap<Mutation.SNV, double[]>();
				
				// --- collect the samples from the tree ---
				if(Parameters.LOCALIZED_SAMPLING) {
					samples = lineageTree.getKLocalizedSamples(numSamples - 1);
				} else { // randomized
					for(int i = 1; i < numSamples; i++) {
						samples.add(lineageTree.getSample());
					}
				}
				if(args.generateSampledDOT) {
					writeOutputFile(treeDir.getAbsolutePath() + "/TREE_s" + numSamples + ".dot", lineageTree.toColoredDOT(samples));
				}
				lineageTree.resetColors();
				
				// --- populate the SNV VAFs for each sample ---
				for(int i = 1; i < numSamples; i++) { // + default normal sample 0
					TumorSample sample = samples.get(i-1);
					HashMap<Mutation.SNV, Double> freqMap = sample.getSNVFrequencies();
					for(Mutation.SNV snv : freqMap.keySet()) {
						if(multiSampleFrequencies.containsKey(snv)) {
							multiSampleFrequencies.get(snv)[i] = freqMap.get(snv);
						} else {
							multiSampleFrequencies.put(snv, new double[numSamples]);
							multiSampleFrequencies.get(snv)[i] = freqMap.get(snv);
						}
					}
					subclones.addAll(sample.cellPopulationCounts.keySet());
				}
				HashMap<Mutation.SNV, String> binaryProfiles = null;
				if(args.outputSampleProfile) {
					binaryProfiles = getBinaryProfile(multiSampleFrequencies, numSamples);
				}
				// --- store true VAFs --- 
				String VAFFileName =  treeDir.getAbsolutePath() + "/VAF_s" + numSamples + "_true.txt";
				writeVAFsToFile(VAFFileName, multiSampleFrequencies, binaryProfiles, numSamples);
				
				// --- generate VAFs with simulated coverage and sequencing error ---
				for(int c = 0; c < Parameters.COVERAGE_ARRAY.length; c++) {
					int coverage = Parameters.COVERAGE_ARRAY[c];
					VAFFileName =  treeDir.getAbsolutePath() + "/VAF_s" + numSamples + "_" + coverage + "X.txt";
					HashMap<Mutation.SNV, double[]> noisyMultiSampleFrequencies = addNoise(multiSampleFrequencies, coverage, numSamples);
					writeVAFsToFile(VAFFileName, noisyMultiSampleFrequencies, binaryProfiles, numSamples);
				}
				// --- store subclone information for evaluation ---
				String lineageFileName =  treeDir.getAbsolutePath() + "/SUBCLONES_s" + numSamples + ".txt";
				writeSubclonesToFile(lineageFileName, subclones);
			}
			if((t+1) % 1 == 0) logger.info("[PROGRESS] Simulated " + (t+1) + " trees.");
		}
		logger.info("[SUMMARY] Simulated " + Parameters.NUM_TREES + " trees. Average number of nodes / tree = " + (double) totalNumNodes/(Parameters.NUM_TREES));
	}
	
	/**
	 * Sample from a binomial with mean = true freq(f) and variance f(1-f)/coverage + sequencing noise 
	 */
	public static HashMap<Mutation.SNV, double[]> addNoise(HashMap<Mutation.SNV, double[]> multiSampleFrequencies, int coverage, int numSamples) {
		HashMap<Mutation.SNV, double[]> noisyMultiSampleFrequencies = new HashMap<Mutation.SNV, double[]>();
		for(Mutation.SNV snv : multiSampleFrequencies.keySet()) {
			noisyMultiSampleFrequencies.put(snv, new double[numSamples]);
			for(int i = 1; i < numSamples; i++) {
				int nReadsSNV = 0;
				if(multiSampleFrequencies.get(snv)[i] > 0) {
					BinomialGenerator b1 = new BinomialGenerator(coverage, multiSampleFrequencies.get(snv)[i], new Random());
					nReadsSNV = b1.nextValue();
				}
				int nReadsRef = coverage - nReadsSNV;
				// add sequencing noise
				int nSNV = 0;
				if(nReadsSNV > 0) {
					BinomialGenerator snvR = new BinomialGenerator(nReadsSNV, 1 - Parameters.SEQUENCING_ERROR, new Random());
					nSNV +=  snvR.nextValue();
				}
				BinomialGenerator flipR = new BinomialGenerator(nReadsRef, ((double) 1/3)*Parameters.SEQUENCING_ERROR, new Random());
				nSNV += flipR.nextValue();
				noisyMultiSampleFrequencies.get(snv)[i] = (double) nSNV/coverage;
			}
		}
		return noisyMultiSampleFrequencies;
	}
	
	public static HashMap<Mutation.SNV, String> getBinaryProfile(HashMap<Mutation.SNV, double[]> multiSampleFrequencies, int numSamples) {
		HashMap<Mutation.SNV, String> snvProfiles = new HashMap<Mutation.SNV, String>();
		for(Mutation.SNV snv : multiSampleFrequencies.keySet()) {
			String profile = "";
			for(int i = 0; i < numSamples; i++) {
				if(multiSampleFrequencies.get(snv)[i] == 0) {
					profile += "0";
				} else {
					profile += "1";
				}
			}
			snvProfiles.put(snv, profile);
		}
		return snvProfiles;
	}
	
	public static void writeVAFsToFile(String fileName, HashMap<Mutation.SNV, double[]> snvToVAFs, HashMap<Mutation.SNV, String> binaryProfiles, int numSamples) {
		String vafs = "";
		vafs += "#chrom\tpos\tdesc";
		if(binaryProfiles != null) {
			vafs += "\tprofile";
		}
		vafs += "\tnormal";
		for(int i = 1; i < numSamples; i++) {
			vafs += "\tsample" + i;
		}
		vafs += "\n";
		DecimalFormat df = new DecimalFormat("#.####");
		for(Mutation.SNV snv : snvToVAFs.keySet()) {
			String v = (snv.chr + 1) + "\t" + snv.position + "\t" + snv.name;
			if(binaryProfiles != null) {
				v += "\t" + binaryProfiles.get(snv);
			}
			for(int i = 0; i < numSamples; i++) {
				v += "\t" + df.format(snvToVAFs.get(snv)[i]);
			}
			v += "\n";
			vafs += v;
		}
		writeOutputFile(fileName, vafs);
	}
	
	public static void writeSubclonesToFile(String fileName, HashSet<CellPopulation> subclones) {
		String s = "";
		for(CellPopulation c : subclones) {
			boolean hasSNVs = false;
			for(Mutation m : c.getMutations()) {
				if(m instanceof Mutation.CNV) continue;
				s += "\t" + m.name;
				hasSNVs = true;
			}
			if(hasSNVs) {
				s += "\n";
			}
		}
		writeOutputFile(fileName, s);
	}
	
	public static void writeOutputFile(String fileName, String data) {
		try {
			FileWriter fw = new FileWriter(fileName);
			fw.write(data);
			fw.close();
		} catch (IOException e) {
			e.printStackTrace();
			System.err.println("Failed to write to the file: " + fileName);
			System.exit(-1);
		}
	}
	
	
	private static final String PROG_NAME = "pomegranate";
	private static final String SIMULATION_DATA_DIR = "simulation_results";

	// ---- LAUNCH ----
	public static void main(String[] args) {
		Options options = new Options(); 
		// commands
		//options.addOption("simulate", false, "Simulate lineage trees");
		//options.addOption("sample", false, "Sample from the simulated trees");
		//options.addOption("evaluate", false, "Evaluate trees");
		
		// tree simulation
		options.addOption("t", "nTrees", true, "Number of trees to simulate (default: 100)");
		options.addOption("i", "nIter", true, "Number of tree growth iterations (default: 50)");
		options.addOption("snv", "probSNV", true, "Per node probablity of generating a descendant cell population with an acquired SNV during a tree growth iteration (default: 0.15)");
		options.addOption("cnv", "probCNV", true, "Per node probablity of generating a descendant cell population with an acquired CNV during a tree growth iteration (default: 0.02)");
		options.addOption("probDeath", true, "Probablity of a cell population death in each tree growth iteration (default: 0.06)");
		options.addOption("maxPopulationSize", true, "Max size of a cell population (default: 1000000)");
		options.addOption("minNodes", true, "Minimum number of undead cell population nodes in a valid tree, tree growth will continue beyond the defined number of iterations until this value is reached (default: 10)");
		options.addOption("maxNodes", true, "Maximum number of undead cell population nodes in a tree, tree growth will stop after the iteration in which this value is reached/first surpassed (default: 1000)");
		
		// sampling
		Option samplesOption = new Option("s", "nSamples", true, "Number of samples to collect, accepts multiple values, e.g. 5 10 15 (default: 5)");
		samplesOption.setArgs(Option.UNLIMITED_VALUES);
		options.addOption(samplesOption);
		Option covOption = new Option("c", "coverage", true, "Simulated coverage to generate the VAFs, accepts multiple values, e.g. 500 1000 (default: 1000)");
		covOption.setArgs(Option.UNLIMITED_VALUES);
		options.addOption(covOption);
		options.addOption("maxSubclones", true, "Max number of subclones per sample (default: 5)");
		options.addOption("sampleSize", true, "Number of cells per sample (default: 100000)");
		options.addOption("e", true, "Sequencing error (default: 0.001)");
		options.addOption("minNC", true, "Minimum percentage of normal contamination per sample; the percentage will be randomly generated from the range [minNC maxNC] for each sample (default: 0)");
		options.addOption("maxNC", true, "Maximum percentage of normal contamination per sample; if maxNC < minNC, maxNC will be automatically set to minNC; the percentage will be randomly generated from the range [minNC maxNC] for each sample (default: 20)");
		//options.addOption("localized", false, "Enable localized sampling (default: random sampling)");
		//options.addOption("mixSubclone", false, "With localized sampling, add an additional subclone from a different subtree to each sample; by default, the sample is localized to a single disjoint subtree");
		
		// input/output/display
		options.addOption("dir", "outputDir", true, "Directory where the output files should be created [required]");
		options.addOption("dot", false, "Produce DOT files for the simulated trees");
		options.addOption("sdot", "sampledDot", false, "Produce DOT files for the simulated trees with indicated samples");		
		options.addOption("sampleProfile", false, "Output VAF file includes an additional column with the binary sample profile for each SNV");		
		
		// other
		options.addOption("v", "verbose", false, "Verbose mode");
		options.addOption("h", "help", false, "Print usage");
			
		// display order
		ArrayList<Option> optionsList = new ArrayList<Option>();
		optionsList.add(options.getOption("dir"));
		optionsList.add(options.getOption("t"));
		optionsList.add(options.getOption("i"));
		optionsList.add(options.getOption("snv"));
		optionsList.add(options.getOption("cnv"));
		optionsList.add(options.getOption("probDeath"));
		optionsList.add(options.getOption("maxPopulationSize"));
		optionsList.add(options.getOption("minNodes"));
		optionsList.add(options.getOption("maxNodes"));
		optionsList.add(options.getOption("s"));
		optionsList.add(options.getOption("c"));
		optionsList.add(options.getOption("maxSubclones"));
		optionsList.add(options.getOption("sampleSize"));
		optionsList.add(options.getOption("e"));
		optionsList.add(options.getOption("minNC"));
		optionsList.add(options.getOption("maxNC"));
		optionsList.add(options.getOption("dot"));
		optionsList.add(options.getOption("sdot"));
		optionsList.add(options.getOption("sampleProfile"));
		optionsList.add(options.getOption("v"));
		optionsList.add(options.getOption("h"));
			
		
		CommandLineParser parser = new BasicParser();
		CommandLine cmdLine = null;
		HelpFormatter hf = new HelpFormatter();
		hf.setOptionComparator(new OptionComarator<Option>(optionsList));
		try {
			cmdLine = parser.parse(options, args);
		} catch (ParseException e) {
			System.err.println(e.getMessage());
			hf.printHelp(PROG_NAME, options);
			System.exit(-1);
		}
		Args params = new Args();	
		if(cmdLine.hasOption("dir")) {
			params.simPath = cmdLine.getOptionValue("dir") + "/" + SIMULATION_DATA_DIR;
		} else {
			System.err.println("Required parameter: output directory path [-dir]");
			hf.printHelp(PROG_NAME, options);
			System.exit(-1);
		}
		if(cmdLine.hasOption("t")) {
			Parameters.NUM_TREES = Integer.parseInt(cmdLine.getOptionValue("t"));
		}
		if(cmdLine.hasOption("i")) {
			Parameters.NUM_ITERATIONS = Integer.parseInt(cmdLine.getOptionValue("i"));
		}
		if(cmdLine.hasOption("snv")) {
			Parameters.PROB_SNV = Double.parseDouble(cmdLine.getOptionValue("snv"));
		}
		if(cmdLine.hasOption("cnv")) {
			Parameters.PROB_CNV = Double.parseDouble(cmdLine.getOptionValue("cnv"));
		}
		if(cmdLine.hasOption("probDeath")) {
			Parameters.PROB_DEATH = Double.parseDouble(cmdLine.getOptionValue("probDeath"));
		}
		if(cmdLine.hasOption("maxPopulationSize")) {
			Parameters.MAX_POPULATION_SIZE = Integer.parseInt(cmdLine.getOptionValue("maxPopulationSize"));
		}
		if(cmdLine.hasOption("minNodes")) {
			Parameters.MIN_NUM_NODES = Integer.parseInt(cmdLine.getOptionValue("minNodes"));
			if (Parameters.MIN_NUM_NODES < 1) {
				System.err.println("Minimum number of nodes [-minNodes] must be at least 1");
				System.exit(-1);
			}
		}
		if(cmdLine.hasOption("maxNodes")) {
			Parameters.MAX_NUM_NODES = Integer.parseInt(cmdLine.getOptionValue("maxNodes"));
			if (Parameters.MAX_NUM_NODES < 1 || Parameters.MAX_NUM_NODES < Parameters.MIN_NUM_NODES) {
				System.err.println("Maximum number of nodes [-maxNodes] must be at least 1 and not less than [-minNodes]");
				System.exit(-1);
			}
		}
		if(cmdLine.hasOption("s")) {
			String[] samples = cmdLine.getOptionValues("s");
			Parameters.NUM_SAMPLES_ARRAY = new int[samples.length];
			for(int i = 0; i < samples.length; i++) {
				Parameters.NUM_SAMPLES_ARRAY[i] = Integer.parseInt(samples[i]);
			}
		}
		if(cmdLine.hasOption("c")) {
			String[] cov = cmdLine.getOptionValues("c");
			Parameters.COVERAGE_ARRAY = new int[cov.length];
			for(int i = 0; i < cov.length; i++) {
				Parameters.COVERAGE_ARRAY[i] = Integer.parseInt(cov[i]);
			}
		}
		if(cmdLine.hasOption("maxSubclones")) {
			Parameters.MAX_NUM_SUBCLONES = Integer.parseInt(cmdLine.getOptionValue("maxSubclones"));
		}
		if(cmdLine.hasOption("sampleSize")) {
			Parameters.NUM_CELLS_PER_SAMPLE = Integer.parseInt(cmdLine.getOptionValue("sampleSize"));
		}
		if(cmdLine.hasOption("e")) {
			Parameters.SEQUENCING_ERROR = Double.parseDouble(cmdLine.getOptionValue("e"));
		}
		if(cmdLine.hasOption("minNC")) {
			Parameters.MIN_PERCENT_NORMAL_CONTAMINATION = Double.parseDouble(cmdLine.getOptionValue("minNC"));
		}
		if(cmdLine.hasOption("maxNC")) {
			Parameters.MAX_PERCENT_NORMAL_CONTAMINATION = Double.parseDouble(cmdLine.getOptionValue("maxNC"));
		}
		if(Parameters.MAX_PERCENT_NORMAL_CONTAMINATION <  Parameters.MIN_PERCENT_NORMAL_CONTAMINATION) {
			Parameters.MAX_PERCENT_NORMAL_CONTAMINATION = Parameters.MIN_PERCENT_NORMAL_CONTAMINATION;
		}
		
		/*if(cmdLine.hasOption("localized")) {
			Parameters.LOCALIZED_SAMPLING = true;
		}
		if(cmdLine.hasOption("mixSubclone")) {
			Parameters.MIX_NBR_SUBTREE_SUBCLONE = true;
		}*/
		
		if(cmdLine.hasOption("dot")) {
			params.generateDOT = true;
		}
		if(cmdLine.hasOption("sampledDot")) {
			params.generateSampledDOT = true;
		}
		if(cmdLine.hasOption("sampleProfile")) {
			params.outputSampleProfile = true;
		}
		if(cmdLine.hasOption("h")) {
			new HelpFormatter().printHelp(" ", options);
		}
		// logger
		ConsoleHandler h = new ConsoleHandler();
		h.setFormatter(new LogFormatter());
		h.setLevel(Level.INFO);
		logger.setLevel(Level.INFO);
		if(cmdLine.hasOption("v")) {
			h.setLevel(Level.FINEST);
			logger.setLevel(Level.FINEST);
		}
		logger.addHandler(h);
		logger.setUseParentHandlers(false);
		
		// validate settings
		if(Parameters.PROB_SNV + Parameters.PROB_CNV + Parameters.PROB_DEATH > 1) {
			System.err.println("The sum of SSNV, CNV, and cell death probabilities cannot exceed 1");
			hf.printHelp(PROG_NAME, options);
			System.exit(-1);
		}
		simulateLineageTrees(params);
	}
	
	protected static class Args {
		String simPath;
		boolean generateDOT = false;
		boolean generateSampledDOT = false;
		boolean outputReadCounts = false;
		boolean outputSampleProfile = false;
		boolean verbose = false;
	}

	protected static class LogFormatter extends Formatter {
		public String format(LogRecord rec) {
			return rec.getMessage() + "\r\n";
		}
	}
	
	protected static class OptionComarator<T extends Option> implements Comparator<T> {
	    protected ArrayList<Option> orderedOptions;
	    public OptionComarator(ArrayList<Option> options) {
	    	orderedOptions = options;
	    }
	    public int compare(T o1, T o2) {
	        return orderedOptions.indexOf(o1) - orderedOptions.indexOf(o2);
	    }
	}

}
