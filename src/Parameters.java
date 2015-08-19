public class Parameters {
	
	// trees
	protected static int NUM_TREES = 100;
	protected static int NUM_ITERATIONS = 50;
	protected static int MIN_NUM_NODES = 10;
	protected static int MAX_NUM_NODES = 1000;
	protected static int MAX_POPULATION_SIZE = 1000000;
	protected static double PROB_SNV = 0.15;
	protected static double PROB_CNV = 0.02;
	protected static double PROB_DEATH = 0.06;
	protected static boolean UP_CNV_EFFECT = false;
	
	// sampling
	protected static int[] NUM_SAMPLES_ARRAY = {5};
	protected static int[] COVERAGE_ARRAY = {1000};
	protected static boolean LOCALIZED_SAMPLING = false;
	protected static int MAX_NUM_SUBCLONES = 5;
	protected static int NUM_CELLS_PER_SAMPLE = 100000;
	protected static double MAX_PERCENT_NORMAL_CONTAMINATION = 20;
	protected static double MIN_PERCENT_NORMAL_CONTAMINATION = 0;
	protected static boolean MIX_NBR_SUBTREE_SUBCLONE = true;
	protected static double SEQUENCING_ERROR = 0.001; // Q30
}
