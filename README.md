POMEGRANATE: Cancer Lineage Tree Simulator
============

### About

POMEGRANATE simulates tumor progression from normal tissue producing a branching hierarchy of monoclonal cell populations in accordance with the branched-tree cancer evolution model. Starting with the GL cell population, the simulator iteratively expands the lineage tree by introducing (with some given probability) new daughter cell populations corresponding to newly acquired somatic SNV or CNV events. In particular, in every iteration, each cell population present in the tree can give rise to a new population of cells (with a randomly generated size) representing a new SNV or CNV event, as well as undergo cell death. The user can specify the number of tree growth iterations, as well the probabilities for each event. Each simulated SNV is randomly associated with a genome location (chromosome and position) and haplotype; the CNVs are associated with a chromosome arm and haplotype and currently correspond to a duplication of this chromosome arm. This process results in lineage trees with an arbitrary number of branches and nodes. 

Multiple samples can then be collected from the produced lineage tree(s). Each sample consists of several cell populations (nodes) of the tree, where each such cell population represents a subclone in the sample. The randomized sampling process selects a random subset of nodes from the tree for each sample. Given a selected subset of cell populations, the sample is then created by obtaining a fraction of the cells from each population by sampling from a multinomial distribution with probabilities corresponding to the cell population sizes. 

The program outputs the produced trees, sampled subclones, and the per sample VAFs of each SNV for specified coverage(s), sequencing accuracy, and normal contamination.

Example of a simulated tree with 9 collected tumor samples indicated by different colors; each node represents a cell population that contains every mutation in its lineage; gray indicates dead cell populations.: 

![tree]( https://github.com/viq854/pomegranate/blob/master/data/examples/t3_s10.png "Simulated tree example")

Example of a simulated tree with 14 tumor samples: 

![tree]( https://github.com/viq854/pomegranate/blob/master/data/examples/t1_s15.png "Simulated tree example")

### Program Parameters

##### TREE SIMULATION

```-t, --nTrees <arg>``` Number of trees to simulate (default: 100)  
```-i, --nIter <arg>``` Number of tree growth iterations (default: 50)  
```-snv, --probSNV <arg>``` Per node probablity of generating a descendant cell population with an acquired SNV during a tree growth iteration (default: 0.15)  
```-cnv, --probCNV <arg>   ``` Per node probablity of generating a descendant cell population with an acquired CNV during a tree growth iteration (default: 0.02)  
```-probDeath <arg>        ``` true, "Probablity of a cell population death in each tree growth iteration; dead cell populations will not be expanded in subsequent iterations and will not be sampled from (default: 0.06)  
```-maxPopulationSize <arg>``` Maximum size of a cell population (default: 1000000)  
```-minNodes <arg>         ``` Minimum number of undead cell population nodes in a valid tree, tree growth will continue beyond the defined number of iterations until this value is reached (default: 10)  
```-maxNodes <arg>         ``` Maximum number of undead cell population nodes in a tree, tree growth will stop after the iteration in which this value is reached/first surpassed (default: 1000)
		
##### SAMPLING

```-s, --nSamples <arg>``` Number of samples to collect, accepts multiple values, e.g. 5 10 15 (default: 5)  
```-c, --coverage <arg>``` Simulated coverage to generate the VAFs, accepts multiple values, e.g. 500 1000 (default: 1000)  
```-maxSubclones <arg>``` "Max number of subclones per sample; the acutal number of subclones is randomly drawn from 1 to maxSubclones (default: 5)  
```-sampleSize <arg>``` Sample size (default: 100000)  
```-e <arg>``` Sequencing error (default: 0.001)  
```-minNC <arg>``` Minimum percentage of normal contamination per sample; the percentage will be randomly generated from the range [minNC maxNC] for each sample (default: 0) (default: 0)  
```-maxNC <arg>``` Maximum percentage of normal contamination per sample; if maxNC < minNC, maxNC will be automatically set to minNC; the percentage will be randomly generated from the range [minNC maxNC] for each sample (default: 0) 
		
##### INPUT/OUTPUT/VISUALIZATION  

```-dir, --outputDir <arg>``` Directory where the output files should be created [required]  
```-dot``` Produce DOT files for the simulated trees  
```-sdot, --sampledDot``` Produce DOT files for the simulated trees with indicated samples  
```-sampleProfile``` Include an additional column with the binary sample profile for each SNV in the VAF output file
		
##### OTHER

```-v, --verbose``` Verbose mode  
```-h, --help``` Print usage  

### How to Run

From the ```release/``` directory:

```
./pomegranate -dir <output> [options]
```

### Examples

(1) Simulate 50 trees; extract 10 and 15 samples; generate VAFs for coverages of 1,000X and 10,000X; produce the dot files for visualization; output files will be created inside the output_dir directory:
```
./pomegranate -dir output_dir -t 50 -s 10 15 -c 1000 10000 -dot -sdot -v
```

(2) Simulate the default 100 trees; extract 5, 10, and 15 samples; set the normal contamination in each sample to a random value between 10% and 40%:
```
./pomegranate -dir output_dir -s 5 10 15 -minNC 10 -maxNC 40 -dot -sdot -v
```

### Output 

A directory called ```simulation_results/``` will be created inside the output directory specified by the user. Data files associated with each simulated tree will be stored inside in separate directories. For instance, the results for example (1) above are shown below. Here the top simulation directory contains 50 subdirectories for each simulated tree. 

```
simulation_results/
├── tree_0
│   ├── TREE_plain.txt
│   ├── TREE.dot
│   ├── TREE_s10.dot
│   ├── TREE_s15.dot
│   ├── VAF_s10_10000X.txt
│   ├── VAF_s10_1000X.txt
│   ├── VAF_s10_true.txt
│   ├── VAF_s15_10000X.txt
│   ├── VAF_s15_1000X.txt
│   └── VAF_s15_true.txt
│   ├── SUBCLONES_s10.txt
│   ├── SUBCLONES_s15.txt
├── tree_1
├── tree_2
├── tree_3
├── tree_4
├── tree_5
├── tree_6
├── tree_7
├── tree_8
└── tree_9
├── tree_10
├── tree_11
├── tree_12
...
├── tree_44
├── tree_45
├── tree_46
├── tree_47
├── tree_48
├── tree_49
```

Several files are produced for each tree. The ```TREE_plain.txt``` file stores the topology of the tree (as a list of edges and nodes) and the information about each simulated mutation. The ```.dot``` files store the tree in DOT format that can be visualized using Graphviz. The ```TREE_s10.dot``` and ```TREE_s15.dot``` display which cell populations were sampled for the 10 and 15 samples experiments, respectively. To visualize these files, the following command can be used (this will produce a pdf file with the same name):

```dot -Tpdf simulation_results/tree_0/TREE_s10.dot -O```  

The ```VAF```-prefixed files contain the multi-sample variant allele frequency (VAF) information for each sampled SNV. Separate files are generated for the different coverages (e.g. 1,000X and 10,000X). The ```VAF_..._true.txt``` contains the actual VAF of each SNV across the samples. For example, this is a fragment of a file containing SNVs from 15 samples at 1000X coverage:

```
#chrom  pos     desc    normal  sample1 sample2 sample3 sample4 sample5 sample6 sample7 sample8 sample9
2       54659465        M114    0       0.0003  0.0005  0.0001  0.0004  0.0632  0.0003  0.0002  0.0001  0.0004
11      106272044       M0      0       0.0708  0.2293  0.2298  0.3559  0.1852  0.1897  0.2192  0.1927  0.3789
22      19631012        M15     0       0.0758  0.0264  0.0561  0.1175  0.0743  0.1923  0.0319  0.007   0.2515
6       22699228        M50     0       0.0006  0.0003  0.0004  0.0846  0.0283  0.0003  0.0004  0.0003  0.1314
12      72038217        M64     0       0.0005  0       0.0002  0.1089  0.0004  0.0002  0.0008  0.0006  0.0003
17      67681707        M9      0       0.0717  0.0254  0.0531  0.1106  0.0751  0.2013  0.0437  0.0651  0.2391
4       180672416       M11     0       0.0007  0.2094  0.0003  0.0004  0.0002  0.0002  0.0002  0.0003  0.0003
7       26582817        M54     0       0.165   0.0005  0.0007  0.0003  0.0003  0.0004  0.0001  0.0001  0
21      33290675        M41     0       0.0003  0.0001  0.1186  0.0001  0       0.0001  0.0749  0.049   0.0004
14      86616439        M151    0       0.0001  0.0002  0.0589  0.0004  0.0002  0.0001  0.0006  0.0002  0.0002
1       128266807       M37     0       0.1723  0.0865  0.0221  0.0003  0.1642  0       0.0006  0.0355  0.0003
13      19210482        M21     0       0.0004  0.0002  0.0003  0.0003  0.0003  0.0003  0.0007  0.0157  0.0003
11      17572309        M63     0       0.0005  0.0251  0.0002  0.0002  0.0003  0.0004  0.0004  0.0003  0.0003
4       122102509       M66     0       0.0008  0.0002  0.0004  0.0002  0.0003  0.0005  0.017   0.0002  0.0001
9       70202460        M10     0       0.0005  0.0002  0.1204  0.0005  0.0002  0.0002  0.0986  0.0696  0.0001
18      33116420        M60     0       0.0705  0.0003  0.0008  0.0008  0.0113  0.0002  0.0003  0.0076  0.0001
3       69139413        M161    0       0.0004  0.0003  0.0078  0.0001  0.0006  0.0003  0.0004  0.0006  0.0005
```

Finally the ```SUBCLONES```-prefixed files contain information about the composition of each cell population sampled from the tree. 

### System Requirements

Java Runtime Environment (JRE) 1.6  
(Optional) Graphviz: for output tree and sampling visualization 

###License

MIT License 

###Support

For help running the program or any questions/suggestions/bug reports, please contact Victoria Popic at viq@stanford.edu
