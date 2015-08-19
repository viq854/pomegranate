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

From the release/ directory:

```
./pomegranate -dir <output>
```

### Examples

Simulate 50 trees; extract 10 and 15 samples; generate VAFs for coverages of 1,000X and 10,000X; produce the dot files for visualization:
```
./pomegranate -dir output_dir -t 50 -s 10 15 -c 1000 10000 -dot -sdot -v
```

Simulate the default 100 trees; extract 5, 10, and 15 samples; set the normal contamination in each sample to a random value between 10% and 40%:
```
./pomegranate -dir output_dir -s 5 10 15 -minNC 10 -maxNC 40 -dot -sdot -v
```

### Output 

### System Requirements

Java Runtime Environment (JRE) 1.6  
(Optional) Graphviz: for output tree and sampling visualization 

###License

MIT License 

###Support

For help running the program or any questions/suggestions/bug reports, please contact viq@stanford.edu
