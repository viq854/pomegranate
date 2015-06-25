POMEGRANATE: Cancer cell Lineage tree simulator
============

### About


### Program Parameters

##### TREE SIMULATION

```-t, --nTrees <arg>``` Number of trees to simulate (default: 100)
```-i, --nIter <arg>``` Number of tree growth iterations (default: 50)
```-snv, --probSNV <arg>``` Per node probablity of producing a descendant cell population with an acquired SSNV in a tree growth iteration (default: 0.15)
```-cnv, --probCNV <arg>``` Per node probablity of producing a descendant cell population with an acquired CNV in a tree growth iteration (default: 0)
```-probDeath <arg>``` true, "Probablity of a cell population death in each tree growth iteration (default: 0.06)
```-maxPopulationSize <arg>``` Max size of a cell population (default: 1000000)
```-minNodes <arg>``` Min number of undead cell population nodes in the tree, tree growth will continue beyond the defined number of iterations until this value is reached (default: 10)
		
##### SAMPLING

```-s, --nSamples <arg>``` Number of samples to collect, accepts multiple values, e.g. 5 10 15 (default: 5)
```-c, --coverage <arg>``` Simulated coverage to generate the VAFs, accepts multiple values, e.g. 500 1000 (default: 1000)
```-maxSubclones <arg>``` "Max number of subclones per sample (default: 5)
```-sampleSize <arg>``` Sample size (default: 100000)
```-e <arg>``` Sequencing error (default: 0.001)
```-nc <arg>``` Percent normal contamination in each sample (default: 0)
```-maxNC <arg>``` Max percent normal contamination per sample; if maxNC > nc, the percentage will be randomly generated from the range [nc maxNC] for each sample (default: 0)
```-localized``` Enable localized sampling (default: random sampling)
```-mixSubclone``` With localized sampling, add an additional subclone from a different subtree to each sample; by default, the sample is localized to a single disjoint subtree
		
##### INPUT/OUTPUT/VISUALIZATION

```-dir, --outputDir <arg>``` Directory where the output files should be created [required]
```-dot``` Produce DOT files for the simulated trees
```-sdot, --sampledDot``` Produce DOT files for the simulated trees with indicated samples		
```-sampleProfile``` "Output VAF/RC file includes an additional column with the binary sample profile for each SNV
		
##### OTHER

```-v, --verbose``` Verbose mode
```-h, --help``` Print usage
