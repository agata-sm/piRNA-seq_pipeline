# piRNA-seq_pipeline
A pipeline for processing sRNA-seq data tailored for piRNA discovery and annotation.

Full version of the pipeline will execute all processing steps and can be used for *de novo* detection of piRNA clusters. 
Minimal version of the pipeline will skip the piRNA cluster detection and instead summarise reads mapped to coordinates of piRNA clusters 
provided in a `gtf` annotation file.

The starting point are `fastq` files, one per sample, merged from different lanes if applicable.

## Workflow description

The output of each step is indicated in parentheses.

### Full version

1. Trimming of the 3' adapter and removal of very short sequences (<17 nt)

2. QC of the processed reads

	1.1. `FastQC` for general read metrics (`html` file for each sample)

	1.2. Distribution of the 1st and the 10th nucleotide (tables, plots for each sample)

3. Read mapping using `bowtie -v 1 -a --best --strata` ; sorting by position and indexing (`bam` files)
	
	(the pipeline will also generate an index for bowtie mapping, if not present already)

4. piRNA cluster detection using `proTRAC` including all preprocessing steps

5. Post-processing of proTRAC output
	
	5.1. Filtering of clusters based on the 1st nt signature and read length distribution

	5.2. Merging overlapping clusters from different samples

	5.3. Adding annotation: tissue / sample, directionality, the 1st and the 10th nt signature, abundance. 

	5.4. Building a comprehensive set of piRNA clusters based on experimental evidence, with annotation from pts 1, 2, 3 (a database) (may be used for a web tool?)

	5.5. Script to parse the piRNA cluster collection from p 5.4 to obtain a `gff` file with clusters of interest to be used directly for read summarisation in p. 6

6. 	Summarisiation of the reads mapped to provied annotation files: all Ensembl `exon` features (reporting by `gene_biotype` for QC purposes); piRNA clusters provided by an already processed annotation `gtf` file or obtained from step 4; merged gtf of different annotations: Ensembl `exon` features and piRNA clusters and repeats from RepeatMasker (simple repeats excluded). (count tables)

	6.1 Generating plots of post-alignment QC metrics (biotype distribution, distribution of the 1st nt in reads mapped to piRNA clusters)



Other steps may be added at a later time:
	
* piPipes https://github.com/bowhan/piPipes

* Sports https://github.com/junchaoshi/sports1.0

* other tools?



### Minimal version

steps 1, 2, 3, 6 of the above

## Installation and running

The pipeline can be run locally or on a HTC system. It is self-contained in that the software requirements are installed as a python environment, which after activation provides the reproducible setting in which the analyses are run.

### Requirements

* python 3.x

* NGS toolbox from **XXX** (currently not available as a conda module, to be implemented)

Python can be installed as a miniconda distribution (recommended), see here **XXXX**

The pipeline has been tested on MacOS (10.14) and Linux (CentOS **XXXX**).

## Installation

To install the pipeline you need to clone the repo:

```
cd /path/to/pipeline/location

git clone https://github.com/agata-sm/piRNA-seq_pipeline.git
```

The first step is to create the python environment which contains all software dependencies. This is done once, and the environment needs to be activated before each run.

```
```

The best way to keep track of the analysis is to keep the directory structure consistent. The recommended way is to have one master copy of the `Snakefile`, and run-specific copies of the `config.yml` and `cluster.yml` (if using an HTC) in the directory of each run.



## Running