# piRNA-seq_pipeline
A pipeline for processing sRNA-seq data tailored for piRNA discovery and annotation.

Full version of the pipeline will execute all processing steps and can be used for *de novo* detection of piRNA clusters. 
Minimal version of the pipeline will skip the piRNA cluster detection and instead summarise reads mapped to coordinates of piRNA clusters 
provided in a `gtf` annotation file.

The starting point are `fastq` files, one per sample, merged from different lanes if applicable.

<!-- 
***OBS!!!*** This is still under active developlent and testing and may not really work yet.

 -->

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

<!-- 	5.3. Adding annotation: tissue / sample, directionality, the 1st and the 10th nt signature, abundance. 

	5.4. Building a comprehensive set of piRNA clusters based on experimental evidence, with annotation from pts 1, 2, 3 (a database) (may be used for a web tool?)

	5.5. Script to parse the piRNA cluster collection from p 5.4 to obtain a `gtf` file with clusters of interest to be used directly for read summarisation in p. 6 -->

6. 	Summarisiation of the reads mapped to provied annotation files: all Ensembl `exon` features (reporting by `gene_biotype` for QC purposes); piRNA clusters provided by an already processed annotation `gtf` file or obtained from step 4; merged gtf of different annotations: Ensembl `exon` features and piRNA clusters and repeats from RepeatMasker (simple repeats excluded). (count tables)

	6.1 Generating plots of post-alignment QC metrics (biotype distribution, distribution of the 1st nt in reads mapped to piRNA clusters)



Following the completion of the pipeline clusters can be further selected, see section 



Other steps may be added at a later time:
	
* [piPipes](https://github.com/bowhan/piPipes)

* [Sports](https://github.com/junchaoshi/sports1.0)

* other tools?



<!-- ### Minimal version

steps 1, 2, 3, 6 of the above
 -->

## Installation and running

The pipeline can be run locally or on a HTC system. It is self-contained in that the software requirements are installed as a python environment, which after activation provides the reproducible setting in which the analyses are run.

### Requirements

* python 2.7

* [NGS toolbox](http://www.smallrnagroup.uni-mainz.de/software.html) (currently not available as a conda module, to be implemented)

* [reallocate.pl](http://www.smallrnagroup.uni-mainz.de/software.html) (currently not available as a conda module, to be implemented)

* [proTRAC 2.4.3](http://www.smallrnagroup.uni-mainz.de/software.html) (currently not available as a conda module, to be implemented)


Python can be installed as a [miniconda](https://docs.conda.io/en/latest/miniconda.html) distribution (recommended)

`NGS toolbox`, `reallocate.pl` and `proTRAC 2.4.3` are installed in the project directory on Rackham, and can be also fetched from there. The path to all proTRAC related software is in `config.yml`.


The pipeline has been tested on MacOS (10.14) and Linux (CentOS).

## Installation

### Get the code

To install the pipeline you need to clone this repository to get the code:

```
cd /path/to/pipeline/location

git clone https://github.com/agata-sm/piRNA-seq_pipeline.git
```
### Configure conda

You need to configure conda to be able to install all necessary packages:

```
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --add channels r
```

Check if the channels have been added correctly:

```
conda config --show channels
```

should result in:
```
channels:
  - conda-forge
  - bioconda
  - defaults
```

### Python environment

The first step is to create a python environment which contains all software dependencies. This is done once (and may take some time). 

Navigate to the directory which contains file `environment.yml` which lists all necessary packages. Create the environment named `piRNA_pipeline`. 

```
conda env create -n piRNA_pipeline -f environment.yml
```

After the environmnet have been created, it needs to be activated:


```
conda activate piRNA_pipeline
```

This environment needs to be activated before each pipeline run.

Sometimes the environment needs to be updated, when new package versions are available or new steps are added to the pipeline:

```
conda activate piRNA_pipeline
conda env update --file environment.yml
```



## Usage

To run the pipeline you need to make some changes in the configuration files. File `config.yml` contains all project-specific information such paths to the input `fastq` files, paths to installed software (if not yet implemented as a conda package), and paths to output directory, if different than the default (the directory where the Snakemake command given below is executed). File `cluster.yml` is used on a cluster, and currently its content is optimised for cluster Rackham at Uppmax. Finally, file `submit-snakemake.sh` is a bash script used on HTC to run the pipepline (each step of the pipeline is submitted as a job to `SLURM` workload manager)

### Configure the pipeline

File `config.yml` contains information on location of important files and directories:

* `datadir` - path to directory with `fastq.gz` files; the lanes have been merged at this point;

* `resdir` - path to results; this is relative to the directory where the pipeline is executed, to keep all logs and output in place

* `logdir` - path to logs; this is relative to the directory where the pipeline is executed, to keep all logs and output in place

* `snkpth` - path to snakefile; no need to change this if the code is executed on Rackham;

* `srcpth` - path to custom scripts in `pipeline/src/`; no need to change this if the code is executed on Rackham;

* `ENS_GENOME_FA` - path to `fasta` file with reference genome sequence (from Ensembl); no need to change this if the code is executed on Rackham;

* `ENS_GTF` - path to `gtf` file with reference genome annotation with gene models (the original from Ensembl with RepeatMasker entries added); no need to change this if the code is executed on Rackham;

* `REPEATMASKER_GTF` - path to `gtf` file with reference genome annotation with RepeatMasker entries (custom processed); no need to change this if the code is executed on Rackham;

* `NGSTB_PTH` - path to directory with `NGS toolbox`, `reallocate.pl` and `proTRAC 2.4.3`; no need to change this if the code is executed on Rackham;

In practice, all entries need to be changed to the relevant paths only before the first run. After that, only the path to `datadir` needs to be set.


### Run the pipeline

#### Local mode

I.e. when the pipeline runs on a local computer or a server without any job scheduling system. Steps will be executed consecutively, and the number of cores will be adjusted automatically based on the system configuration.


Navigate to the directory where you will run the pipeline, and where the results of this run will be saved:

```
cd /proj/piRNA_pipeline/MySample_1
```

Activate the conda environment:

```
conda activate piRNA_pipeline
```

Copy the config file to the current directory:

```
cp /proj/piRNA-seq_pipeline/pipeline/config.yml .
```

Confirm that the path to data `datadir` in file `config.yml` is correct! 

Start the pipeline:

```
snakemake \
    --snakefile /proj/piRNA-seq_pipeline/pipeline/Snakefile_piRNA_full \
    --rerun-incomplete
```

You can also run the pipeline in the background. This protects the run from accidental session interruption - for example when you connect remotely to the server and the session disconnects.

You can use several programs to achieve this, in this example we use [screen](https://linux.die.net/man/1/screen), which is usually already installed in your Linux distribtion.

First, start the program by typing:

```
screen 
```

A new terminal appears. You can start a process in it, disconnect from it, then reconnect at any time.

To start a new screen press `Ctrl-a`, then `c`. Type:

```
conda activate piRNA_pipeline

snakemake --snakefile /proj/piRNA-seq_pipeline/pipeline/Snakefile_piRNA_full --rerun-incomplete
```
then press `Ctrl-a`, then `d`.

To reconnect type `screen -r`







#### HTC mode

I.e. when the pipeline runs on a HTC server with SLURM job scheduling system (for example Uppmax). Steps will be queued and executed in parallel whenever possible, and the number of cores will be adjusted based on information in file `cluster.yml`.


Navigate to the directory where you will run the pipeline, and where the results of this run will be saved:

```
cd /proj/piRNA_pipeline/MySample_1
```

Activate the conda environment:

```
conda activate piRNA_pipeline
```

Copy the config files to the current directory:

```
cp /proj/piRNA-seq_pipeline/pipeline/config.yml .
cp /proj/piRNA-seq_pipeline/pipeline/cluster.yml .
cp /proj/piRNA-seq_pipeline/pipeline/submit-snakemake.sh .
```

Confirm that the path to data `datadir` in file `config.yml` is correct! 


**IMPORTANT!** File `submit-snakemake.sh` is used to actually start the pipeline - it is therefore necessary to control that the path to the snakefile `Snakefile_piRNA_full` is correct! In this pipeline version the path is on Rackham - therefore no need to change anything if the piepline is run on Rackham.


Start the pipeline:

```
bash submit-snakemake.sh
```

It is ***strongly*** recommended that you use `screen` whenever running the pipeline remotely.


Start `screen` by typing:

```
screen 
```

A new terminal appears. You can start a process in it, disconnect from it, then reconnect at any time.

To start a new screen press `Ctrl-a`, then `c`. Type:

```
conda activate piRNA_pipeline

bash submit-snakemake.sh
```
then press `Ctrl-a`, then `d`.

To reconnect type `screen -r`





### Comments

The best way to keep track of the analysis is to keep the directory structure consistent. The recommended way is to have one master copy of the `Snakefile`, and run-specific copies of the `config.yml` and `cluster.yml` (if using an HTC) in the directory of each run.


#### Conda on Rackham

Conda is installed in the project directory:

`/proj/snic2017-7-171/conda`

There is a chunk of code that needs to be added to your `.bashrc` if you want to use this installation:

```
#### Miniconda2 ####
# >>> conda init >>>
# !! Contents within this block are managed by 'conda init' !!
__conda_setup="$(CONDA_REPORT_ERRORS=false '/proj/snic2017-7-171/conda/bin/conda' shell.bash hook 2> /dev/null)"
if [ $? -eq 0 ]; then
    \eval "$__conda_setup"
else
    if [ -f "/proj/snic2017-7-171/conda/etc/profile.d/conda.sh" ]; then
        . "/proj/snic2017-7-171/conda/etc/profile.d/conda.sh"
        CONDA_CHANGEPS1=false conda activate base
    else
        \export PATH="/proj/snic2017-7-171/conda/bin:$PATH"
    fi
fi
unset __conda_setup
# <<< conda init <<<
```



## Pipeline Output

The structure of the output from executing the pipeline is:

```
pipeline running directory
|-- logdir
|-- out_slurm
`-- results
    |-- bwt_idx
    |-- cutadapt
    |-- fastqc
    |-- mapped
    |   `-- bowtie_v
    |       |-- counted_by_gene_biotype
    |       |   `-- plots
    |       |-- counted_by_gene_id
    |       `-- counted_by_repeat_id
    |-- multiqc
    |   `-- multiqc_data
    |-- nt_distribution
    |   `-- plots
    |-- proTRAC
    |   |-- merged
    |   |-- processed
    |   |   `-- plots
    |   |-- runs
    |   |   |-- proTRAC_<SAMPLEID>
    |   |   |   `-- proTRAC_<SAMPLEID>.collapsed.no-dust.map.weighted-10000-1000-b-1_<TIMESTAMP>
    |   `-- samples
    `-- temporary
```


*logdir* - logs from software used in the pipeline run; used for generation of `MultiQC` report

*out_slurm* - logs from `slurm` job scheduling system

***results*** - results

*bwt_idx* - index of the reference genome used for read mapping

*cutadapt* - trimmed reads

*fastqc* - results of read focused QC by `FastQC`

*mapped* - read alignments 

*bowtie_v* - read alignments performed by `bowtie` in `-v` mode

*counted_by_gene_biotype* - reads summarised to gene biotypes (`-g gene_biotype`)

*counted_by_gene_id* - reads summarised to gene ids (`-g gene_id`)

*counted_by_reteat_id* - reads summarised to repeat ids (`-g gene_id`) - in the original pipeline it uses a annotation file which contains only repeats

*multiqc* - results of `MultiQC` run

*nt_distribution* - tables and plots pertaining to distribution of the 1st and the 10th nucleotide in reads (input files from the `cutadapt` directory)

*proTRAC* - results of piRNA cluster detection by `proTRAC`

*runs* - unprocessed proTRAC output files

*samples* - processed results,  by sample

*processed* - processed results,  by sample

*merged* - all detected clusters merged (minimal overlap `1 bp`)

*temporary* - temporary files, should be empty after the run, may be deleted

*plots* - diagnostic plots with graphical summary of metrics related to the section where the directory is present



## Selection of Clusters of Interest


Please note that these scripts require `perl 5.18.4` and perl module `List::MoreUtils`. These are not included in the `conda` environment. Instructions below show how to use these scripts on `rackham` Uppmax cluster.




<!-- 
\begin{itemize}
  \item\texttt{src/clusterSel/sumClusters\_proTRAC.1.1.pl} collects cluster characteristics from the dataset by parsing the sample specific cluster information in \texttt{<SAMPLE>.proTRAC.tab} files; the final output file of this script links clusters in the "master file" (for example \texttt{proTRAC\_merged\_clusters.bed}) to clusters detected in individual samples; careful selection of the "master file" is recommended, more comments on this in the piRNA-seq pipeline repository;
  \texttt{src/clusterSel/getClusters\_gff.1.0.pl} uses the output file of \item\texttt{sumClusters\_proTRAC.1.1.pl} to select clusters with specifc properties; requires sample metadata file (details in the piRNA-seq pipeline repository).
\end{itemize}

 -->


## References

David Rosenkranz & Hans Zischler. 2012. proTRAC - a software for probabilistic piRNA cluster detection, visualization and analysis. doi: [10.1186/1471-2105-13-5](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-13-5); 