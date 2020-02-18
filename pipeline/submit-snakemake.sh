#!/bin/bash

# Run snakemake
snakemake \
    --snakefile /proj/snic2019-30-8/nbis4067/piRNA-seq_pipeline/pipeline/Snakefile_piRNA_full \
    --rerun-incomplete \
    --jobs 50 \
    --cluster-config cluster.yml \
    --cluster "sbatch \
                  -A {cluster.account} \
                  -t {cluster.time} \
                  -p {cluster.partition} \
                  -n {cluster.N} \
                  -J {cluster.jobname}"


