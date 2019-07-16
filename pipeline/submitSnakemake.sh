#!/bin/bash

# Run snakemake
snakemake \
    --snakefile /proj/sllstore2017028/nobackup/nbis4067/scripts/cat_lanes/Snakefile_cat_lanes \
    --rerun-incomplete \
    --jobs 50 \
    --cluster-config cluster.yml \
    --cluster "sbatch \
                  -A {cluster.account} \
                  -t {cluster.time} \
                  -p {cluster.partition} \
                  -n {cluster.N} \
                  -J {cluster.jobname}"


