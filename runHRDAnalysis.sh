#!/bin/bash

# RUN HRD ANALYSIS BASED ON THE DATA

module load anaconda/3
source activate /ifs/data/research/projects/ltg/software/miniconda3/envs/R.3.6.0/

cd /ifs/data/research/projects/Samhita/srcFinal/

Rscript main_HRD.R \
--path <OUTPUT_DIRECTORY>/<COHORT_NAME> \
--sampleList <SAMPLE_LIST>

source deactivate