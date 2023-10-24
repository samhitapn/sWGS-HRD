#!/bin/bash

# TASK 1: MAKE ALL THE REQUIRED DIRECTORIES
mkdir ${1}
mkdir ${1}/${2}
mkdir ${1}/${2}/ichorCNA
mkdir ${1}/${2}/Genotyping
mkdir ${1}/${2}/lohVaf
mkdir ${1}/${2}/lohRatios
mkdir ${1}/${2}/hrdResults
mkdir ${1}/${2}/hrdPlots

mkdir ${1}/${2}/jobFiles
mkdir ${1}/${2}/err
mkdir ${1}/${2}/out
mkdir ${1}/${2}/temp

# TASK 2: MAKE AND RUN SAMPLE WISE JOB FILES FOR ICHORCNA RUNS
while read p; do   
	echo ${p}; 
	S=$(basename $p);
    S=$(basename $p | sed 's:.bam::g')
	echo ${S};
    cp "/ifs/data/research/projects/Samhita/srcFinal/template_ichorCNA.job" "${1}/${2}/jobFiles/${S}_ichorCNA.job";
    sed -i "s:<OUTPUTDIRECTORY>:${1}/${2}/ichorCNA:g" "${1}/${2}/jobFiles/${S}_ichorCNA.job";
    sed -i "s:<DIROUT>:${1}/${2}/out:g" "${1}/${2}/jobFiles/${S}_ichorCNA.job";
    sed -i "s:<DIRERR>:${1}/${2}/err:g" "${1}/${2}/jobFiles/${S}_ichorCNA.job";

    sed -i "s:<SAMPLE>:${S}:g" "${1}/${2}/jobFiles/${S}_ichorCNA.job";
    sed -i "s:<BAMFILE>:${p}:g" "${1}/${2}/jobFiles/${S}_ichorCNA.job";
    sed -i "s:<PLOIDY1>:'c(2)':g" "${1}/${2}/jobFiles/${S}_ichorCNA.job";
    sed -i "s:<PLOIDY2>:'c(3)':g" "${1}/${2}/jobFiles/${S}_ichorCNA.job";
    sed -i "s:<PLOIDY3>:'c(4)':g" "${1}/${2}/jobFiles/${S}_ichorCNA.job";
    sed -i "s:<PLOIDY4>:'c(2,3,4)':g" "${1}/${2}/jobFiles/${S}_ichorCNA.job";
    #sbatch "${1}/${2}/jobFiles/${S}_ichorCNA.job";
done < ${3}


# TASK 3: MAKE AND RUN SAMPLE WISE JOB FILES FOR GENOTYPING RUNS
while read p; do   
	echo ${p}; 
	S=$(basename $p);
    S=$(basename $p | sed 's:.bam::g')
	echo ${S};
    for c in {1..22}; do
        cp "/ifs/data/research/projects/Samhita/srcFinal/template_splitChr_mpileup_hg38.job" "${1}/${2}/jobFiles/splitChr_mpileup_${S}_${c}.job";
	sed -i "s:<DIROUT>:${1}/${2}/out:g" "${1}/${2}/jobFiles/splitChr_mpileup_${S}_${c}.job";
    	sed -i "s:<DIRERR>:${1}/${2}/err:g" "${1}/${2}/jobFiles/splitChr_mpileup_${S}_${c}.job";
        sed -i "s:<SAMPLE>:${S}:g" "${1}/${2}/jobFiles/splitChr_mpileup_${S}_${c}.job";
        sed -i "s:<BAM>:${p}:g" "${1}/${2}/jobFiles/splitChr_mpileup_${S}_${c}.job";
        sed -i "s:<chr>:chr${c}:g" "${1}/${2}/jobFiles/splitChr_mpileup_${S}_${c}.job";
        sed -i "s:<PATH>:${1}/${2}/temp/${S}_${c}:g" "${1}/${2}/jobFiles/splitChr_mpileup_${S}_${c}.job";
        #sbatch "${1}/${2}/jobFiles/splitChr_mpileup_${S}_${c}.job";
    done
done < ${3}

# RUN HRD ANALYSIS BASED ON THE DATA
# source activate /ifs/data/research/projects/ltg/software/miniconda3/envs/R.3.6.0/

# cd /ifs/data/research/projects/Samhita/srcFinal/

# Rscript main_HRD.R \
# --path ${1}/${2} \
# --sampleList ${4}

# REMOVE THE TEMP DIRECTORIES AFTER ALL ANALYSIS
# rm ${1}/${2}/jobFiles
# rm ${1}/${2}/err
# rm ${1}/${2}/out
# rm ${1}/${2}/temp
