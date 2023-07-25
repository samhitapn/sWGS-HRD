# sWGS-HRD

Purpose: Pipeline overview for sWGS – HRD runs
(Optimised for run in in-house setting)

### Softwares Required:
1. R, version 3.6 or higher
2. HMM-Utils
3. ichorCNA

### R packages required
1. Data.table
2. reshape2
3. dplyr
4. ggplot2
5. hash
6. optparse
7. HMMcopy
8. GenomicRanges
9. GenomeInfoDb
10. heatmap3

Note : To run directly using in-house softwares, all the above softwares and r packages are installed in 
    "/ifs/data/research/projects/ltg/software/"

### Launch the CNA and Genotyping jobs - automatically
In the file /ifs/data/research/projects/ltg/software/sWGS_HRD/hrd_sWGS.job, modify the necessary parameters:

    <OUTPUT_DIRECTORY> : Provide the directory where the results need to be saved
    
    <COHORT_NAME> : Provide the name of the cohort 
    
    <BAMPATHs_FILE> : Provide the path to the file that contains the list of all BAM files to be used in the run
    
    <SAMPLE_LIST> : Provide the path to the file that contains the list of sample names and their corresponding pathological TC if available. The sample names should match the sample names as in the BAM files.

    In addition modify the output and error paths of the job file.
    
Launch the following job file after modifying the necessary input paths:

    sbatch "/ifs/data/research/projects/ltg/software/sWGS_HRD/hrd_sWGS.job"

### Launch the CNA and Genotyping jobs - manually
The above job will call the runJobs_sWGS_HRD.sh file which does the following tasks:
(Note : The tasks in this file, can also be chosen to run manually)

    TASK1 : CREATE ALL THE NECESSARY DIRECTORIES/ SUB-DIRECTORIES
    
    ${1} : <OUTPUT_DIRECTORY>
    
    ${2} : <COHORT_NAME>

    ${3} : <BAMPATHs_FILE>

    ${4} : <SAMPLE_LIST>

    (Refer Section "Launch the CNA and Genotyping jobs - automatically" for definition of these parameters)


    TASK 2 : CREATES ALL THE JOB FILES NECESASRY FOR ICHORCNA RUNS USING THE BAM PATH FILES PROVIDED

    These jobs use the hmm-utils and the ichorCNA tools to run the ichorCNA and generate results based on different parameters. These jobs additionally call the /ifs/data/research/projects/ltg/software/sWGS_HRD/ichorCNA.sh script to run ichorCNA. The main parameters to be noted are : 

    1. The initial normal content, provided by the NORM parameter in the ichorCNA.sh : "c(0.2,0.35,0.5,0.6,0.7,0.8)"
    2. The initial ploidy status, provided by PLOIDY parameter in the ichorCNA.sh, multiple ichorCNA runs will be carried out. These changes are provided by the different <PLOIDY#> calls in the ichorCNA.sh file.
    
    TASK 3: CREATES ALL THE JOB FILES NECESASRY FOR GENOTYPING THE SNPs USING THE BAM PATH FILES PROVIDED

    These jobs first use samtools1.5 to split the genome into individual chromosomes and then for each chromosome, uses bcftools1.15 to genotype at the desired locations as provided in the common SNP file. 
    
### Run the HRD-Analysis
Once all the above jobs are run and the results are in the desired folders, use the following to run the HRD-Analysis, update the values in <> as per definitions described in teh previous sections:

    source activate /ifs/data/research/projects/ltg/software/miniconda3/envs/R.3.6.0/

    cd /ifs/data/research/projects/Samhita/srcFinal/

    Rscript main_HRD.R \
    --path <OUTPUT_DIRECTORY>/<COHORT_NAME> \
    --sampleList <SAMPLE_LIST>
    

