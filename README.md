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

Note: 
a) To run directly using in-house softwares, all the above softwares are available in:

    /ifs/data/research/projects/ltg/software/
    
b) The R version and the necessary pacakges can all be accessed as part of a miniconda environment which can be accessed when required using: 

    source activate /ifs/data/research/projects/ltg/software/miniconda3/envs/R.3.6.0/

### Launch the CNA and Genotyping jobs - automatically
In the file /ifs/data/research/projects/Samhita/srcFinal/hrd_sWGS.job, modify the necessary parameters:

    <OUTPUT_DIRECTORY> : Provide the directory where the results need to be saved
    <COHORT_NAME> : Provide the name of the cohort 
    <BAMPATHs_FILE> : Provide the path to the file that contains the list of all BAM files to be used in the run
    <SAMPLE_LIST> : Provide the path to the file that contains the list of sample names and their corresponding pathological TC if available. The sample names should match the sample names as in the BAM files.
    In addition modify the output and error paths of the job file.

Launch the following job file after modifying the necessary input paths:

    sbatch "/ifs/data/research/projects/Samhita/srcFinal/hrd_sWGS.job"

### Launch the CNA and Genotyping jobs - manually
The above job will call the runJobs_sWGS_CNA_Genotyping.sh file which does the following tasks:
(Note : The tasks in this file, can also be chosen to run manually)

    TASK1 : CREATE ALL THE NECESSARY DIRECTORIES/ SUB-DIRECTORIES
    ${1} : <OUTPUT_DIRECTORY>
    ${2} : <COHORT_NAME>
    ${3} : <BAMPATHs_FILE>
    ${4} : <SAMPLE_LIST>
    (Refer Section "Launch the CNA and Genotyping jobs - automatically" for definition of these parameters)

    TASK 2 : CREATES ALL THE JOB FILES NECESASRY FOR ICHORCNA RUNS USING THE BAM PATH FILES PROVIDED
    These jobs use the hmm-utils and the ichorCNA tools to run the ichorCNA and generate results based on different parameters. These jobs additionally call the /ifs/data/research/projects/Samhita/srcFinal/ichorCNA.sh script to run ichorCNA. The main parameters to be noted are : 
    1. The initial normal content, provided by the NORM parameter in the ichorCNA.sh : "c(0.2,0.35,0.5,0.6,0.7,0.8)"
    2. The initial ploidy status, provided by PLOIDY parameter in the ichorCNA.sh, multiple ichorCNA runs will be carried out. These changes are provided by the different <PLOIDY#> calls in the ichorCNA.sh file.
    
    TASK 3: CREATES ALL THE JOB FILES NECESASRY FOR GENOTYPING THE SNPs USING THE BAM PATH FILES PROVIDED
    These jobs first use samtools1.5 to split the genome into individual chromosomes and then for each chromosome, uses bcftools1.15 to genotype at the desired locations as provided in the common SNP file. 
    
### Run the HRD-Analysis
Once all the above jobs are run and the results are in the desired folders, use the following to run the HRD-Analysis, update the values within <> as per definitions described in the previous sections:

    source activate /ifs/data/research/projects/ltg/software/miniconda3/envs/R.3.6.0/

    cd /ifs/data/research/projects/Samhita/srcFinal/
    Rscript main_HRD.R \
    --path <OUTPUT_DIRECTORY>/<COHORT_NAME> \
    --sampleList <SAMPLE_LIST>

    source deactivate

The above code is also in runHRDAnalysis.sh

### Explanation of the Directories structure

    1. ichorCNA : All the results from all the ichorCNA runs
    
    2. Genotyping : Stores the results generated from bcftools parsing; there will be individual .txt files per sample per chromosome

    3. lohVAF : Stores the processed genotype results with VAF calculated, that will be used further; there will be individual .txt file per sample

    4. lohRatios : Calculated ratios per sample per 10Mb bin with an overlap of 5Mb is created; there will be individual .txt file per sample

    5. hrdResults : All the calcualted HRD metrics : LST, LOH, tAI based on each ichorCNA ploidy run along with ichorCNA reported ploidy and purity at each run and additional parameters are reported. There will be 3 files per sample, one named <SAMPLE>_selectedResults.txt which will have the above output for only the ichorCNA selected result - per ploidy mode run; the other named <SAMPLE>_allResults.txt will have all the above output for all different combination of ichorCAN runs.

    6. hrpPlots : There will be again 2 files per sample. One named <SAMPLE>_selectedSolutions.pdf will have 4 genome wide CNA plots along with the LOH ratio plot. One named <SAMPLE>_allSOlutions.pdf will have 18 differnt genome wide CNA plots and one ratio plot (each plot per page) summarising/ displaying all the HRD results.

    7. jobFiles : Directory to store all the job files created for ichorCNA and genotyping runs

    8. err and out : directories to store .err and .out files for each of the runs

    9. temp : To temporarily store the per chromosome bam and bai files. This directory should be empty at the end of the genotyping.

(7,8 and 9 directories can be deleted once the run is complete)

### Hard coded paths that may need to be changed manually

    File : runJobs_sWGS_HRD.R
    Lines :
        1) #24 (Path to template_ichorCNA.job)
        2) #46 (Path to template_splitChr_mpileup_hg38.job)
       
    File : main_HRD.R
    Lines : 
        1) #9 (Path to the HRD Functions R file) 
        2) #27 (Path to bins for LOH calculaltion)
        3) #32 (Path to the reference LOH cohort)
        
    File : Functions_HRD.R
    Lines :
        1) #326 (Path to the additional chrS windows for a dummy chrX in the ratio plots)

    File : ichorCNA.sh
    Lines :
        1) #33 (Path to the .rds own reference file required for generating panel of Normals for ichorCNA)

### Script file flow

    hrd_sWGS.job => runJobs_sWGS_CNA_Genotyping.sh => runHRDAnalysis.sh
