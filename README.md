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

(Note:

a) To run directly using in-house softwares, all the above softwares are available in:

    /ifs/data/research/projects/ltg/software/
    
b) The R version and the necessary pacakges can all be accessed as part of a miniconda environment which can be accessed when required using: 

    module load anaconda/3
    source activate /ifs/data/research/projects/ltg/software/miniconda3/envs/R.3.6.0/
)

### Launch the ichorCNA and Genotyping jobs - automatically
In the file /ifs/data/research/projects/Samhita/srcFinal/hrd_sWGS.job, modify the necessary parameters:

    <OUTPUT_DIRECTORY> : Provide the directory where the results need to be saved
    <COHORT_NAME> : Provide the name of the cohort 
    <BAMPATHs_FILE> : Provide the path to the file that contains the list of all BAM files to be used in the run
    <SAMPLE_LIST> : Provide the path to the file that contains the list of sample names and their corresponding pathological TC if available. 
    The sample names should match the sample names as in the BAM files. 
    
    Example File as below (Without header, tab separated. If no Pathology TC is available, please provide as NA):

    Sample1    90
    Sample2    56
    Sample3    10
    Sample4    NA
    
    In addition modify the output and error paths of the job file.

Launch the following job file after modifying the necessary input parameters:

    sbatch "/ifs/data/research/projects/Samhita/srcFinal/hrd_sWGS.job"

### Launch the ichorCNA and Genotyping jobs - manually
The above job will call the runJobs_sWGS_CNA_Genotyping.sh file which does the following tasks:
(Note : The tasks in this file, can also be chosen to run manually)

    ${1} : <OUTPUT_DIRECTORY>
    ${2} : <COHORT_NAME>
    ${3} : <BAMPATHs_FILE>
    ${4} : <SAMPLE_LIST>
    (Refer Section "Launch the ichorCNA and Genotyping jobs - automatically" for definition of these parameters)
    
    TASK1 : CREATE ALL THE NECESSARY DIRECTORIES/ SUB-DIRECTORIES

    TASK 2 : CREATES ALL THE JOB FILES NECESASRY FOR ICHORCNA RUNS USING THE BAM PATH FILES PROVIDED
    These jobs use the hmm-utils and the ichorCNA tools to run the ichorCNA and generate results based on different parameters. These jobs additionally call the /ifs/data/research/projects/Samhita/srcFinal/ichorCNA.sh script to run ichorCNA. The main parameters to be noted are : 
    1. The initial normal content, provided by the NORM parameter in the ichorCNA.sh : "c(0.2,0.35,0.5,0.6,0.7,0.8)"
    2. The initial ploidy status, provided by PLOIDY parameter in the ichorCNA.sh, multiple ichorCNA runs will be carried out. These changes are provided by the different <PLOIDY#> calls in the ichorCNA.sh file.
    
    TASK 3: CREATES ALL THE JOB FILES NECESASRY FOR GENOTYPING THE SNPs USING THE BAM PATH FILES PROVIDED
    These jobs first use samtools1.5 to split the genome into individual chromosomes and then for each chromosome, uses bcftools1.15 to genotype at the desired locations as provided in the common SNP file. 
    
### Run the HRD-Analysis
Once all the above jobs are run and the results are in the desired folders, use the following to run the HRD-Analysis, update the values within <> as per definitions described in the previous sections:
    
    module load anaconda/3
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

    File : template_ichorCNA.job
    Lines :
        1) 26 (Path to the src folder to access the ichorCNA.sh file)
        
    File : ichorCNA.sh
    Lines :
        1) #33 (Path to the .rds own reference file required for generating panel of Normals for ichorCNA)
    
    File : main_HRD.R
    Lines : 
        1) #9 (Path to the HRD Functions R file) 
        2) #27 (Path to bins for LOH calculaltion)
        3) #32 (Path to the reference LOH cohort)
        
    File : Functions_HRD.R
    Lines :
        1) #326 (Path to the additional chrX windows for a dummy chrX in the ratio plots)

### Script file flow

    hrd_sWGS.job => runJobs_sWGS_CNA_Genotyping.sh => runHRDAnalysis.sh

### HG19 changes
The above scripts/ descriptions are optimised for hg38. In case of use of hg19, please make the following changes:

    File : runJobs_sWGS_HRD.R
    Lines :
        1) #24 (Path to template_ichorCNA.job) => /ifs/data/research/projects/Samhita/srcFinal/template_ichorCNA_hg19.job
        2) #46 (Path to template_splitChr_mpileup_hg38.job) => /ifs/data/research/projects/Samhita/srcFinal/template_splitChr_mpileup_hg19.job

### Flow and Schematic explaining the LOH, LST, and tAI calling
#### 1. Calling LOH
A. Describes the flow of LOH calling along with the necessary functions and scripts involved
![LOH_Flow_Functions](https://github.com/samhitapn/sWGS-HRD/blob/main/LOH_Flow_Functions.png)

B. Describes 3 different scenarios to illustrate how the LOH merge step is done
![LOH_Schematic_Description](https://github.com/samhitapn/sWGS-HRD/blob/main/LOH_Schematic_Description.png)

#### 2. Calling LST and tAI
A. Describes the flow of CNV calling followed by LST and tAI calculations
![LST_tAI_Flow_Functions](https://github.com/samhitapn/sWGS-HRD/blob/main/LST_tAI_Flow_Functions.png)

B. Describes the schema of how the LSTs and tAIs are called
![LST_tAI_Schematic](https://github.com/samhitapn/sWGS-HRD/blob/main/LST_tAI_Schematic.png)

### Details of the various R functions included in the file Functions_HRD.R

| Function Name  | Input Params | Description | Return Value |
| ------------- | ------------- | ------------- | ------------- |
| `getFiles`  | `fileList` | List of all the files from which to filter the neceesary files based on the sample name | Selected list of files |
| `getFiles`  | `sample` | Sample Name whose files need to be filtered | Selected list of files | 
| `getSNPData`  | `fileList` | List of all the files from which to filter the neceesary files based on the sample name | Resulting SNP data extracted from the genotyped files, with the ALT_VAF |
| `getSNPData`  | `sample` | Sample Name whose files need to be filtered | Resulting SNP data extracted from the genotyped files, with the ALT_VAF |
| `getRatios` | `data` | SNP Data | Data Table with per bin details of number of SNPs and the resulting ratios |
| `getRatios` | `x` | Bin Definitions | Data Table with per bin details of number of SNPs and the resulting ratios |
| `getRatios` | `s` | VAF lower limit between which the selected SNPs need to be counted | Data Table with per bin details of number of SNPs and the resulting ratios |
| `getRatios` | `s` | VAF lower limit between which the selected SNPs need to be counted | Data Table with per bin details of number of SNPs and the resulting ratios |
| `getLogCohort_ratio` | `ratios` | The ratios data as resulted from the getRatios function for the whole cohort, including the reference cohort | Data table with cohort wise ratio calculations |
| `getLOHColumns` | `data` | The cohort wise calculated data as resulted from getLogCohort_ratio function | Data table with cohort wise ratio calculations and the individual bins marked as LOH or Non-LOH |
| `getMergedGroups` | `data` | Ratio data marked with cohort wise metrics and the bin-wise LOH status | Data table with individual bins merged to call LOH segments |
| `getCummulativeCounts` | `data` | The merged data table with LOH segment calls | Number of HRD-LOH present in the sample |
| `getCummulativeCounts` | `limit` | Lower limit defined for the HRD-LOH definition; default = 15 | Data table with sample and number of HRD-LOH |
| `getLST` | `data` | ichorCNA bin wise output data with columns (sample, chr, start, end and copy number); excluding X chromosome | Data table with sample and number of LSTs |
| `gettAI` | `data` | ichorCNA segment wise output data with columns (sample, chr, start, end and copy number); excluding X chromosome | Data table with sample and number of tAIs |
| `getHRDScore` | `bins` | ichorCNA bin wise data | Data table with sample, number of LSTs, number of tAIs, number of HRD-LOH and total HRD-score (addition of all three) |
| `getHRDScore` | `cna` | ichorCNA segment wise data | Data table with sample, number of LSTs, number of tAIs, number of HRD-LOH and total HRD-score (addition of all three) |
| `getHRDScore` | `hrdLOH` | Data table with sample and number of HRD-LOH | Data table with sample, number of LSTs, number of tAIs, number of HRD-LOH and total HRD-score (addition of all three) |
| `getHRDScore` | `sample` | Sample name | Data table with sample, number of LSTs, number of tAIs, number of HRD-LOH and total HRD-score (addition of all three) |
| `getSelectedSol_HRDPlots` | `rData` | Respective .RData file from the ichorCNA output | Returns the CNV plot with LOH segments highlighted |
| `getSelectedSol_HRDPlots` | `sName` | Sample name | Returns the CNV plot with LOH segments highlighted |
| `getSelectedSol_HRDPlots` | `mergedLOH` | Data table with LOH segment calls | Returns the CNV plot with LOH segments highlighted |
| `getSelectedSol_HRDPlots` | `hrdLOH` | Data table with sample and number of HRD-LOH | Returns the CNV plot with LOH segments highlighted |
| `getSelectedSol_HRDPlots` | `patTC` | Pathological Tumour Content if provided; default NA | Returns the CNV plot with LOH segments highlighted |
| `getSelectedSol_HRDPlots` | `mode` | ichorCNA ploidy run mode. One of [c(2),c(3),c(4),c(2,3,4)] | Returns the CNV plot with LOH segments highlighted |
| `getSelectedSol_HRDPlots` | `resDir` | Output directory to store the tabulated results | Returns the CNV plot with LOH segments highlighted |
| `getAllSol_HRDPlots` | `rData` | Respective .RData file from the ichorCNA output | Returns the CNV plot with LOH segments highlighted |
| `getAllSol_HRDPlots` | `sName` | Sample name | Returns the CNV plot with LOH segments highlighted |
| `getAllSol_HRDPlots` | `mergedLOH` | Data table with LOH segment calls | Returns the CNV plot with LOH segments highlighted |
| `getAllSol_HRDPlots` | `hrdLOH` | Data table with sample and number of HRD-LOH | Returns the CNV plot with LOH segments highlighted |
| `getAllSol_HRDPlots` | `patTC` | Pathological Tumour Content if provided; default NA | Returns the CNV plot with LOH segments highlighted |
| `getAllSol_HRDPlots` | `resDir` | Output directory to store the tabulated results | Returns the CNV plot with LOH segments highlighted |
| `writeHRDResults` | `sName` | Sample name | Writes the output to a file (selected solution or all solutions): includes HRD and ichorCNA metrics |
| `writeHRDResults` | `results` | ichorCNA run output parameters | Writes the output to a file (selected solution or all solutions): includes HRD and ichorCNA metrics |
| `writeHRDResults` | `hrdScore` | Data table with sample, number of LSTs, number of tAIs, number of HRD-LOH and total HRD-score (addition of all three) | Writes the output to a file (selected solution or all solutions): includes HRD and ichorCNA metrics |
| `writeHRDResults` | `patTC` | Pathological Tumour Content if provided; default NA | Writes the output to a file (selected solution or all solutions): includes HRD and ichorCNA metrics |
| `writeHRDResults` | `mode` | ichorCNA ploidy run mode. One of [c(2),c(3),c(4),c(2,3,4)]; Default = NA | Writes the output to a file (selected solution or all solutions): includes HRD and ichorCNA metrics |
| `writeHRDResults` | `selected` | File mode if it contains the results of selected results of ichorCNA or all; default = "All" | Writes the output to a file (selected solution or all solutions): includes HRD and ichorCNA metrics |
| `writeHRDResults` | `outDir` | Output directory to store the tabulated results | Writes the output to a file (selected solution or all solutions): includes HRD and ichorCNA metrics |
| `plotHRDSolutions` | `sample` | Sample name | No return value; saves the plots and tabulated results at their respective locations |
| `plotHRDSolutions` | `rDataFileList_all` | List of all RData files in the cohort run at the initial ploidy mode : c(2,3,4) | No return value; saves the plots and tabulated results at their respective locations |
| `plotHRDSolutions` | `rDataFileList_2` | List of all RData files in the cohort run at the initial ploidy mode : c(2) | No return value; saves the plots and tabulated results at their respective locations |
| `plotHRDSolutions` | `rDataFileList_3` | List of all RData files in the cohort run at the initial ploidy mode : c(3) | No return value; saves the plots and tabulated results at their respective locations |
| `plotHRDSolutions` | `rDataFileList_4` | List of all RData files in the cohort run at the initial ploidy mode : c(4) | No return value; saves the plots and tabulated results at their respective locations |
| `plotHRDSolutions` | `patTC` | Pathological Tumour Content if provided; default NA | No return value; saves the plots and tabulated results at their respective locations |
| `plotHRDSolutions` | `lohCalls_bins` | LOH cohort wise bin out data table | No return value; saves the plots and tabulated results at their respective locations |
| `plotHRDSolutions` | `outDir` | Output directory to store the plots and results | No return value; saves the plots and tabulated results at their respective locations |
| `getCNVPlot_seg_bin_altSolutions` | `data` | ichorCNA CNA data/ segment data | Return the CNV plot with LOH segments highlighted  |
| `getCNVPlot_seg_bin_altSolutions` | `cnvBins` | ichorCNA bin wise data | Return the CNV plot with LOH segments highlighted |
| `getCNVPlot_seg_bin_altSolutions` | `results` | The ichorCNA parametric output | Return the CNV plot with LOH segments highlighted |
| `getCNVPlot_seg_bin_altSolutions` | `hrdScore` | Data table with sample, number of LSTs, number of tAIs, number of HRD-LOH and total HRD-score (addition of all three) | Return the CNV plot with LOH segments highlighted |
| `getCNVPlot_seg_bin_altSolutions` | `patTC` | Pathological Tumour Content if provided; default NA | Return the LOH-ratio plot with LOH segments highlighted |
| `getCNVPlot_seg_bin_altSolutions` | `lohSegData` | LOH cohort wise bin out data table |  Return the LOH-ratio plot with LOH segments highlighted |
| `getCNVPlot_seg_bin_altSolutions` | `selected` | File mode if it contains the results of selected results of ichorCNA or all |  Return the LOH-ratio plot with LOH segments highlighted |
| `getRatioPlot` | `data` |  Data Table with per bin details of number of SNPs and the resulting ratios | Return the LOH-ratio plot with LOH segments highlighted  |
| `getRatioPlot` | `extra` | Any extra information to be added to the title of the plot; default = "" |  Return the LOH-ratio plot with LOH segments highlighted |
| `getRatioPlot` | `segData` | LOH cohort wise bin out data table  |  Return the LOH-ratio plot with LOH segments highlighted |
| `getRatioPlot` | `hrdScore` | Data table with sample, number of LSTs, number of tAIs, number of HRD-LOH and total HRD-score (addition of all three) |  Return the LOH-ratio plot with LOH segments highlighted |
| `getRatioPlot` | `patTC` | Pathological Tumour Content if provided; default NA |  Return the LOH-ratio plot with LOH segments highlighted |
| `getRatioPlot` | `selected` | File mode if it contains the results of selected results of ichorCNA or all |  Return the LOH-ratio plot with LOH segments highlighted |







