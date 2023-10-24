# file:   main_HRD.R
# authors: Samhita Pamidimarri Naga
# contact: <samhitapn.96@gmail.com>

# date created: 20.July.2023
# description: Main steps to perform all the HRD calculations

# Call the Functions file
source("/ifs/data/research/projects/Samhita/srcFinal/Functions_HRD.R")
options(scipen = 100)

# User option parser
options = list(
  make_option(c("--path"), type = "character",help = "Path to output directory"),
  make_option(c("--sampleList"), type = "character", help = "Path to file with list of sample names and pathology TC. If no pathology TC available, provide as NA"))

optionsParse = parse_args(OptionParser(option_list = options))
outDir <- optionsParse$path

sampleList = fread(optionsParse$sampleList,
                    stringsAsFactors = FALSE, header = FALSE)
names(sampleList) <- c("SAMPLE", "patTC")

# GET THE LIST OF GENOTYPED FILES
genotypedFiles <- list.files(paste0(outDir,"/Genotyping/",pattern = "*mpileup.txt",full.names=TRUE))

# GET THE FILE WITH 10MB:5MB WINDOW DEFINITIONS
windows_w2 = fread("/ifs/data/research/projects/Samhita/srcFinal/helperFiles/windows_hg38_10Mb_5Mb_excl.txt",
    header = TRUE,
    stringsAsFactors = FALSE)

# GET THE LOH REFERENCE COHORT
lohCohort <- fread("/ifs/data/research/projects/Samhita/srcFinal/helperFiles/lohRatioCohort.txt",
        stringsAsFactors = FALSE, 
        header = TRUE)

# ratiosCohort <- lapply(ratiosFiles, function(x){
#   data <- fread(x, stringsAsFactors = FALSE, header = TRUE)
# })

# FOR EACH OF THE SAMPLE IN THE PROVIDED SAMPLE LIST, PARSE THE GENOTYPED RESUTLS AND THEN CALCULATE PER WINDOW RATIOS AND SAVE EACH OF THESE RESULTS TO SEPARATE FILES
ratiosCohort <- lapply(sampleList$SAMPLE,function(x){
    print(x)
    snpdata_All = getSNPData(genotypedFiles,x)
    snpdata_All <- snpdata_All[!is.na(ID),]
    vafData <- snpdata_All[,c("ID","CHROM","START","END","REF","ALT","DP","ALT_AD","ALT_VAF"), with = FALSE]
    path <- paste0(outDir,"/lohVaf/",x,"_vaf.txt")
    writeOutput(vafData,path)

    ratios <- getRatios(vafData, windows_w2, 0.35, 0.65)
    path <- paste0(outDir,"/lohRatios/",x,"_ratios.txt")
    writeOutput(ratios,path)

    return(ratios)
})

# MIX THE RESULTS OF THE CURRENT COHORT WITH THE REFERENCE COHORT TO GET THE FINAL LOH COHORT
ratiosCohort <- rbindlist(ratiosCohort)
ratiosCohort1 <- data.table::copy(ratiosCohort)
ratiosCohort1 <- ratiosCohort1[,c("SAMPLE","Chr") := tstrsplit(ID, "_c", fixed=TRUE)][,list(SAMPLE,CHROM,i.START,i.END,snps_all,snps_sel,ratio)]
N <- length(unique(ratiosCohort1[,SAMPLE]))
ratiosCohort1 <- rbind(ratiosCohort1,lohCohort[!SAMPLE %in% unique(lohCohort[,SAMPLE])[1:N]])

# PERFORM LOH STEPS
ratiosCohort1 <- getLogCohort_ratio(ratiosCohort1)
lohCalls_bins <- getLOHColumns(ratiosCohort1)

# GET VARIOUS ICHORCNA RDATA FILES TO EXTRACT THE DESIRED RESULTS
rDataFileList_all <- list.files(paste0(outDir,"/ichorCNA/"),pattern = "*_1MB_c\\(2,3,4\\).RData",full.names=TRUE, recursive = TRUE)
rDataFileList_2 <- list.files(paste0(outDir,"/ichorCNA/"),pattern = "*_1MB_c\\(2\\).RData",full.names=TRUE, recursive = TRUE)
rDataFileList_3 <- list.files(paste0(outDir,"/ichorCNA/"),pattern = "*_1MB_c\\(3\\).RData",full.names=TRUE, recursive = TRUE)
rDataFileList_4 <- list.files(paste0(outDir,"/ichorCNA/"),pattern = "*_1MB_c\\(4\\).RData",full.names=TRUE, recursive = TRUE)

# FOR EACH OF THE SAMPLE PROVIDED, GENERATE ALL HRD RESULTS AND PLOT ALL THE SOLUTIONS
lapply(sampleList$SAMPLE, function(s){
        print(s)
        # DELETE THE RESULTS AND PLOT FILES IF ALREADY EXISTS
        unlink(paste0(outDir,"/hrdResults/",s,"_selectedResults.txt"))
        unlink(paste0(outDir,"/hrdResults/",s,"_allResults.txt"))
        unlink(paste0(outDir,"/hrdPlots/",s,"_selectedSolutions.pdf"))
        unlink(paste0(outDir,"/hrdPlots/",s,"_allSolutions.pdf"))

        patTC <- sampleList[SAMPLE == s,patTC]
        plotHRDSolutions(s,rDataFileList_all, rDataFileList_2, rDataFileList_3, rDataFileList_4, patTC = patTC, lohCalls_bins, outDir)
})
