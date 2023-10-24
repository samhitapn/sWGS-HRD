# file:   Functions_HRD.R
# authors: Samhita Pamidimarri Naga
# contact: <samhitapn.96@gmail.com>

# date created: 20.July.2023
# description: Definition of all the functions necessary to calculate, plot and save the results for HRD phenotype using the CNA and LOH results from sWGS


library(data.table)
library(reshape2)
library(dplyr)
library(ggplot2)
library(hash)
library(optparse)

# GET FILES
getFiles = function(fileList, sample){
    filelist = unlist(fileList)
    term = paste0(sample,"_")
    selects = grep(term, fileList)
    files = fileList[selects]
    return(files)
}

############## LOH ANALYSIS ###############
# GET SNP DATA FROM GENOTYPED RESULTS
getSNPData = function(fileList, sample){
    files = getFiles(fileList, sample)
    print(files)
    data <- lapply(files,function(f){
        print(f)
        s <- basename(f)
        s <- gsub(".txt","",s)
   
        snps <- fread(f,
                header = FALSE,stringsAsFactors = FALSE)
        setnames(snps,c("chr","pos","REF","ALT","DP","AD"))

            snps[DP != 0,ALT := gsub(",<*>","",ALT,fixed = TRUE)]
            snps[,ALT := ifelse(ALT == "<*>", REF,ALT)]
            snps[nchar(REF) == 1 & nchar(ALT) == 1,]
            snps[,c("REF_AD","ALT_AD") := tstrsplit(AD, ",", fixed=TRUE,keep=c(1,2))]
            snps[,`:=` (DP = as.numeric(DP), REF_AD = as.numeric(REF_AD), ALT_AD = as.numeric(ALT_AD))]
            snps[, `:=` (ALT_VAF = ALT_AD/DP, CHROM = chr, chr = gsub("chr","",chr), START = pos, END = pos, ID = s)]
            snps[,c("SAMPLE","Chr","ExtraName","Ex1") := tstrsplit(s, "_", fixed=TRUE)]
            snps <- snps[,c("ID","SAMPLE","CHROM","START","END","REF","ALT","DP","REF_AD","ALT_AD","ALT_VAF")]
        
        return(snps)
    })
    snpData = rbindlist(data,fill = TRUE)
    return(snpData)
}

# CALCUALTE RATIOS FOR LOH PER SAMPLE
getRatios = function(data, x, s, e){
    
    y <- data
    setkey(y,CHROM,START,END)
    allCount <- foverlaps(x,y, type = "any")
    allCount[,snps_all := .N, by = .(CHROM,i.START, i.END)]
    y1 <- y[between(as.numeric(ALT_VAF),s,e),]
    setkey(y1,CHROM,START,END)
    selCount <- foverlaps(x,y1,type = "any")
    selCount[,snps_sel := .N, by = .(CHROM,i.START, i.END)]
    snpResult <- merge(allCount,selCount)
    #snpResult[!is.na(ID),]
    snpResult[,ratio := as.numeric(snps_sel)/as.numeric(snps_all)]
    snpResult <- unique(snpResult[,c("ID","CHROM","i.START","i.END","snps_all","snps_sel","ratio")])
    snpResult <- snpResult[!is.na(ID),]
   
    return(snpResult)
}

# CALCULATE COHORT WISE LOH METRICS
getLogCohort_ratio <- function(ratios){
  ratioData_logNorm <- data.table::copy(ratios)
  ratioData_logNorm <- ratioData_logNorm[,median := median(ratio),by = SAMPLE][,logR := log(ratio/median,2)][, `:=` (median_s = median(logR),mad1_s = mad(logR),mad2_s = 2*mad(logR))]
  ratioData_logNorm <- ratioData_logNorm[,m_mad1 := ifelse((median_s - mad1_s)<=logR,"In_mad1","Out_mad1")]
  ratioData_logNorm <- ratioData_logNorm[,m_mad2 := ifelse((median_s - mad2_s)<=logR,"In_mad2","Out_mad2")]
  ratioData_logNorm <- ratioData_logNorm[,madRange := ifelse(m_mad2 == "In_mad2" & m_mad1 == "In_mad1","InRange",
                                                             ifelse(m_mad2 == "Out_mad2" & m_mad1 == "Out_mad1","Out_mad2",
                                                               ifelse(m_mad1 == "Out_mad1" & m_mad2 == "In_mad2","Out_mad1","Others")))]
  return(ratioData_logNorm)      
}

# MARK THE LOH PER BIN BASED ON COHORT METRICS
getLOHColumns <- function(data){
  data <- data[, `:=` (loh_MAD1 = ifelse(logR <= median_s - mad1_s,"LOH","Non-LOH"),
                       loh_MAD2 = ifelse(logR <= median_s - mad2_s,"LOH","Non-LOH"))]
  data <- data[,list(SAMPLE,CHROM,i.START,i.END,logR,median_s,mad1_s,mad2_s,loh_MAD1,loh_MAD2,madRange), which = FALSE]
  return(data)
}

# GET THE GROUP NAMES IN ORDER, ALSO USED FOR LST CALCULATIONS
getGroups = function(data){
  rle = rle(data$gp)
  rle1 = data.table(ca = rle$values,lengths = rle$lengths,rowNum = cumsum(rle$lengths) - rle$lengths + 1)
  groups = rep(paste0(rle1$ca,"_",rle1$rowNum),rle1$lengths)
  return(groups)
}

# GET THE MERGED GROUPS BASED ON LOH BINS
getMergedGroups <- function(data){
  manualSeg <- data[,gp := paste0(SAMPLE,"_",CHROM,"_",loh_MAD1)]
  groups <- getGroups(manualSeg)
  manualSeg1 <- cbind(manualSeg,groups)
  setDT(manualSeg1)
  
  manualSeg2 <- manualSeg1[,count := .N, by = groups]
  
  manualSeg2 <- manualSeg2[,LOHNew := ifelse(loh_MAD1 == "LOH" & count >= 3,"LOH","Non-LOH")]
  
  groups <- getGroups(manualSeg2)
  manualSeg2 <- cbind(manualSeg2,"groups1" = groups)
  manualSeg2 <- manualSeg2[,count1 := .N, by = groups1][,numRow := 1:.N]
  manualSeg2 <- manualSeg2[,LOHNew1 := ifelse(numRow %in% c(1,.N) & count1 == 1 & LOHNew == "Non-LOH","LOH",
                                              ifelse(numRow %in% c(1,.N) & count1 == 1 & LOHNew == "LOH","Non-LOH",LOHNew))]
  manualSeg2 <- manualSeg2[count1 != 1 | LOHNew1 != "Non-LOH",]
  manualSeg2 <- manualSeg2[,gp := paste0(SAMPLE,"_",CHROM,"_",LOHNew1)]
  
  
  groups <- getGroups(manualSeg2)
  manualSeg3 <- cbind(manualSeg2,"groupesNew" = groups)
  setDT(manualSeg3)
  
  manualSeg3 <- manualSeg3[,newStart := .SD[1][,i.START],by = groupesNew]
  manualSeg3 <- manualSeg3[,newEnd := .SD[.N][,i.END],by = groupesNew]
  manualSeg3 <- unique(manualSeg3[,list(SAMPLE,CHROM,newStart,newEnd,LOHNew1)])
  manualSeg3 <- manualSeg3[,check := 1:.N,by = .(SAMPLE,CHROM)]
  manualSeg3 <- manualSeg3[,start := ifelse(check != 1,as.numeric(newStart) + 5000000 + 1,as.numeric(newStart)),by = .(SAMPLE,CHROM)]
  manualSeg3 <- manualSeg3[,end := ifelse(check < .N,as.numeric(newEnd) - 5000000,as.numeric(newEnd)),by = .(SAMPLE,CHROM)]
  manualSeg3 <- manualSeg3[,length := end - start][,list(SAMPLE,CHROM,start,end,length,LOHNew1)]
  names(manualSeg3)[c(3,4,5,6)] <- c("START","END","LENGTH","LOH") 

  return(manualSeg3)
}

# COUNT THE NUMBER OF LOH PER SAMPLE
getCummulativeCounts <- function(data,limit = 15){
  data <- data[,countLOH := .N, by = .(SAMPLE,CHROM)]
  data <- data[,(paste0("LOH_at_",limit)) := ifelse(LENGTH >= limit*1000000 & LOH == "LOH" & countLOH > 1,"Yes","No"), by = .(SAMPLE,CHROM)]
  data <- data[,(paste0("HRD-LOH_",limit)) := .N, by = c("SAMPLE",(paste0("LOH_at_",limit)))]
  temp <- unique(data[get(names(data)[8]) == "Yes",list(SAMPLE,`HRD-LOH_15`)])
  return(temp)
}

############## LST ANALYSIS ###############
# GET ARM
centromere_Start = hash(c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,"X","Y"),
                        c(121700000,91800000,87800000,48200000,46100000,58500000,58100000,43200000,
                        42200000,38000000,51000000,33200000,16500000,16100000,17500000,35300000,
                        22700000,15400000,24200000,25700000,10900000,13700000,58100000,10300000))

centromere_End = hash(c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,"X","Y"),
                      c(125100000,96000000,94000000,51800000,51400000,62600000,62100000,
                      47200000,45500000,41600000,55800000,37800000,18900000,18200000,20500000,
                      38400000,27400000,21500000,28100000,30400000,13000000,17400000,63800000,10600000))

getArm_new = function(row){
  s = centromere_Start[[as.character(row[[1]])]]
  e = centromere_End[[as.character(row[[1]])]]
  arm = ifelse(as.numeric(row[[2]]) <= s & as.numeric(row[[3]]) < e, "p", "q")
  return(arm)
}

# CALCUALTE THE NUMBER OF LSTs
getLST = function(data){
  # UPDATE COLNAMES
  #Make sure that the column names are as follows (if using the code exactly as below) : sampleid, copyNumber, seqnames (chrom) +
  #I did not include Y chromosome for the analysis
  
  
  # GET ARM
  data$arm = apply(data,1,function(x) getArm_new(x[c("seqnames","start",'end')]))
    
  # ROUND CN VALUES
  data = data %>%
    group_by(sampleid,seqnames,arm) %>%
    mutate(cn = round(copyNumber,0))
  
  # GROUP EQUAL ONES - TAKE 1
  
  data = data %>%
    mutate(gp = paste0(seqnames,arm,cn)) 
  
  groups = getGroups(data)
  
  data = cbind(data,group = groups)
  
  data = data %>%
    group_by(sampleid,seqnames,arm,group) %>%
    mutate(newStart = start[row_number() == 1],newEnd = end[row_number() == n()]) %>%
    distinct(sampleid,newStart,newEnd,cn,group)
  
  # REMOVE LESS THAN 3MB
  data = data %>%
    mutate(newWidth = newEnd - newStart) %>%
    filter(newWidth > 3000000)
  
  # COMBINE AGAIN 
  data = data %>%
    group_by(sampleid,seqnames,arm) %>%
    mutate(gp = paste0(seqnames,arm,cn)) 
  
  groups = getGroups(data)
  
  data = cbind(data,groupNew = groups)
  
  data = data %>%
    group_by(sampleid,seqnames,arm,groupNew) %>%
    mutate(newStart = newStart[row_number() == 1],newEnd = newEnd[row_number() == n()]) %>%
    distinct(sampleid,newStart,newEnd,cn,groupNew)
  
  # GREATER THAN 10MB
  data = data %>%
    group_by(sampleid,seqnames,arm,groupNew) %>%
    mutate(newLength = newEnd - newStart) %>%
    mutate(gtr_10Mb = ifelse(newLength > 10000000, 1, 0)) %>%
    ungroup() %>%
    group_by(sampleid,seqnames,arm) %>%
    mutate(countLST = ifelse(row_number() != 1 & gtr_10Mb == 1 & (lag(gtr_10Mb) == 1),
                             1, 0))
  # COUNT LST
  lst = data %>%
    group_by(sampleid,countLST) %>% tally() %>% filter(countLST == 1) %>%
    select(sampleid,n) 
  names(lst) <- c("sampleid","lst")
  return(lst)
}

########## tAI CALCULATIONS #############
gettAI <- function(taiData){
    centromere = data.frame(c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,"X"),
                                c(1000001,1000001,1000001,2000001,1000001,1000001,1000001,1000001,
                                1000001,1000001,1000001,1000001,20000001,20000001,25000001,1000001,
                                1000001,1000001,1000001,1000001,14000001,19000001,4000001), 
                                c(121700000,91800000,87800000,48200000,46100000,58500000,58100000,43200000,
                                42200000,38000000,51000000,33200000,16500000,16100000,17500000,35300000,
                                22700000,15400000,24200000,25700000,10900000,13700000,58100000), 
                                c(125100000,96000000,94000000,51800000,51400000,62600000,62100000,
                                47200000,45500000,41600000,55800000,37800000,18900000,18200000,20500000,
                                38400000,27400000,21500000,28100000,30400000,13000000,17400000,63800000),
                                c(248000000,242000000,198000000,189000000,181000000,170000000,159000000,
                                145000000,138000000,133000000,135000000,132000000,113000000,105000000,101000000,
                                90000000,83000000,80000000,58000000,64000000,46000000,49000000,154000000))

      colnames(centromere) = c("Chr", "start_chr", "start_cen", "end_cen", "end_chr")
      
      names(taiData) <- c( "sampleid","Chr", "start", "end", "copyNumber")
      joined <- left_join(taiData, centromere, by =  "Chr")

      #selecting instable elements at telomeres not spanning the centromere and counting the number of TAI
      joined_TAI <- joined %>%
                  mutate(TAI = ifelse(((start == start_chr & end <= start_cen) | (end == end_chr & start >= end_cen)) & (end-start >= 2000000) & (copyNumber != "NEUT"), 1, 0)) 
      tai = setDT(joined_TAI %>%
          group_by(sampleid,TAI) %>% tally() %>% filter(TAI == 1) %>%
          select(sampleid,n))
      
      return(tai)
}

######## GET ALL HRD SCORE ##############
getHRDScore <- function(bins, cna, hrdLOH, sample){
      lstData <- data.table::copy(bins[,c(12,2,3,4,9)])
      names(lstData) <- c("sampleid","seqnames", "start", "end", "copyNumber")
      lst <- setDT(getLST(lstData[seqnames != "X"]))
      if(dim(lst)[1] == 0){
          lst <- data.table(sampleid = sample,lst = 0)
      }
      names(lst)[2] <- "LST"
      names(lst)[1] <- "SAMPLE"

      taiData <- data.table::copy(cna[,c(12,1,2,3,11)])
      tai <- setDT(gettAI(taiData[chr != "X",]))
      if(dim(tai)[1] == 0){
          tai <- data.table(sampleid = sample,tai = 0)
      }
      names(tai)[2] <- c("tAI")
      names(tai)[1] <- c("SAMPLE")

      hrdScore <- merge(merge(lst,tai),hrdLOH)

      hrdScore <- hrdScore[,hrdScore := LST + tAI + `HRD-LOH_15`]
      return(hrdScore)
}

# SCRIPT TO PLOT THE LOH-RATIO PLOT WITH LOH SEGMENTS
getRatioPlot <- function(data,extra = "",segData,hrdScore, selected){
    plotDataTest = unique(data[,c("SAMPLE","CHROM","i.START","i.END","logR","median_s","mad1_s","mad2_s","madRange"), with = FALSE])
    medianLines = plotDataTest %>%
        select(SAMPLE, median_s) %>%
        distinct()
    sample <- unique(medianLines[,SAMPLE])
    segData <- segData[LOH == "LOH",][,chr := gsub("chr","",CHROM)]
    hrdLOHC <- hrdScore[,`HRD-LOH_15`]
    if(selected == FALSE){
        title <- paste0(unique(data[,SAMPLE])," : ","Bin wise ratios (0.35-0.65/All)\n","HRD-LOH: ",hrdLOHC)
    }else{
      title <- paste0(unique(data[,SAMPLE])," : ","Bin wise ratios (0.35-0.65/All)\t => ","HRD-LOH: ",hrdLOHC)
    }
    chrXWindows <- fread("/ifs/data/research/projects/Samhita/srcFinal/helperFiles/chrX_windows_hg38.txt",
                    header = FALSE, stringsAsFactors = FALSE)
    names(chrXWindows) <- c("CHROM","i.START","i.END")
    chrXWindows <- chrXWindows[,`:=` (logR = NA, median_s = 0, mad1_s = 0, mad2_s = 0, madRange = NA, chr = "X", SAMPLE = sample)][,c("SAMPLE","CHROM","i.START","i.END","logR","median_s","mad1_s","mad2_s","madRange", "chr")]
    plotDataTest <- rbind(plotDataTest,chrXWindows, fill = TRUE)
    p = plotDataTest %>%
        mutate(chr = gsub("chr","",CHROM)) %>%
        ggplot(aes(y = logR,col = madRange)) +
        geom_point(aes(x = (i.START + i.END)/2), size = 1.2) +
        facet_grid(.~factor(chr,levels = c(seq(1,22,1),"X")),
                scales = "free",switch = "x",space = "free_x") +
        theme_test() +
        labs(x = NULL, y = "log2(binRatio/cohortMedian)", col = NULL, title =  title) +
        ylim(-3,3) +
        scale_x_continuous(limits = c(0,NA)) +
        geom_hline(data = medianLines, aes(yintercept = median_s), col = "black",linetype = "dashed",size = 0.3) +
        theme(legend.position = "none",
            legend.text = element_text(size = 10),
            legend.title = element_blank(),
            title = element_text(size = 10),
            panel.spacing = unit(0, "lines"),
            panel.border = element_rect(fill = NA, color = "grey", linetype = "solid",size = 0.2),
            axis.text.x = element_blank(),
            axis.ticks = element_blank(),
            axis.title.x = element_blank(),
            #axis.text.y = element_blank(),
            axis.title.y = element_text(size = 9),
            strip.placement = "inside",
            strip.text = element_text(size = 9),
            strip.text.y = element_blank(),
            #strip.background = element_rect(fill = "white"),
            strip.background = element_blank(),
            plot.margin = unit(c(0.1, 0.2, 0, 0.5), "cm")) +
        scale_color_manual(values=c("InRange" = 'grey', "Out_mad1" = "steelblue", "Out_mad2" = "steelblue","Others"="grey"),
                        labels=c("InRange" = 'Around Median', "Out_mad1" = "Median +/- MAD", "Out_mad2" = "Median +/- 2*MAD")) +
        scale_shape_manual(values = c("LOH" = 16,"Non-LOH" = 21)) 
    if(dim(segData)[1] > 0) {
        p = p + 
            geom_rect(data = segData, inherit.aes = FALSE, 
                    aes(xmin = as.numeric(START), xmax = as.numeric(END),ymin = -3, ymax = 3),
                    alpha = 0.2, fill = "orange", color = "transparent")
    }
    return(p)
}

# SCRIPT TO PLOT THE CNV RESULTS WITH HRD RESULTS
getCNVPlot_seg_bin_altSolutions <- function(data, cnvBins, results,hrdScore,pTC = NA, lohSegData, selected){
  cnCol <- c("lightskyblue4","#56B4E9","grey","plum2","tan3",rep("tan3", 25))
  names(cnCol) <- c("HOMD","HETD","NEUT","GAIN","AMP",c("HLAMP",paste0(rep("HLAMP", 8), 2:25)))
  cnvBins <- cnvBins[,chrom := chr]
  data <- data[,chrom := chr]

  purity <- 1 - as.numeric(results[,n_est])
  ploidy <- as.numeric(results[,phi_est])
  n <- results[,n]
  p <- results[,p]
  loglik <- results[,loglik]
  lstC <- hrdScore[,LST]
  taiC <- hrdScore[,tAI]
  hrdLOHC <- hrdScore[,`HRD-LOH_15`]
  hrdScoreC <- hrdScore[,hrdScore]
  ploidyC <- (1 - as.numeric(results[,n_est])) * ploidy + as.numeric(results[,n_est]) * 2
  if(selected == FALSE){
    title <- paste0(unique(data[,SAMPLE]),"\nn: ",n,"; p: ",p,"\nPurity_Est: ",purity,"; Ploidy_Est: ",ploidy,"; Pathology TC(%): ",pTC,"\n","LST: ",lstC,"; tAI: ", taiC,"; HRD-LOH: ",hrdLOHC,"; HRD-Score: ",hrdScoreC,"\n","loglik: ",loglik)
  }else{
    title <- paste0(unique(data[,SAMPLE])," => n: ",n,"; p: ",p,"; \tPurity_Est: ",purity,"; Ploidy_Est: ",ploidy,"; Pathology TC(%): ",pTC,"; \t","LST: ",lstC,"; tAI: ", taiC,"; HRD-LOH: ",hrdLOHC,"; HRD-Score: ",hrdScoreC,"; \t","loglik: ",loglik)
  }
  #data <- data[]
  p = cnvBins %>%
    ggplot(aes(col = Corrected_Call)) +
    geom_point(aes(x = (as.numeric(start) + as.numeric(end))/2, y = (as.numeric(logR)) + log2(ploidyC/2)), size = 0.5) + 
    labs(title =  title,
              y = "logR (copy number)") +
    geom_linerange(data = data[chrom != "X",],aes(xmin = start, xmax = end,y = (as.numeric(median) + log2(ploidyC/2))), size = 0.8) +
    facet_grid(.~factor(chrom,levels = c(seq(1,22,1),"X")),
               scales = "free",switch = "x",space = "free_x") +
    theme_test() +
    scale_color_manual(values = cnCol) + 
    ylim(-2,2) + 
    scale_x_continuous(limits = c(0,NA)) +
    geom_hline(yintercept = 0,col = "black",linetype = "dashed",size = 0.4) +
    theme(legend.position = "none",
          legend.text = element_text(size = 10),
          legend.title = element_blank(),
          title = element_text(size = 10),
          panel.spacing = unit(0, "lines"),
          panel.border = element_rect(fill = NA, color = "grey", linetype = "solid",size = 0.2),
          axis.text.x = element_blank(),
          axis.ticks = element_blank(),
          axis.title.x = element_blank(),
          #axis.text.y = element_blank(),
          axis.title.y = element_text(size = 9),
          strip.placement = "inside",
          strip.text = element_text(size = 9),
          strip.text.y = element_blank(),
          #strip.background = element_rect(fill = "white"),
          strip.background = element_blank(),
          plot.margin = unit(c(0.1, 0.2, 0, 0.5), "cm"))
    lohSegData <- lohSegData[LOH == "LOH",][,chrom := gsub("chr","",CHROM)]
    if(dim(lohSegData)[1] > 0) {
        p = p + 
            geom_rect(data = lohSegData, inherit.aes = FALSE, 
                    aes(xmin = as.numeric(START), xmax = as.numeric(END),ymin = -2, ymax = 2),
                    alpha = 0.2, fill = "orange", color = "transparent")
    }
  return(p)
}

############# GET THE PLOTS FOR SELECTED ichorCNA SOLUTIONS AT EACH PLOIDY RUN MODE  ###################
getSelectedSol_HRDPlots <- function(rData, sName, mergedLOH, hrdLOH, patTC = NA, mode, resDir){
    load(rData)
    loglik_man <- data.table::copy(setDT(loglik))[,c("n","p") := tstrsplit(init, "-", fixed=TRUE)]
    indSel <- order(loglik_man[, "loglik"], decreasing=TRUE)  # Get the solution with the highest loglikelihood within the mode
    n <- indSel[1]

    bins <- setDT(results[[n]]$cna[[1]])
    bins <- bins[,SAMPLE := sName]

    cna <- setDT(results[[n]]$results$segs[[1]])
    cna <- cna[,SAMPLE := sName]

    hrdScore <- getHRDScore(bins, cna, hrdLOH, sName)
    writeHRDResults(sName, loglik_man[n],hrdScore, patTC,selected = "selected", mode = mode, outDir = resDir)
    p <- getCNVPlot_seg_bin_altSolutions(cna,bins, loglik_man[n],hrdScore,patTC,mergedLOH, selected = TRUE)
  
    return(p)
}

################## GET THE PLOTS FOR ALL ichorCNA SOLUTIONS ##################
getAllSol_HRDPlots <- function(rData, sName, mergedLOH, hrdLOH, patTC = NA, resDir){
    load(rData)
    loglik_man <- data.table::copy(setDT(loglik))[,c("n","p") := tstrsplit(init, "-", fixed=TRUE)]
    indNew <- order(as.numeric(loglik_man[, factor(p)]), decreasing=FALSE) # Get the list of all the solutions from the run c(2,3,4) ordered from lowest to highest on the initial normal content estimate per initial ploidy estimate
    
  
    for(n in indNew){
      bins <- setDT(results[[n]]$cna[[1]])
      bins <- bins[,SAMPLE := sName]

      cna <- setDT(results[[n]]$results$segs[[1]])
      cna <- cna[,SAMPLE := sName]

      hrdScore <- getHRDScore(bins, cna, hrdLOH, sName)
      writeHRDResults(sName, loglik_man[n],hrdScore, patTC,selected = "all", outDir = resDir)
      p <- getCNVPlot_seg_bin_altSolutions(cna,bins, loglik_man[n],hrdScore,patTC,mergedLOH, selected = FALSE)
    
      print(p)
    }
}

########## WRITE ALL THE OUTPUT BOTH SELECTED SOLUTIONS AND ALL SOLUTIONS TO A FILE PER SAMPLE ###########
writeHRDResults <- function(sName, results, hrdScore, patTC = NA, mode = NA, selected = "all", outDir){
    purity <- 1 - as.numeric(results[,n_est])
    ploidy <- as.numeric(results[,phi_est])
    n <- results[,n]
    p <- results[,p]
    loglik <- results[,loglik]
    lstC <- hrdScore[,LST]
    taiC <- hrdScore[,tAI]
    hrdLOHC <- hrdScore[,`HRD-LOH_15`]
    hrdScoreC <- hrdScore[,hrdScore]
    
    data <- data.table(sample = sName, 
                           pathologicalTC = patTC,
                           ichorCNA_mode = mode, 
                           n0 = n, 
                           p0 = p, 
                           purityEstimated = purity, 
                           ploidyEstimated = ploidy, 
                           LST = lstC,
                           tAI = taiC,
                           `HRD-LOH` = hrdLOHC,
                           `HRD-Score` = hrdScoreC)
    path <- paste0(outDir,"/hrdResults_N/",sName,"_",selected,"Results.txt")
    writeOutput(data, path, append = TRUE)
}

############ PLOT ALL THE SOLUTIONS - SELECTED AND ALL VERSIONS #############
plotHRDSolutions <- function(sample,rDataFileList_all, rDataFileList_2, rDataFileList_3, rDataFileList_4, patTC = NA, lohCalls_bins, outDir){
        #sample = sampleList[s]
        
        ratios <- lohCalls_bins[SAMPLE == sample]
        mergedLOH <- getMergedGroups(ratios)
        hrdLOH <- getCummulativeCounts(mergedLOH,15)
        if(dim(hrdLOH)[1] == 0){
            hrdLOH <- data.table(sampleid = sample,hrdLOH_15 = 0)
        }
        names(hrdLOH)[2] <- "HRD-LOH_15"
        names(hrdLOH)[1] <- "SAMPLE"

        # Get c(2) solution
        rDataFile <- getFiles(rDataFileList_2,sample)
        print(outDir)
        p1 <- getSelectedSol_HRDPlots(rDataFile, sample, mergedLOH, hrdLOH, patTC = patTC, mode = "c(2)", resDir = outDir)

        # Get c(3) solution
        rDataFile <- getFiles(rDataFileList_3,sample)
        p2 <- getSelectedSol_HRDPlots(rDataFile, sample, mergedLOH, hrdLOH, patTC = patTC, mode = "c(3)", resDir = outDir)

        # Get c(4) solution
        rDataFile <- getFiles(rDataFileList_4,sample)
        p3 <- getSelectedSol_HRDPlots(rDataFile, sample, mergedLOH, hrdLOH, patTC = patTC, mode = "c(4)", resDir = outDir)

        # Get c(2,3,4) solution
        rDataFile <- getFiles(rDataFileList_all,sample)
        p4 <- getSelectedSol_HRDPlots(rDataFile, sample, mergedLOH, hrdLOH, patTC = patTC, mode = "c(2,3,4)", resDir = outDir)

        pRatios <- getRatioPlot(ratios,"0.35-0.65/All",mergedLOH,hrdLOH, selected = TRUE)
        
        pAll <- cowplot::plot_grid(p1,p2,p3,p4,pRatios,nrow = 5,align = "v")

        pdf(paste0(outDir,"/hrdPlots_N/",sample,"_selectedSolutions.pdf"),
        onefile = TRUE,height = 14, width = 19)

        print(pAll)

        dev.off()

        # Print all solutions
        pdf(paste0(outDir,"/hrdPlots_N/",sample,"_allSolutions.pdf"),
        onefile = TRUE,height = 4, width = 14)
        
        rDataFile <- getFiles(rDataFileList_all,sample)
        getAllSol_HRDPlots(rDataFile, sample, mergedLOH, hrdLOH, patTC = patTC, resDir = outDir)
        pRatios <- getRatioPlot(ratios,"0.35-0.65/All",mergedLOH,hrdLOH, selected = FALSE)
        print(pRatios)

        dev.off()
}

# WRITE OUTPUT
writeOutput = function(data, path, append = FALSE){
    fwrite(data,
            path,
            quote = FALSE, row.names = FALSE, sep = "\t", append = append)
}
