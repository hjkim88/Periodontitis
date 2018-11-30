###
#   File name : TestAffyData.R
#   Author    : Hyunjin Kim
#   Date      : Nov 28, 2018
#   Email     : hk2990@cumc.columbia.edu
#   Purpose   : Because there are so many DE genes with significantly small FDRs,
#               we want to do sanity checks with the data
#
#   * 1. What if I do not consider patient info on linear modeling?
#   * 2. How about DE genes after combining all the replicates?
#   * 3. How about DE analysis with random sample info?
#   * 4. Does the size of sample affects FDRs of DE genes?
#   * 5. What if I personally perform my own t-test on the data?
#   * 6. What if I pre-process and normalize the raw data on my own? Will there be any difference?
#
#   Instruction
#               1. Source("TestRNASEQData.R")
#               2. Run the function "sanityChecks" - specify an input file path (normalized microarray)
#               3. All the results will be presented in the console
#
#   Example
#               > source("The_directory_of_TestRNASEQData.R/TestRNASEQData.R")
#               > sanityChecks(rawDataDir="E:/Panos/Affy/CEL/",
#                              rawDataSampleInfoPath="E:/Panos/Affy/PatientsSamplesBiopsiesExperimentsSummary_rf2.txt",
#                              normDataPath="./data/Affy/panos_affy_data.rda")
###

sanityChecks <- function(rawDataDir="//isilon.c2b2.columbia.edu/ifs/archive/shares/bisr/Papapanou/R21/OriginalData/mRNA/LinkTomRNAData/",
                         rawDataSampleInfoPath="//isilon.c2b2.columbia.edu/ifs/archive/shares/bisr/Papapanou/R21/rich_aim_1/PatientsSamplesBiopsiesExperimentsSummary_rf2.txt",
                         normDataPath="//isilon.c2b2.columbia.edu/ifs/archive/shares/bisr/Papapanou/Sept_2018/PreprocessedData/Affy/panos_affy_data.rda") {
  
  ### load dataset
  load(normDataPath)
  
  
  ### Now the samples have replicates. Make them one sample for one person
  
  ### get indicies of periodontitis and control
  idx_perio <- which(affy_norm_ge_sample_info$Condition == "Affected")
  idx_control <- which(affy_norm_ge_sample_info$Condition == "Control")
  
  ### get unique patients
  unique_patients_perio <- unique(affy_norm_ge_sample_info$Patient[which(affy_norm_ge_sample_info$Condition == "Affected")])
  unique_patients_control <- unique(affy_norm_ge_sample_info$Patient[which(affy_norm_ge_sample_info$Condition == "Control")])
  
  ### average the replicates from the same patient - periodontitis
  affy_norm_ge_perio <- matrix(NA, nrow(affy_norm_ge), length(unique_patients_perio))
  for(i in 1:length(unique_patients_perio)) {
    affy_norm_ge_perio[,i] <- apply(affy_norm_ge[intersect(idx_perio, which(affy_norm_ge_sample_info$Patient == unique_patients_perio[i]))+1],
                                    1, mean)
  }
  rownames(affy_norm_ge_perio) <- rownames(affy_norm_ge)
  colnames(affy_norm_ge_perio) <- paste0("P", unique_patients_perio)
  affy_norm_ge_perio <- data.frame(affy_norm_ge_perio, stringsAsFactors = FALSE, check.names = FALSE)
  
  ### average the replicates from the same patient - healthy
  affy_norm_ge_control <- matrix(NA, nrow(affy_norm_ge), length(unique_patients_control))
  for(i in 1:length(unique_patients_control)) {
    affy_norm_ge_control[,i] <- apply(affy_norm_ge[intersect(idx_control, which(affy_norm_ge_sample_info$Patient == unique_patients_control[i]))+1],
                                    1, mean)
  }
  rownames(affy_norm_ge_control) <- rownames(affy_norm_ge)
  colnames(affy_norm_ge_control) <- paste0("C", unique_patients_control)
  affy_norm_ge_control <- data.frame(affy_norm_ge_control, stringsAsFactors = FALSE, check.names = FALSE)
  
  ### combine the periodontitis and control samples - duplicate-averaged
  affy_norm_ge_per_patient <- data.frame(affy_norm_ge_perio, affy_norm_ge_control,
                                         stringsAsFactors = FALSE, check.names = FALSE)
  
  ### make a new sample info for the duplicate-averaged data
  affy_norm_ge_per_patient_sample_info <- data.frame(Patient=c(unique_patients_perio, unique_patients_control),
                                                     Penotype=c(rep("Periodontitis", length(unique_patients_perio)),
                                                                rep("Control", length(unique_patients_control))),
                                                     stringsAsFactors = FALSE, check.names = FALSE)
  
  
  ### A function to perform repetitive DE analysis with limma
  limmaWithComparisons <- function(normCnt, grp, exp_class, ctrl_class, bat_eff=NULL) {
    
    ### load library
    if(!require(limma)) {
      source("https://bioconductor.org/biocLite.R")
      biocLite("limma")
      library(limma)
    }
    
    ### sometimes, there are some variables which can not be transformed into R variable names
    ### so, just to be safe, change all the variables to R-usuable ones
    grp <- make.names(grp)
    exp_class <- make.names(exp_class)
    ctrl_class <- make.names(ctrl_class)
    
    ### make a design matrix for DE analysis
    sampleType <- relevel(as.factor(grp), ref = ctrl_class)
    if(is.null(bat_eff)) {
      design <- model.matrix(~0+sampleType)
      colnames(design) <- levels(sampleType)
    } else {
      bat_eff <- make.names(bat_eff)
      bat_eff <- as.factor(bat_eff)
      design <- model.matrix(~0+sampleType+bat_eff)
      colnames(design) <- c(levels(sampleType), levels(bat_eff)[-1])
    }
    
    ### fir the linear model
    fit <- lmFit(normCnt, design)
    
    ### extract specific comparison of interest
    contrastMat <- makeContrasts(contrasts=paste(exp_class,ctrl_class,sep="-"), levels=design)
    
    ### fit the contrasts
    fit2 <- contrasts.fit(fit, contrastMat)
    fit2 <- eBayes(fit2)
    
    ### get the differentially expressed genes
    result <- topTable(fit2, adjust.method="BH", number=Inf)
    
    ### order based on adj.p.val
    result <- result[order(result$adj.P.Val),]
    
    return(result)
  }
  
  
  ### Original DE analysis
  deresult1 <- limmaWithComparisons(normCnt = affy_norm_ge[,-1],
                                    grp = affy_norm_ge_sample_info$Condition,
                                    exp_class = "Affected",
                                    ctrl_class = "Control",
                                    bat_eff = affy_norm_ge_sample_info$Patient)
  
  
  ### 1. What if I do not consider patient info on linear modeling?
  ### DE analysis without considering patient info
  deresult2 <- limmaWithComparisons(normCnt = affy_norm_ge[,-1],
                                    grp = affy_norm_ge_sample_info$Condition,
                                    exp_class = "Affected",
                                    ctrl_class = "Control",
                                    bat_eff = NULL)
  
  
  ### 2. How about DE genes after combining all the replicates?
  ### DE analysis after merging replicates into one
  deresult3 <- limmaWithComparisons(normCnt = affy_norm_ge_per_patient,
                                    grp = affy_norm_ge_per_patient_sample_info$Penotype,
                                    exp_class = "Periodontitis",
                                    ctrl_class = "Control",
                                    bat_eff = affy_norm_ge_per_patient_sample_info$Patient)
  
  
  ### 2. How about DE genes after combining all the replicates?
  ### DE analysis after merging replicates into one
  ### And without considering patient info
  deresult4 <- limmaWithComparisons(normCnt = affy_norm_ge_per_patient,
                                    grp = affy_norm_ge_per_patient_sample_info$Penotype,
                                    exp_class = "Periodontitis",
                                    ctrl_class = "Control",
                                    bat_eff = NULL)
  
  
  ### 3. How about DE analysis with random sample info?
  set.seed(1234)
  deresult5 <- limmaWithComparisons(normCnt = affy_norm_ge[,-1],
                                    grp = affy_norm_ge_sample_info$Condition[sample(nrow(affy_norm_ge_sample_info), nrow(affy_norm_ge_sample_info))],
                                    exp_class = "Affected",
                                    ctrl_class = "Control",
                                    bat_eff = NULL)
  
  
  ### 4. Does the size of sample affects FDRs of DE genes?
  num_for_each_class <- 5
  new_affy_norm_ge <- affy_norm_ge[,(c(idx_perio[sample(length(idx_perio), num_for_each_class)],
                                       idx_control[sample(length(idx_control), num_for_each_class)])+1)]
  deresult6 <- limmaWithComparisons(normCnt = new_affy_norm_ge,
                                    grp = c(rep("Affected", num_for_each_class),
                                            rep("Control", num_for_each_class)),
                                    exp_class = "Affected",
                                    ctrl_class = "Control",
                                    bat_eff = NULL)
  
  num_for_each_class <- 10
  new_affy_norm_ge <- affy_norm_ge[,(c(idx_perio[sample(length(idx_perio), num_for_each_class)],
                                      idx_control[sample(length(idx_control), num_for_each_class)])+1)]
  deresult7 <- limmaWithComparisons(normCnt = new_affy_norm_ge,
                                    grp = c(rep("Affected", num_for_each_class),
                                            rep("Control", num_for_each_class)),
                                    exp_class = "Affected",
                                    ctrl_class = "Control",
                                    bat_eff = NULL)
  
  num_for_each_class <- 20
  new_affy_norm_ge <- affy_norm_ge[,(c(idx_perio[sample(length(idx_perio), num_for_each_class)],
                                       idx_control[sample(length(idx_control), num_for_each_class)])+1)]
  deresult8 <- limmaWithComparisons(normCnt = new_affy_norm_ge,
                                    grp = c(rep("Affected", num_for_each_class),
                                            rep("Control", num_for_each_class)),
                                    exp_class = "Affected",
                                    ctrl_class = "Control",
                                    bat_eff = NULL)
  
  num_for_each_class <- 30
  new_affy_norm_ge <- affy_norm_ge[,(c(idx_perio[sample(length(idx_perio), num_for_each_class)],
                                       idx_control[sample(length(idx_control), num_for_each_class)])+1)]
  deresult9 <- limmaWithComparisons(normCnt = new_affy_norm_ge,
                                    grp = c(rep("Affected", num_for_each_class),
                                            rep("Control", num_for_each_class)),
                                    exp_class = "Affected",
                                    ctrl_class = "Control",
                                    bat_eff = NULL)
  
  num_for_each_class <- 40
  new_affy_norm_ge <- affy_norm_ge[,(c(idx_perio[sample(length(idx_perio), num_for_each_class)],
                                       idx_control[sample(length(idx_control), num_for_each_class)])+1)]
  deresult10 <- limmaWithComparisons(normCnt = new_affy_norm_ge,
                                    grp = c(rep("Affected", num_for_each_class),
                                            rep("Control", num_for_each_class)),
                                    exp_class = "Affected",
                                    ctrl_class = "Control",
                                    bat_eff = NULL)
  
  num_for_each_class <- 50
  new_affy_norm_ge <- affy_norm_ge[,(c(idx_perio[sample(length(idx_perio), num_for_each_class)],
                                       idx_control[sample(length(idx_control), num_for_each_class)])+1)]
  deresult11 <- limmaWithComparisons(normCnt = new_affy_norm_ge,
                                    grp = c(rep("Affected", num_for_each_class),
                                            rep("Control", num_for_each_class)),
                                    exp_class = "Affected",
                                    ctrl_class = "Control",
                                    bat_eff = NULL)
  
  num_for_each_class <- 60
  new_affy_norm_ge <- affy_norm_ge[,(c(idx_perio[sample(length(idx_perio), num_for_each_class)],
                                       idx_control[sample(length(idx_control), num_for_each_class)])+1)]
  deresult12 <- limmaWithComparisons(normCnt = new_affy_norm_ge,
                                     grp = c(rep("Affected", num_for_each_class),
                                             rep("Control", num_for_each_class)),
                                     exp_class = "Affected",
                                     ctrl_class = "Control",
                                     bat_eff = NULL)
  
  
  ### the number of DE genes filtered with the cut-offs
  ge_num_0.05 <- c(length(which(deresult6$adj.P.Val < 0.05)),
                   length(which(deresult7$adj.P.Val < 0.05)),
                   length(which(deresult8$adj.P.Val < 0.05)),
                   length(which(deresult9$adj.P.Val < 0.05)),
                   length(which(deresult10$adj.P.Val < 0.05)),
                   length(which(deresult11$adj.P.Val < 0.05)),
                   length(which(deresult12$adj.P.Val < 0.05)),
                   length(which(deresult2$adj.P.Val < 0.05)))
  names(ge_num_0.05) <- c("5_vs_5",
                          "10_vs_10",
                          "20_vs_20",
                          "30_vs_30",
                          "40_vs_40",
                          "50_vs_50",
                          "60_vs_60",
                          "The_Original")
  ge_num_0.01 <- c(length(which(deresult6$adj.P.Val < 0.01)),
                   length(which(deresult7$adj.P.Val < 0.01)),
                   length(which(deresult8$adj.P.Val < 0.01)),
                   length(which(deresult9$adj.P.Val < 0.01)),
                   length(which(deresult10$adj.P.Val < 0.01)),
                   length(which(deresult11$adj.P.Val < 0.01)),
                   length(which(deresult12$adj.P.Val < 0.01)),
                   length(which(deresult2$adj.P.Val < 0.01)))
  names(ge_num_0.01) <- c("5_vs_5",
                          "10_vs_10",
                          "20_vs_20",
                          "30_vs_30",
                          "40_vs_40",
                          "50_vs_50",
                          "60_vs_60",
                          "The_Original")
  
  ### draw bar plots
  barplot(ge_num_0.05)
  barplot(ge_num_0.01)
  
  
  ### 5. What if I personally perform my own t-test on the data?
  
  ### A function to perform t-test for DEA
  ### normGE = normalized gene expression matrix, row = genes, column = samples
  ### idx1 = column indices of condition1 (e.g., tumor)
  ### idx2 = column indices of condition2 (e.g., control)
  simple_t_test <- function(normGE, idx1, idx2) {
    
    ### compute mean values, fold changes, and t-statistic values
    mean <- 0
    lfc <- 0
    t <- 0
    for(i in 1:nrow(normGE)) {
      mean[i] <- mean(as.numeric(normGE[i,]))
      lfc[i] <- mean(as.numeric(normGE[i,idx1])) - mean(as.numeric(normGE[i,idx2]))
      t[i] <- lfc[i] / sqrt((sd(normGE[i,idx1])^2)/length(idx1) + (sd(normGE[i,idx2])^2)/length(idx2))
    }
    
    ### compute p-value
    p <- 2*pt(-abs(t), df=length(union(idx1, idx2))-2)
    
    ### correct the p-value with Bonferroni correction
    bon <- 1-((1-p)^nrow(normGE))
    
    ### A function to correct p-values with Benjamini-Hochberg approach
    correct_bh <- function(pvs) {
      
      temp <- cbind(p=pvs, No=seq(1:length(pvs)))
      
      if(length(which(is.na(temp[,"p"]))) > 0) {
        temp[which(is.na(temp[,"p"])),"p"] <- 1
      }
      
      temp <- temp[order(temp[,"p"]),]
      temp <- cbind(temp, Rank=seq(1:length(pvs)), BH=1)
      
      
      temp[length(pvs), "BH"] <- temp[length(pvs), "p"]
      for(i in (length(pvs)-1):1) {
        temp[i,"BH"] <- min(temp[i+1, "BH"], temp[i,"p"]*length(pvs)/temp[i,"Rank"])
      }
      
      temp <- temp[order(temp[,"No"]),]
      
      return(as.numeric(temp[,"BH"]))
    }  
    
    ### correct the p-value with Benjamini-Hochberg correction (FDR)
    bh <- correct_bh(p)
    
    ### combine mean, fold changes, t values, p-values, bonferroni, and benjamini-hochberg
    result <- data.frame(Gene=rownames(normGE), mean=mean, log2FC=lfc, t=t, pVal=p, bon_pVal=bon, bh_pVal=bh,
                         stringsAsFactors = FALSE, check.names = FALSE)
    
    ### sort the result in order of bh p-value
    result <- result[order(result$bh_pVal, result$pVal),]
    
    ### change numeric columns as numeric
    result[2:ncol(result)] <- lapply(result[2:ncol(result)], function(x) as.numeric(x))
    
    return(result)
  }
  
  ### do DE analysis with simple t-test
  deresult13 <- simple_t_test(affy_norm_ge[,-1], idx_perio, idx_control)
  rownames(deresult13) <- deresult13$Gene
  
  ### plot FDR correlation
  common_probes <- intersect(rownames(deresult2), rownames(deresult13)) 
  plot(deresult2[common_probes, "adj.P.Val"], deresult13[common_probes, "bh_pVal"],
       xlab = "Limma 235 vs 69 - adj.P.Val",
       ylab = "Simple t-test 235 vs 69 - P.Val(BH corrected)")
  
  
  ### 6. What if I pre-process and normalize the raw data on my own? Will there be any difference?
  
  ### load library
  if(!require(affy)) {
    source("https://bioconductor.org/biocLite.R")
    biocLite("affy")
    library(affy)
  }
  if(!require(sva)) {
    source("https://bioconductor.org/biocLite.R")
    biocLite("sva")
    library(sva)
  }
  
  
  ### read CEL files
  data <- ReadAffy(celfile.path = rawDataDir)
  
  
  ### RMA normalization
  eset <- rma(data)
  
  
  ### filter the probes
  eset <- featureFilter(eset, require.entrez=FALSE, remove.dupEntrez=FALSE)
  
  
  ### get the normalized gene expressions as a data frame
  df <- data.frame(exprs(eset), stringsAsFactors = FALSE, check.names = FALSE)
  
  
  ### load sample info
  sample_info <- read.table(file = rawDataSampleInfoPath, header = TRUE, sep = "\t",
                            stringsAsFactors = FALSE, check.names = FALSE)
  
  
  ### keep the samples only in the sample info
  df <- df[,sample_info$mRNA]
  
  
  ### A function to perform 2D PCA and save a plot
  pca_plot <- function(normalizedMat, grp, title, component=c("PC1&PC2", "PC2&PC3"), filePath=NULL) {
    
    ### load library
    if(!require(ggfortify)) {
      install.packages("ggfortify")
      library(ggfortify)
    }
    
    ### PCA
    pca_result <- prcomp(t(normalizedMat))
    pca_group <- data.frame(pca_result$x, group=grp)
    
    colors = rainbow(length(unique(grp)))
    names(colors) = unique(grp)
    
    ### plot
    if(!is.null(filePath)) {
      png(filename=filePath, width = 1000, height = 800)
    }
    if(component[1] == "PC1&PC2") {
      print(ggplot(pca_group,aes(x=PC1,y=PC2,col=group)) +
              labs(title=title) +
              geom_text(aes(label=colnames(normalizedMat)),hjust=0, vjust=0) +
              scale_color_manual(values = colors) +
              theme_classic(base_size = 16))
    } else if(component[1] == "PC2&PC3") {
      print(ggplot(pca_group,aes(x=PC2,y=PC3,col=group)) +
              labs(title=title) +
              geom_text(aes(label=colnames(normalizedMat)),hjust=0, vjust=0) +
              scale_color_manual(values = colors) +
              theme_classic(base_size = 16))
    } else {
      stop("\"component\" parameter should be \"PC1&PC2\" or \"PC2&PC3\"")
    }
    if(!is.null(filePath)) {
      dev.off()
    }
  }
  
  
  ### pca plot
  pca_plot(normalizedMat = df, grp = sample_info$Condition, title = "Before Batch Correction - Phenotype")
  
  
  ### DE analysis
  deresult14 <- limmaWithComparisons(normCnt = df,
                                     grp = sample_info$Condition,
                                     exp_class = "Affected",
                                     ctrl_class = "Control",
                                     bat_eff = NULL)
  
  
  ### add scan date to the sample info
  scan_date <- eset@protocolData@data$ScanDate
  scan_date <- strsplit(scan_date, split = " ", fixed = TRUE)
  scan_date <- sapply(scan_date, function(x) x[1])
  scan_date <- strsplit(scan_date, split = "/", fixed = TRUE)
  scan_date <- sapply(scan_date, function(x) x[3])
  scan_date <- as.numeric(scan_date)
  names(scan_date) <- colnames(exprs(eset))
  sample_info <- data.frame(sample_info, ScanYear=scan_date[sample_info$mRNA],
                            stringsAsFactors = FALSE, check.names = FALSE)
  
  
  ### add batch info based on the scan year (4/5: Batch1, 6/7/8: Batch2)
  sample_info <- data.frame(sample_info, Batch="", stringsAsFactors = FALSE, check.names = FALSE)
  sample_info$Batch[which(sample_info$ScanYear %in% c(4, 5))] <- "Batch1"
  sample_info$Batch[which(sample_info$ScanYear %in% c(6, 7, 8))] <- "Batch2"
  
  
  ### pca plot
  pca_plot(normalizedMat = df, grp = as.character(sample_info$ScanYear),
           title = "Before Batch Correction - ScanYear")
  pca_plot(normalizedMat = df, grp = as.character(sample_info$Batch),
           title = "Before Batch Correction - Batch")
  
  
  ### batch effect correction
  new_df <- ComBat(dat = df, batch = as.character(sample_info$Batch),
                   mod = model.matrix(~sample_info$Condition))
  
  
  ### pca plot after batch effect correction
  pca_plot(normalizedMat = new_df, grp = sample_info$Condition, title = "After Batch Correction - Phenotype")
  pca_plot(normalizedMat = new_df, grp = as.character(sample_info$Batch),
           title = "After Batch Correction - Batch")
  
  
  ### DE analysis after batch effect correction
  deresult15 <- limmaWithComparisons(normCnt = new_df,
                                     grp = sample_info$Condition,
                                     exp_class = "Affected",
                                     ctrl_class = "Control",
                                     bat_eff = NULL)
  
}


