###
#   File name : LevelCorrelations2.R
#   Author    : Hyunjin Kim
#   Date      : Nov 21, 2018
#   Email     : hk2990@cumc.columbia.edu
#   Purpose   : Make a correlation plot based on levels between methylation and gene expression
#
#   * THIS CODE IS DIFFERENT FROM "LevelCorrelation.R", SINCE THIS INVOLVES WITH
#     BOTH OLD (450K METHYLATION & AFFY CHIP) AND NEW (EPIC & RNA-SEQ) DATASETS
#     THE PLOTS GENERATED HERE HAVE TWO PLOTS FROM THE OLD & THE NEW
#
#   When matching CpGs to genes, there may be many CpGs associated with one gene
#   In that case, choose one which has the most negative spearman correlation
#   Therefore, the comparison is done with the matching genes
#
#   Instruction
#               1. Source("LevelCorrelations2.R")
#               2. Run the function "levelCor2" - specify the input files (levels of methylation and gene expression) and output directory
#               3. A correlation plot will be generated in the output directory
#
#   Example
#               > source("The_directory_of_LevelCorrelations2.R/LevelCorrelations2.R")
#               > levelCor2(epicLevelPath="./data/Methylation/beta_values_gene.txt",
#                           rnaseqLevelPath="./data/RNA_Seq/raw_counts.txt",
#                           methyl450kLevelPath="./data/450k/450k_m_values_gene.txt",
#                           affyLevelPath="./data/Affy/panos_affy_data.rda",
#                           newSampleInfoPath="./data/RNA_Seq/sampleInfo.txt",
#                           importantGenePath="./results/Combined/List_-log10(FDR)_DMR_Periodontitis_Healthy_all.xlsx",
#                           fdrThreshold=0.05,
#                           outputDir="./results/Cor/Combined/")
###

levelCor2 <- function(epicLevelPath="//isilon.c2b2.columbia.edu/ifs/archive/shares/bisr/Papapanou/Sept_2018/PreprocessedData/Methylation/m_values_gene.txt",
                      rnaseqLevelPath="//isilon.c2b2.columbia.edu/ifs/archive/shares/bisr/Papapanou/Sept_2018/PreprocessedData/RNA_Seq/raw_counts.txt",
                      methyl450kLevelPath="//isilon.c2b2.columbia.edu/ifs/archive/shares/bisr/Papapanou/Sept_2018/PreprocessedData/450k/450k_m_values_gene.txt",
                      affyLevelPath="//isilon.c2b2.columbia.edu/ifs/archive/shares/bisr/Papapanou/Sept_2018/PreprocessedData/Affy/panos_affy_data.rda",
                      newSampleInfoPath="//isilon.c2b2.columbia.edu/ifs/archive/shares/bisr/Papapanou/Sept_2018/PreprocessedData/RNA_Seq/sampleInfo.txt",
                      importantGenePath="//isilon.c2b2.columbia.edu/ifs/archive/shares/bisr/Papapanou/Sept_2018/Combined/FDR_Plot/List_-log10(FDR)_DMR_Periodontitis_Healthy_all.xlsx",
                      fdrThreshold=0.05,
                      outputDir="//isilon.c2b2.columbia.edu/ifs/archive/shares/bisr/Papapanou/Sept_2018/Combined/Level_Correlation/Affy_vs_450k/") {
  
  ### load data
  epicLev <- read.table(file = epicLevelPath, header = TRUE, sep = "\t", check.names = FALSE, row.names = 1)
  rnaseqLev <- read.table(file = rnaseqLevelPath, header = TRUE, sep = "\t", check.names = FALSE)
  newSampleInfo <- read.table(file = newSampleInfoPath, header = TRUE, sep = "\t", check.names = FALSE, row.names = 1)
  methyl450kLev <- read.table(file = methyl450kLevelPath, header = TRUE, sep = "\t", check.names = FALSE, row.names = 1, stringsAsFactors = FALSE)
  load(affyLevelPath)
  
  
  ### A function to transform RNA-Seq data with VST in DESeq2 package
  normalizeRNASEQwithVST <- function(readCount, filtering=TRUE) {
    
    ### load library
    if(!require(DESeq2)) {
      source("https://bioconductor.org/biocLite.R")
      biocLite("DESeq2")
      library(DESeq2)
    }
    
    ### make a design matrix for DESeq2 data
    condition <- data.frame(factor(rep("OneClass", ncol(readCount))))
    
    ### Data preparation for DESeq2 format
    deSeqData <- DESeqDataSetFromMatrix(countData=readCount, colData=condition, design= ~0)
    
    if(filtering == TRUE) {
      ### Remove rubbish rows - this will decrease the number of rows
      deSeqData <- deSeqData[rowSums(counts(deSeqData))>1,]
    }
    
    ### VST
    vsd <- vst(deSeqData)
    transCnt <- data.frame(assay(vsd), check.names = FALSE)
    
    return (transCnt)
    
  }
  
  
  ### normalize the raw counts
  rnaseqLev <- normalizeRNASEQwithVST(rnaseqLev, filtering = FALSE)
  
  
  ### make two datasets having the same format
  ### NEW
  common_genes <- intersect(rownames(epicLev), rownames(rnaseqLev))
  common_samples <- intersect(colnames(epicLev), colnames(rnaseqLev))
  rnaseqLev <- rnaseqLev[common_genes,common_samples]
  epicLev <- epicLev[common_genes,common_samples]
  ### OLD
  common_genes <- intersect(rownames(methyl450kLev), rownames(affy_norm_ge))
  common_samples <- intersect(colnames(methyl450kLev), colnames(affy_norm_ge))
  methyl450kLev <- methyl450kLev[common_genes,common_samples]
  affy_norm_ge <- affy_norm_ge[common_genes,common_samples]
  
  
  ### only select Periodontitis samples - for the NEW data
  ### old methylation (450k) only has periodontitis samples, thus no need for OLD data
  common_samples <- intersect(colnames(epicLev), colnames(rnaseqLev))
  newSampleInfo <- newSampleInfo[common_samples,]
  rnaseqLev_periodontitis <- rnaseqLev[,which(newSampleInfo$Phenotype == "Periodontitis")]
  epic_periodontitis <- epicLev[,which(newSampleInfo$Phenotype == "Periodontitis")]
  
  
  ### correlation measuring
  
  ### load library
  if(!require(xlsx)) {
    install.packages("xlsx")
    library(xlsx)
  }
  if(!require(ggplot2)) {
    install.packages("ggplot2")
    library(ggplot2)
  }
  if(!require(gridExtra)) {
    install.packages("gridExtra")
    library(gridExtra)
  }
  
  ### get important genes
  gList <- read.xlsx2(file = importantGenePath, sheetIndex = 1, stringsAsFactors = FALSE)
  gList[,3:4] <- apply(gList[,3:4], 2, as.numeric)
  gList <- gList$Gene_Name[intersect(which(gList$FDR_Expression < fdrThreshold),
                                     which(gList$FDR_Methylation < fdrThreshold))]
  gList <- unique(gList)
  
  ### make a background distribution for p-value generation
  
  ### NEW
  new_cor_dist <- NULL
  for(i in 1:nrow(rnaseqLev_periodontitis)) {
    new_cor_dist <- c(new_cor_dist, cor(as.numeric(rnaseqLev_periodontitis[i,]),
                                        as.numeric(epic_periodontitis[i,]),
                                        use = "pairwise.complete.obs"))
  }
  names(new_cor_dist) <- rownames(rnaseqLev_periodontitis)
  
  ### set NA to 0
  naIdx <- which(is.na(new_cor_dist))
  if(length(naIdx) > 0) {
    new_cor_dist[naIdx] <- 0
  }
  
  ### sort in ascending order
  new_cor_dist <- new_cor_dist[order(new_cor_dist)]
  
  ### OLD - Gene Symbol base
  old_cor_dist <- NULL
  old_unique_genes <- unique(affy_norm_ge$Gene_Symbol)
  for(i in 1:length(old_unique_genes)) {
    ### get all the rows with the gene name
    dupsIdx <- which(affy_norm_ge$Gene_Symbol == old_unique_genes[i])
    
    ### combine all the values with the same gene name to measure correlation
    temp_cor_data <- NULL
    for(j in 1:length(dupsIdx)) {
      temp_cor_data <- rbind(temp_cor_data, cbind(t(affy_norm_ge[dupsIdx[j],-1]), t(methyl450kLev[dupsIdx[j],-1])))
    }
    
    ### correlation calculation
    old_cor_dist <- c(old_cor_dist, cor(temp_cor_data[,1],
                                        temp_cor_data[,2],
                                        use = "pairwise.complete.obs"))
  }
  names(old_cor_dist) <- old_unique_genes
  
  ### set NA to 0
  naIdx <- which(is.na(old_cor_dist))
  if(length(naIdx) > 0) {
    old_cor_dist[naIdx] <- 0
  }
  
  ### sort in ascending order
  old_cor_dist <- old_cor_dist[order(old_cor_dist)]
  
  ### a function to retreive one-tale empirical p-value (negative: small p-value, positive: large p-value)
  getEmpiricalPV <- function(target_cor, cor_dist) {
    return(length(which(cor_dist < target_cor)) / length(cor_dist))
  }
  
  ### iteratively performes generating correlation plots
  for(i in 1:length(gList)) {
    ### make a source data for a plot
    ### NEW
    cor_data1 <- cbind(t(rnaseqLev_periodontitis[gList[i],]), t(epic_periodontitis[gList[i],]))
    colnames(cor_data1) <- c("RNASEQ_NormCnt_Perio", "EPIC_Mvalue_Perio")
    ### OLD
    dupsIdx <- which(affy_norm_ge$Gene_Symbol == gList[i])
    ### since the important genes are from the new data, it may not exist in the old data
    if(length(dupsIdx) > 0) {
      cor_data2 <- NULL
      for(j in 1:length(dupsIdx)) {
          cor_data2 <- rbind(cor_data2, cbind(t(affy_norm_ge[dupsIdx[j],-1]), t(methyl450kLev[dupsIdx[j],-1])))
      }
      colnames(cor_data2) <- c("AFFY_NormExp_Perio", "Methyl450K_Mvalue_Perio")
      
      ### create a plot for periodontitis samples
      fName <- paste0(gList[i], "_Periodontitis_Correlation_Old_New.png")
      oldPlot <- ggplot(data = as.data.frame(cor_data2), aes(x=AFFY_NormExp_Perio, y=Methyl450K_Mvalue_Perio)) +
                 geom_point(color = "black", size = 1) +
                 labs(title=paste0(substr(fName, 1, nchar(fName)-11), "Old"),
                      subtitle=sprintf("P.Cor = %s, Empirical.p.val = %s",
                                  round(old_cor_dist[gList[i]], 5),
                                  round(getEmpiricalPV(old_cor_dist[gList[i]], old_cor_dist), 5))) +
                 xlab("Expression Levels") +
                 ylab("Methylation Levels") +
                 geom_smooth(method = lm, color="blue", se=FALSE) +
                 theme_classic(base_size = 16)
      newPlot <- ggplot(data = as.data.frame(cor_data1), aes(x=RNASEQ_NormCnt_Perio, y=EPIC_Mvalue_Perio)) +
                 geom_point(color = "black", size = 1) +
                 labs(title=paste0(substr(fName, 1, nchar(fName)-11), "New"),
                      subtitle=sprintf("P.Cor = %s, Empirical.p.val = %s",
                                  round(new_cor_dist[gList[i]], 5),
                                  round(getEmpiricalPV(new_cor_dist[gList[i]], new_cor_dist), 5))) +
                 xlab("Expression Levels") +
                 ylab("Methylation Levels") +
                 geom_smooth(method = lm, color="blue", se=FALSE) +
                 theme_classic(base_size = 16)
      ggsave(filename = paste0(outputDir, fName), arrangeGrob(oldPlot, newPlot, ncol = 2), width = 20, height = 10)
    }
  }
  
  ### correlation of all the genes - for a reference
  ### NEW
  cor_data1 <- NULL
  for(i in 1:nrow(rnaseqLev_periodontitis)) {
    cor_data1 <- rbind(cor_data1, cbind(t(rnaseqLev_periodontitis[i,]), t(epic_periodontitis[i,])))
  }
  colnames(cor_data1) <- c("RNASEQ_NormCnt_Perio", "EPIC_Mvalue_Perio")
  ### OLD
  cor_data2 <- NULL
  for(i in 1:nrow(affy_norm_ge)) {
    cor_data2 <- rbind(cor_data2, cbind(t(affy_norm_ge[i,-1]), t(methyl450kLev[i,-1])))
  }
  colnames(cor_data2) <- c("AFFY_NormExp_Perio", "Methyl450K_Mvalue_Perio")
  
  ### make a plot 
  fName <- "All_Genes_Periodontitis_Correlation_Old_New.png"
  oldPlot <- ggplot(data = as.data.frame(cor_data2), aes(x=AFFY_NormExp_Perio, y=Methyl450K_Mvalue_Perio)) +
             geom_point(color = "black", size = 1) +
             labs(title=paste0(substr(fName, 1, nchar(fName)-11), "Old"),
                  subtitle=sprintf("P.Cor = %s, p-value = %s",
                                   round(cor(cor_data2[,"AFFY_NormExp_Perio"], cor_data2[,"Methyl450K_Mvalue_Perio"], use = "pairwise.complete.obs"), 5),
                                   signif(cor.test(cor_data2[,"AFFY_NormExp_Perio"], cor_data2[,"Methyl450K_Mvalue_Perio"])$p.value, 5))) +
             xlab("Expression Levels") +
             ylab("Methylation Levels") +
             geom_smooth(method = lm, color="blue", se=FALSE) +
             theme_classic(base_size = 16)
  newPlot <- ggplot(data = as.data.frame(cor_data1), aes(x=RNASEQ_NormCnt_Perio, y=EPIC_Mvalue_Perio)) +
             geom_point(color = "black", size = 1) +
             labs(title=paste0(substr(fName, 1, nchar(fName)-11), "New"),
                  subtitle=sprintf("P.Cor = %s, p-value = %s",
                                   round(cor(cor_data1[,"RNASEQ_NormCnt_Perio"], cor_data1[,"EPIC_Mvalue_Perio"], use = "pairwise.complete.obs"), 5),
                                   signif(cor.test(cor_data1[,"RNASEQ_NormCnt_Perio"], cor_data1[,"EPIC_Mvalue_Perio"])$p.value, 5))) +
             xlab("Expression Levels") +
             ylab("Methylation Levels") +
             geom_smooth(method = lm, color="blue", se=FALSE) +
             theme_classic(base_size = 16)
  ggsave(filename = paste0(outputDir, fName), arrangeGrob(oldPlot, newPlot, ncol = 2), width = 20, height = 10)
  
}
