###
#   File name : LevelCorrelations.R
#   Author    : Hyunjin Kim
#   Date      : Nov 5, 2018
#   Email     : hk2990@cumc.columbia.edu
#   Purpose   : Make a correlation plot based on levels between methylation and gene expression
#
#   When matching CpGs to genes, there may be many CpGs associated with one gene
#   In that case, choose one which has the most negative spearman correlation
#   Therefore, the comparison is done with the matching genes
#
#   Instruction
#               1. Source("LevelCorrelations.R")
#               2. Run the function "levelCor" - specify the input files (levels of methylation and gene expression) and output directory
#               3. A correlation plot will be generated in the output directory
#
#   Example
#               > source("The_directory_of_LevelCorrelations.R/LevelCorrelations.R")
#               > levelCor(methylLevelPath="./data/Methylation/beta_values_gene.txt",
#                          gexpLevelPath="./data/RNA_Seq/raw_counts.txt",
#                          sampleInfoPath="./data/RNA_Seq/sampleInfo.txt",
#                          usingMethylMix=FALSE,
#                          importantGenePath="./results/Combined/List_-log10(FDR)_DMR_Periodontitis_Healthy_all.xlsx",
#                          fdrThreshold=0.05,
#                          outputDir="./results/Cor/")
###

levelCor <- function(methylLevelPath="//isilon.c2b2.columbia.edu/ifs/archive/shares/bisr/Papapanou/Sept_2018/PreprocessedData/Methylation/beta_values_gene.txt",
                     gexpLevelPath="//isilon.c2b2.columbia.edu/ifs/archive/shares/bisr/Papapanou/Sept_2018/PreprocessedData/RNA_Seq/raw_counts.txt",
                     sampleInfoPath="//isilon.c2b2.columbia.edu/ifs/archive/shares/bisr/Papapanou/Sept_2018/PreprocessedData/RNA_Seq/sampleInfo.txt",
                     usingMethylMix=FALSE,
                     importantGenePath="//isilon.c2b2.columbia.edu/ifs/archive/shares/bisr/Papapanou/Sept_2018/Combined/List_-log10(FDR)_DMR_Periodontitis_Healthy_all.xlsx",
                     fdrThreshold=0.05,
                     outputDir="//isilon.c2b2.columbia.edu/ifs/archive/shares/bisr/Papapanou/Sept_2018/Combined/") {
  
  ### load data
  methylLev <- read.table(file = methylLevelPath, header = TRUE, sep = "\t", check.names = FALSE, row.names = 1)
  gexpLev <- read.table(file = gexpLevelPath, header = TRUE, sep = "\t", check.names = FALSE)
  sampleInfo <- read.table(file = sampleInfoPath, header = TRUE, sep = "\t", check.names = FALSE, row.names = 1)
  
  
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
  gexpLev <- normalizeRNASEQwithVST(gexpLev, filtering = FALSE)
  
  
  ### make two datasets having the same format
  common_genes <- intersect(rownames(methylLev), rownames(gexpLev))
  common_samples <- intersect(colnames(methylLev), colnames(gexpLev))
  gexpLev <- gexpLev[common_genes,common_samples]
  methylLev <- methylLev[common_genes,common_samples]
  
  
  ### separate samples into three periontal categories
  sampleInfo <- sampleInfo[common_samples,]
  gexp_periodontitis <- gexpLev[,which(sampleInfo$Phenotype == "Periodontitis")]
  gexp_gingivitis <- gexpLev[,which(sampleInfo$Phenotype == "Gingivitis")]
  gexp_healthy <- gexpLev[,which(sampleInfo$Phenotype == "Healthy")]
  methyl_periodontitis <- methylLev[,which(sampleInfo$Phenotype == "Periodontitis")]
  methyl_gingivitis <- methylLev[,which(sampleInfo$Phenotype == "Gingivitis")]
  methyl_healthy <- methylLev[,which(sampleInfo$Phenotype == "Healthy")]
  
  
  ### using MethylMix package
  if(usingMethylMix == TRUE) {
    
    ### load library
    if(!require(MethylMix)) {
      source("https://bioconductor.org/biocLite.R")
      biocLite("MethylMix")
      library(MethylMix)
    }
    
    
    ### MethylMix plots
    MethylMixResults <- MethylMix(as.matrix(methyl_periodontitis),
                                  as.matrix(gexp_periodontitis),
                                  as.matrix(methyl_healthy))
    plots <- MethylMix_PlotModel("KCNAB2", MethylMixResults,
                                 as.matrix(methyl_periodontitis),
                                 as.matrix(gexp_periodontitis),
                                 as.matrix(methyl_healthy))
    plots$MixtureModelPlot
    plots$CorrelationPlot
  
  ### using my own codes
  } else {
    
    ### load library
    if(!require(xlsx)) {
      install.packages("xlsx")
      library(xlsx)
    }
    if(!require(ggplot2)) {
      install.packages("ggplot2")
      library(ggplot2)
    }
    
    ### get important genes
    gList <- read.xlsx2(file = importantGenePath, sheetIndex = 1, stringsAsFactors = FALSE)
    gList[,3:4] <- apply(gList[,3:4], 2, as.numeric)
    gList <- gList$Gene_Name[intersect(which(gList$FDR_Expression < fdrThreshold),
                                       which(gList$FDR_Methylation < fdrThreshold))]
    gList <- unique(gList)
    
    ### iteratively performes generating correlation plots
    for(i in 1:length(gList)) {
      ### make a source data for a plot
      cor_data <- cbind(t(gexp_periodontitis[gList[i],]), t(methyl_periodontitis[gList[i],]),
                        t(gexp_gingivitis[gList[i],]), t(methyl_gingivitis[gList[i],]),
                        t(gexp_healthy[gList[i],]), t(methyl_healthy[gList[i],]))
      colnames(cor_data) <- c("GEXP_PT", "METHYL_PT",
                              "GEXP_GV", "METHYL_GV",
                              "GEXP_HT", "METHYL_HT")
      
      ### create a plot for periodontitis samples
      fName <- paste0(gList[i], "_Periodontitis_Correlation.png")
      ggplot(data = as.data.frame(cor_data), aes(x=GEXP_PT, y=METHYL_PT)) +
        geom_point(color = "black", size = 1) +
        labs(title=substr(fName, 1, nchar(fName)-4),
             subtitle=sprintf("P.Cor = %s, p-value = %s",
                              round(cor(cor_data[,"GEXP_PT"], cor_data[,"METHYL_PT"], use = "pairwise.complete.obs"), 5),
                              signif(cor.test(cor_data[,"GEXP_PT"], cor_data[,"METHYL_PT"])$p.value, 5))) +
        xlab("Expression Levels") +
        ylab("Methylation Levels") +
        geom_smooth(method = lm, color="blue", se=FALSE) +
        theme_classic(base_size = 16)
      ggsave(filename = paste0(outputDir, fName), width = 10, height = 10)
      
      ### create a plot for gingivitis samples
      fName <- paste0(gList[i], "_Gingivitis_Correlation.png")
      ggplot(data = as.data.frame(cor_data), aes(x=GEXP_GV, y=METHYL_GV)) +
        geom_point(color = "black", size = 1) +
        labs(title=substr(fName, 1, nchar(fName)-4),
             subtitle=sprintf("P.Cor = %s, p-value = %s",
                              round(cor(cor_data[,"GEXP_GV"], cor_data[,"METHYL_GV"], use = "pairwise.complete.obs"), 5),
                              signif(cor.test(cor_data[,"GEXP_GV"], cor_data[,"METHYL_GV"])$p.value, 5))) +
        xlab("Expression Levels") +
        ylab("Methylation Levels") +
        geom_smooth(method = lm, color="blue", se=FALSE) +
        theme_classic(base_size = 16)
      ggsave(filename = paste0(outputDir, fName), width = 10, height = 10)
      
      ### create a plot for healthy samples
      fName <- paste0(gList[i], "_Healthy_Correlation.png")
      ggplot(data = as.data.frame(cor_data), aes(x=GEXP_PT, y=METHYL_PT)) +
        geom_point(color = "black", size = 1) +
        labs(title=substr(fName, 1, nchar(fName)-4),
             subtitle=sprintf("P.Cor = %s, p-value = %s",
                              round(cor(cor_data[,"GEXP_HT"], cor_data[,"METHYL_HT"], use = "pairwise.complete.obs"), 5),
                              signif(cor.test(cor_data[,"GEXP_HT"], cor_data[,"METHYL_HT"])$p.value, 5))) +
        xlab("Expression Levels") +
        ylab("Methylation Levels") +
        geom_smooth(method = lm, color="blue", se=FALSE) +
        theme_classic(base_size = 16)
      ggsave(filename = paste0(outputDir, fName), width = 10, height = 10)
      
    }
    
    ### correlation of all the genes - for a reference
    cor_data <- NULL
    for(i in 1:nrow(gexp_periodontitis)) {
      ###
      cor_data <- rbind(cor_data, cbind(t(gexp_periodontitis[i,]), t(methyl_periodontitis[i,]),
                                        t(gexp_gingivitis[i,]), t(methyl_gingivitis[i,]),
                                        t(gexp_healthy[i,]), t(methyl_healthy[i,])))
    }
    colnames(cor_data) <- c("GEXP_PT", "METHYL_PT",
                            "GEXP_GV", "METHYL_GV",
                            "GEXP_HT", "METHYL_HT")
    
    ### periodontitis
    fName <- "All_Genes_Periodontitis_Correlation.png"
    ggplot(data = as.data.frame(cor_data), aes(x=GEXP_PT, y=METHYL_PT)) +
      geom_point(color = "black", size = 1) +
      labs(title=substr(fName, 1, nchar(fName)-4),
           subtitle=sprintf("P.Cor = %s, p-value = %s",
                            round(cor(cor_data[,"GEXP_PT"], cor_data[,"METHYL_PT"], use = "pairwise.complete.obs"), 5),
                            signif(cor.test(cor_data[,"GEXP_PT"], cor_data[,"METHYL_PT"])$p.value, 5))) +
      xlab("Expression Levels") +
      ylab("Methylation Levels") +
      geom_smooth(method = lm, color="blue", se=FALSE) +
      theme_classic(base_size = 16)
    ggsave(filename = paste0(outputDir, fName), width = 10, height = 10)
    
    ### gingivitis
    fName <- "All_Genes_Gingivitis_Correlation.png"
    ggplot(data = as.data.frame(cor_data), aes(x=GEXP_GV, y=METHYL_GV)) +
      geom_point(color = "black", size = 1) +
      labs(title=substr(fName, 1, nchar(fName)-4),
           subtitle=sprintf("P.Cor = %s, p-value = %s",
                            round(cor(cor_data[,"GEXP_GV"], cor_data[,"METHYL_GV"], use = "pairwise.complete.obs"), 5),
                            signif(cor.test(cor_data[,"GEXP_GV"], cor_data[,"METHYL_GV"])$p.value, 5))) +
      xlab("Expression Levels") +
      ylab("Methylation Levels") +
      geom_smooth(method = lm, color="blue", se=FALSE) +
      theme_classic(base_size = 16)
    ggsave(filename = paste0(outputDir, fName), width = 10, height = 10)
    
    ### healthy
    fName <- "All_Genes_Healthy_Correlation.png"
    ggplot(data = as.data.frame(cor_data), aes(x=GEXP_HT, y=METHYL_HT)) +
      geom_point(color = "black", size = 1) +
      labs(title=substr(fName, 1, nchar(fName)-4),
           subtitle=sprintf("P.Cor = %s, p-value = %s",
                            round(cor(cor_data[,"GEXP_HT"], cor_data[,"METHYL_HT"], use = "pairwise.complete.obs"), 5),
                            signif(cor.test(cor_data[,"GEXP_HT"], cor_data[,"METHYL_HT"])$p.value, 5))) +
      xlab("Expression Levels") +
      ylab("Methylation Levels") +
      geom_smooth(method = lm, color="blue", se=FALSE) +
      theme_classic(base_size = 16)
    ggsave(filename = paste0(outputDir, fName), width = 10, height = 10)
    
  }
  
}
