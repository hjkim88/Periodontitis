###
#   File name : DEAnalysis.R
#   Author    : Hyunjin Kim
#   Date      : Oct 14, 2018
#   Email     : hk2990@cumc.columbia.edu
#   Purpose   : Perfrom DE analysis using DESeq2 on Panos' RNA-Seq data
#
#   Instruction
#               1. Source("DEAnalysis.R")
#               2. Run the function "DEA_Panos" - specify an input file path (raw counts), a comparison list, and output directory
#               3. The differential expression analysis results will be generated under the output directory
#
#   Example
#               > source("The_directory_of_DEAnalysis.R/DEAnalysis.R")
#               > DEA_Panos(rawCntPath="./data/RNA_Seq/panos_rna_seq_data.rda",
#                           comparisonList=list(c("Gingivitis", "Healthy"),
#                                               c("Periodontitis", "Healthy"),
#                                               c("Periodontitis", "Gingivitis")),
#                           fdrThreshold=0.05,
#                           outputDir="./results/DEA/")
###

DEA_Panos <- function(rawCntPath="//isilon.c2b2.columbia.edu/ifs/archive/shares/bisr/Papapanou/Sept_2018/PreprocessedData/RNA_Seq/panos_rna_seq_data.rda",
                      comparisonList=list(c("Gingivitis", "Healthy"),
                                          c("Periodontitis", "Healthy"),
                                          c("Periodontitis", "Gingivitis")),
                      fdrThreshold=0.05,
                      outputDir="//isilon.c2b2.columbia.edu/ifs/archive/shares/bisr/Papapanou/Sept_2018/DifferentialExpression/") {
  
  ### load library
  if(!require(xlsx)) {
    install.packages("xlsx")
    library(xlsx)
  }
  
  ### load dataset
  load(rawCntPath)
  
  ### A function to perform repetitive DE analysis with DESeq2
  deseqWithComparisons <- function(rCnt, grp, exp_class, ctrl_class, bat_eff=NULL) {
    
    ### load library
    if(!require(DESeq2)) {
      source("https://bioconductor.org/biocLite.R")
      biocLite("DESeq2")
      library(DESeq2)
    }
    
    ### make a design matrix for DE analysis
    sampleType <- as.character(grp)
    
    if(is.null(bat_eff)) {
      Coldata <- data.frame(sampleType)
    } else {
      batch_eff <- as.character(bat_eff)
      Coldata <- data.frame(sampleType, batch_eff)
    }
    
    rownames(Coldata) <- colnames(rCnt)
    Coldata$sampleType <- relevel(Coldata$sampleType, ref = ctrl_class)
    
    ### data preparation for DE analysis
    if(is.null(bat_eff)) {
      deSeqData <- DESeqDataSetFromMatrix(countData=rCnt, colData=Coldata, design= ~sampleType)
    } else {
      deSeqData <- DESeqDataSetFromMatrix(countData=rCnt, colData=Coldata, design= ~sampleType+batch_eff)
    }
    
    deSeqData <- deSeqData[rowSums(counts(deSeqData))>1,]
    
    ### run DE analysis
    dea <- DESeq(deSeqData)
    deresults <- results(dea, contrast = c("sampleType", exp_class, ctrl_class))
    if(length(which(is.na(deresults$pvalue))) > 0) {
      deresults$pvalue[which(is.na(deresults$pvalue))] <- 1
    }
    if(length(which(is.na(deresults$padj))) > 0) {
      deresults$padj[which(is.na(deresults$padj))] <- 1
    }
    deresults <- deresults[order(deresults$padj, na.last = TRUE),]
    
    return(deresults)
  }
  
  
  ### A function to print volcano plot of DE analysis with DESeq2 result
  volPlotWithDeseq <- function(deresult, outputFilePath, pvalue=0.05) {
    
    ### load library
    if(!require(ggplot2)) {
      install.packages("ggplot2")
      library(ggplot2)
    }
    
    deresult$padj[which(is.na(deresult$padj))] <- 1
    volcanoData <- as.data.frame(cbind(deresult$log2FoldChange, -log10(deresult$padj), as.character(as.factor(deresult$padj < pvalue))))
    colnames(volcanoData) <- c("logFC", "logFDR", "Significance")
    volcanoData$logFC <- as.numeric(as.character(volcanoData$logFC))
    volcanoData$logFDR <- as.numeric(as.character(volcanoData$logFDR))
    
    s <- strsplit(outputFilePath, "/")
    f1 <- s[[1]][length(s[[1]])]
    
    ggplot(data=volcanoData, aes(x=logFC, y=logFDR, colour=Significance)) + xlab("log2_Fold_Change") + ylab("-log10(FDR)") + ggtitle(paste("Significant (adj.p < ", pvalue, " ) DE genes -", sum(volcanoData$Significance == "TRUE"))) + theme_classic(base_size = 16) + geom_point(alpha=0.4)
    ggsave(filename = outputFilePath)
  }
  
  ### iteratively perform DE analysis for every comparison in the list
  for(i in 1:length(comparisonList)) {
    ### DE analysis with DESeq2
    de <- deseqWithComparisons(rCnt = rawCnt,
                               grp = sampleInfo$Phenotype,
                               exp_class = comparisonList[[i]][1],
                               ctrl_class = comparisonList[[i]][2],
                               bat_eff = sampleInfo$Batch)
    
    ### set file name
    fileName <- paste0("DE_Results_", comparisonList[[i]][1], "_vs_", comparisonList[[i]][2])
    
    ### put gene symbol column
    de <- cbind(Gene_Symbol=rownames(de), data.frame(de))
    
    ### numerize the factor columns
    idx <- sapply(de, is.factor)
    de[idx] <- lapply(de[idx], function(x) as.character(x))
    
    ### write out the DE result
    write.xlsx2(de, file = paste0(outputDir, fileName, ".xlsx"), sheetName = fileName, row.names = FALSE)
    
    # Volcano plot
    volPlotWithDeseq(de, paste0(outputDir, fileName, "_volPlot.png"), pvalue = fdrThreshold)
  }
  
}





