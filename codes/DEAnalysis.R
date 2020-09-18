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
  
  #####################################################
  ### A function to perform DE analysis with DESeq2 ###
  #####################################################
  #' @title deseqWithComparisons
  #' @param rCnt raw count matrix (rows: genes, columns: samples)
  #' @param grp a character vector of class info of the samples
  #'            or a list with two character vectors representing
  #'            two different independent classes of the samples
  #' @param exp_class a string of the experiment group's name
  #'                  or a vector with two strings that are from
  #'                  the grp[[1]] for a comparison
  #' @param ctrl_class a string of the control group's name
  #'                  or a vector with two strings that are from
  #'                  the grp[[2]] for a comparison
  #'                  
  #' * If exp_class = "A" and ctrl_class = "B", the direction of
  #'   differential expression is A - B. If exp_class = c("A", "B") and
  #'   and ctrl_class = c("C", "D"), then the direction of
  #'   differential expression is (AC - AD) - (BC - BD).
  #'   
  #' @param bat_eff a character vector of batch effect info of the samples
  #' @param thresh numeric. Filters out from the results genes with adjusted
  #' 		           p-value larger than this value
  #' @return data.frame
  #' @export
  #' @author Hyunjin Kim
  ####################################################
  deseqWithComparisons <- function(rCnt, grp, exp_class, ctrl_class, bat_eff=NULL, thresh = 1) {
    
    ### load library
    if(!require(DESeq2, quietly = TRUE)) {
      source("https://bioconductor.org/biocLite.R")
      biocLite("DESeq2")
      require(DESeq2, quietly = TRUE)
    }
    if(!require("checkmate", quietly = TRUE)) {
      install.packages("checkmate")
      require("checkmate", quietly = TRUE)
    }
    
    ### argument checking
    assertDataFrame(rCnt)
    assert(checkVector(grp), checkList(grp))
    if(testList(grp)) {
      assert(length(grp) == 2)
      assert(length(exp_class) == 2)
      assert(length(ctrl_class) == 2)
      assert(length(which(grp[[1]] == exp_class[1])) > 0)
      assert(length(which(grp[[1]] == exp_class[2])) > 0)
      assert(length(which(grp[[2]] == ctrl_class[1])) > 0)
      assert(length(which(grp[[2]] == ctrl_class[2])) > 0)
    } else {
      assert(length(exp_class) == 1)
      assert(length(ctrl_class) == 1)
      assert(length(which(grp == exp_class)) > 0, length(which(grp == ctrl_class)) > 0)
    }
    assertCharacter(exp_class)
    assertCharacter(ctrl_class)
    assert(checkNull(bat_eff), checkCharacter(bat_eff), checkFactor(bat_eff))
    assertNumeric(thresh)
    
    ### sometimes, there are some variables which can not be transformed into R variable names
    ### so, just to be safe, change all the variables to R-usuable ones
    exp_class <- make.names(exp_class)
    ctrl_class <- make.names(ctrl_class)
    
    ### there are two different grp options
    if(testList(grp)) {
      sampleType1 <- make.names(as.character(grp[[1]]))
      sampleType2 <- make.names(as.character(grp[[2]]))
      
      ### there are two options: considering batch effect or not
      if(is.null(bat_eff)) {
        ### make a data frame for design matrix
        Coldata <- data.frame(sampleType1, sampleType2)
        
        ### set names and the refererence group
        rownames(Coldata) <- colnames(rCnt)
        Coldata$sampleType1 <- relevel(Coldata$sampleType1, ref = exp_class[2])
        Coldata$sampleType2 <- relevel(Coldata$sampleType2, ref = ctrl_class[2])
        
        ### data preparation for DE analysis
        deSeqData <- DESeqDataSetFromMatrix(countData=rCnt, colData=Coldata, design= ~sampleType1*sampleType2)
      } else {
        ### make a data frame for design matrix
        batch_eff <- as.character(bat_eff)
        Coldata <- data.frame(sampleType1, sampleType2, batch_eff)
        
        ### set names and the refererence group
        rownames(Coldata) <- colnames(rCnt)
        Coldata$sampleType1 <- relevel(Coldata$sampleType1, ref = exp_class[2])
        Coldata$sampleType2 <- relevel(Coldata$sampleType2, ref = ctrl_class[2])
        
        ### data preparation for DE analysis
        deSeqData <- DESeqDataSetFromMatrix(countData=rCnt, colData=Coldata, design= ~sampleType1*sampleType2+batch_eff)
      }
      
      ### filtering some useless genes out
      deSeqData <- deSeqData[rowSums(counts(deSeqData))>1,]
      
      ### run DE analysis
      dea <- DESeq(deSeqData)
      
      ### get DE genes using the contrast
      deresults <- results(dea, name = paste0("sampleType1", exp_class[1], ".sampleType2", ctrl_class[1]))
    } else {
      sampleType <- make.names(as.character(grp))
      
      ### there are two options: considering batch effect or not
      if(is.null(bat_eff)) {
        ### make a data frame for design matrix
        Coldata <- data.frame(sampleType)
        
        ### set names and the refererence group
        rownames(Coldata) <- colnames(rCnt)
        Coldata$sampleType <- relevel(Coldata$sampleType, ref = ctrl_class)
        
        ### data preparation for DE analysis
        deSeqData <- DESeqDataSetFromMatrix(countData=rCnt, colData=Coldata, design= ~sampleType)
      } else {
        ### make a data frame for design matrix
        batch_eff <- as.character(bat_eff)
        Coldata <- data.frame(sampleType, batch_eff)
        
        ### set names and the refererence group
        rownames(Coldata) <- colnames(rCnt)
        Coldata$sampleType <- relevel(Coldata$sampleType, ref = ctrl_class)
        
        ### data preparation for DE analysis
        deSeqData <- DESeqDataSetFromMatrix(countData=rCnt, colData=Coldata, design= ~sampleType+batch_eff)
      }
      
      ### filtering some useless genes out
      deSeqData <- deSeqData[rowSums(counts(deSeqData))>1,]
      
      ### run DE analysis
      dea <- DESeq(deSeqData)
      
      ### get DE genes using the contrast
      deresults <- results(dea, contrast = c("sampleType", exp_class, ctrl_class))
    }
    
    ### change p-values that have NA values to 1
    if(length(which(is.na(deresults$pvalue))) > 0) {
      deresults$pvalue[which(is.na(deresults$pvalue))] <- 1
    }
    if(length(which(is.na(deresults$padj))) > 0) {
      deresults$padj[which(is.na(deresults$padj))] <- 1
    }
    
    ### ordering and filtering the DE result
    deresults <- deresults[order(deresults$padj, na.last = TRUE), ,drop = FALSE]
    deresults <- deresults[deresults$padj <= thresh, ,drop = FALSE]
    
    ### there are two different grp options
    if(testList(grp)) {
      ### add baseMean for each group
      nCnt <- counts(dea, normalized=TRUE, replaced=TRUE)
      exp1_ctrl1_rowMeans <- apply(nCnt[,intersect(which(Coldata$sampleType1 == exp_class[1]),
                                                   which(Coldata$sampleType2 == ctrl_class[1])),
                                        drop=FALSE], 1, mean)
      exp1_ctrl2_rowMeans <- apply(nCnt[,intersect(which(Coldata$sampleType1 == exp_class[1]),
                                                   which(Coldata$sampleType2 == ctrl_class[2])),
                                        drop=FALSE], 1, mean)
      exp2_ctrl1_rowMeans <- apply(nCnt[,intersect(which(Coldata$sampleType1 == exp_class[2]),
                                                   which(Coldata$sampleType2 == ctrl_class[1])),
                                        drop=FALSE], 1, mean)
      exp2_ctrl2_rowMeans <- apply(nCnt[,intersect(which(Coldata$sampleType1 == exp_class[2]),
                                                   which(Coldata$sampleType2 == ctrl_class[2])),
                                        drop=FALSE], 1, mean)
      deresults <- data.frame(baseMean=deresults[,1],
                              V1=exp1_ctrl1_rowMeans[rownames(deresults)],
                              V2=exp1_ctrl2_rowMeans[rownames(deresults)],
                              V3=exp2_ctrl1_rowMeans[rownames(deresults)],
                              V4=exp2_ctrl2_rowMeans[rownames(deresults)],
                              deresults[,2:6],
                              stringsAsFactors = FALSE, check.names = FALSE)
      colnames(deresults)[2:5] <- c(paste0("baseMean_", exp_class[1], "_", ctrl_class[1]),
                                    paste0("baseMean_", exp_class[1], "_", ctrl_class[2]),
                                    paste0("baseMean_", exp_class[2], "_", ctrl_class[1]),
                                    paste0("baseMean_", exp_class[2], "_", ctrl_class[2]))
    } else {
      ### add baseMean for each group
      nCnt <- counts(dea, normalized=TRUE, replaced=TRUE)
      exp_rowMeans <- apply(nCnt[,which(Coldata$sampleType == exp_class), drop=FALSE], 1, mean)
      ctrl_rowMeans <- apply(nCnt[,which(Coldata$sampleType == ctrl_class), drop=FALSE], 1, mean)
      deresults <- data.frame(baseMean=deresults[,1],
                              V1=exp_rowMeans[rownames(deresults)],
                              V2=ctrl_rowMeans[rownames(deresults)],
                              deresults[,2:6],
                              stringsAsFactors = FALSE, check.names = FALSE)
      colnames(deresults)[2:3] <- c(paste0("baseMean_", exp_class), paste0("baseMean_", ctrl_class))  
    }
    
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





