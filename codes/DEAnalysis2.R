###
#   File name : DEAnalysis2.R
#   Author    : Hyunjin Kim
#   Date      : Nov 23, 2018
#   Email     : hk2990@cumc.columbia.edu
#   Purpose   : Perfrom DE analysis using DESeq2 on Panos' Affy chip data
#
#   * THIS CODE IS DIFFERENT FROM "DEAnalysis.R" since this is for the Affy chip data,
#     and the other one is for the RNA-Seq data of Panos.
#
#   Instruction
#               1. Source("DEAnalysis2.R")
#               2. Run the function "DEA_Panos2" - specify an input file path (normalized gexp) and output directory
#               3. The differential expression analysis results will be generated under the output directory
#
#   Example
#               > source("The_directory_of_DEAnalysis2.R/DEAnalysis2.R")
#               > DEA_Panos2(normGEXPPath="./data/Affy/panos_affy_data.rda",
#                            fdrThreshold=0.05,
#                            outputDir="./results/DEA/")
###

DEA_Panos2 <- function(normGEXPPath="//isilon.c2b2.columbia.edu/ifs/archive/shares/bisr/Papapanou/Sept_2018/PreprocessedData/Affy/panos_affy_data.rda",
                       fdrThreshold=0.05,
                       outputDir="//isilon.c2b2.columbia.edu/ifs/archive/shares/bisr/Papapanou/Sept_2018/DifferentialExpression/") {
  
  ### load library
  if(!require(xlsx)) {
    install.packages("xlsx")
    library(xlsx)
  }
  
  
  ### load data
  load(normGEXPPath)
  
  
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
  
  
  ### A function to print volcano plot of DE analysis with limma result
  volPlotWithLimma <- function(deresult, outputFilePath, pvalue=0.05) {
    
    ### load library
    if(!require(ggplot2)) {
      install.packages("ggplot2")
      library(ggplot2)
    }
    
    deresult$adj.P.Val[which(is.na(deresult$adj.P.Val))] <- 1
    volcanoData <- as.data.frame(cbind(deresult$logFC, -log10(deresult$adj.P.Val), as.character(as.factor(deresult$adj.P.Val < pvalue))))
    colnames(volcanoData) <- c("logFC", "logFDR", "Significance")
    volcanoData$logFC <- as.numeric(as.character(volcanoData$logFC))
    volcanoData$logFDR <- as.numeric(as.character(volcanoData$logFDR))
    
    s <- strsplit(outputFilePath, "/")
    f1 <- s[[1]][length(s[[1]])]
    
    ggplot(data=volcanoData, aes(x=logFC, y=logFDR, colour=Significance)) + xlab("log2_Fold_Change") + ylab("-log10(FDR)") + ggtitle(paste("Significant (adj.p < ", pvalue, " ) DE genes -", sum(volcanoData$Significance == "TRUE"))) + theme_classic(base_size = 16) + geom_point(alpha=0.4)
    ggsave(filename = outputFilePath)
  }
  
  
  ### DE analysis with limma
  deresult <- limmaWithComparisons(normCnt = affy_norm_ge[,-1],
                                   grp = affy_norm_ge_sample_info$Condition,
                                   exp_class = "Affected",
                                   ctrl_class = "Control",
                                   bat_eff = affy_norm_ge_sample_info$Patient)
  
  
  ### annotate gene symbols
  deresult <- data.frame(Gene_Symbol=affy_norm_ge[rownames(deresult), "Gene_Symbol"],
                         deresult, stringsAsFactors = FALSE)
  
  
  ### write out the DE results in Excel
  write.xlsx2(deresult, file = paste0(outputDir, "limma_affy_periodontitis_vs_healthy.xlsx"))
  
  
  ### volcano plot
  volPlotWithLimma(deresult = deresult,
                   outputFilePath = paste0(outputDir, "volPlot_affy_periodontitis_vs_healthy.png"),
                   pvalue = 1e-20)
  
}
