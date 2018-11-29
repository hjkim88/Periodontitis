###
#   File name : CpGToGenes.R
#   Author    : Hyunjin Kim
#   Date      : Oct 25, 2018
#   Email     : hk2990@cumc.columbia.edu
#   Purpose   : Change rows of methylation data to genes
#
#   There may be many CpGs associated with one gene
#   In that case, choose one which has the most negative spearman correlation 
#
#   Instruction
#               1. Source("CpGToGenes.R")
#               2. Run the function "changeToGenes" - specify the input file (beta/m) of methylation and output file path
#               3. Gene version of methylation data will be generated in the output file path
#
#   Example
#               > source("The_directory_of_CpGToGenes.R/CpGToGenes.R")
#               > changeToGenes(inputPath="./data/Methylation/beta_values.txt",
#                               rCntPath="./data/RNA_Seq/raw_counts.txt",
#                               outputPath="./data/Methylation/beta_values_gene.txt")
###

changeToGenes <- function(inputPath="./data/Methylation/beta_values.txt",
                          rCntPath="./data/RNA_Seq/raw_counts.txt",
                          outputPath="./data/Methylation/beta_values_gene.txt") {
  
  ### load library
  if(!require(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)) {
    source("https://bioconductor.org/biocLite.R")
    biocLite("IlluminaHumanMethylationEPICanno.ilm10b2.hg19")
    library(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)
  }
  
  
  ### load data
  methylData <- read.table(file = inputPath, header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)
  rCnt <- read.table(file = rCntPath, header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)
  
  
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
  normCnt <- normalizeRNASEQwithVST(rCnt)
  
  
  ### get Illumina EPIC hg19 annotation data
  anno <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)
  anno <- anno[match(rownames(methylData), anno$Name), "GencodeCompV12_NAME"]
  names(anno) <- rownames(methylData)
  
  
  ### get unique gene set
  unique_gene_set <- NULL
  unique_anno <- unique(anno)
  for(i in 1:length(unique_anno)) {
    temp <- strsplit(unique_anno[i], ";", fixed = TRUE)[[1]]
    unique_gene_set <- unique(c(unique_gene_set, temp))
    
    if(i %% 1000 == 0) {
      writeLines(paste(i, "/", length(unique_anno)))
    }
  }
  unique_gene_set <- intersect(unique_gene_set, rownames(normCnt))
  
  
  ### make an empty result matrix
  methyl_gene <- matrix(NA, length(unique_gene_set), ncol(methylData))
  rownames(methyl_gene) <- unique_gene_set
  colnames(methyl_gene) <- colnames(methylData)
  
  
  ### get a CpG with the most negative cor value with gene expression
  for(i in 1:length(unique_gene_set)) {
    ### get indicies of associated CpGs
    idx <- grep(unique_gene_set[i], anno)
    
    ### calculate correlations
    corr <- NULL
    for(j in 1:length(idx)) {
      corr <- c(corr, cor(as.numeric(methylData[idx[j],]),
                          as.numeric(normCnt[unique_gene_set[i],]),
                          use = "pairwise.complete.obs",
                          method = "spearman"))
    }
    
    ### index of the most negative correlation
    idx <- idx[which(corr == min(corr))]
    
    ### tie breaker with pearson correlation
    if(length(idx) > 1) {
      corr <- NULL
      for(j in 1:length(idx)) {
        corr <- c(corr, cor(as.numeric(methylData[idx[j],]),
                            as.numeric(normCnt[unique_gene_set[i],]),
                            use = "pairwise.complete.obs",
                            method = "pearson"))
      }
      idx <- idx[which(corr == min(corr))][1]
    }
    
    ### put methylation levels
    methyl_gene[i,] <- as.numeric(methylData[idx,])
    
    ### print progress
    if(i %% 500 == 0) {
      writeLines(paste(i, "/", length(unique_gene_set)))
    }
  }
  
  
  ### write out the result
  methyl_gene <- data.frame(Gene_Symbol=rownames(methyl_gene), methyl_gene)
  write.table(methyl_gene, file = outputPath, sep = "\t", row.names = FALSE)
  
}


