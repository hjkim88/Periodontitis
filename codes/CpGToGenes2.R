###
#   File name : CpGToGenes2.R
#   Author    : Hyunjin Kim
#   Date      : Nov 19, 2018
#   Email     : hk2990@cumc.columbia.edu
#   Purpose   : Change rows of methylation data to genes
#
#   There may be many CpGs associated with one gene
#   In that case, choose one which has the most negative spearman correlation 
#
#   This code is different from "CpGToGenes.R", since this is for
#   Panos' old datasets (450k and Affy chip)
#
#   Instruction
#               1. Source("CpGToGenes2.R")
#               2. Run the function "changeToGenes2" - specify the input file (beta/m) of methylation and output file path
#               3. Gene version of methylation data will be generated in the output file path
#
#   Example
#               > source("The_directory_of_CpGToGenes2.R/CpGToGenes2.R")
#               > changeToGenes2(inputPath="./data/450k/panos_450k_data.rda",
#                                rCntPath="./data/Affy/panos_affy_data.rda",
#                                outputPath="./data/450k/450k_m_values_gene.txt")
###

changeToGenes2 <- function(inputPath="//isilon.c2b2.columbia.edu/ifs/archive/shares/bisr/Papapanou/Sept_2018/PreprocessedData/450k/panos_450k_data.rda",
                           rCntPath="//isilon.c2b2.columbia.edu/ifs/archive/shares/bisr/Papapanou/Sept_2018/PreprocessedData/Affy/panos_affy_data.rda",
                           outputPath="//isilon.c2b2.columbia.edu/ifs/archive/shares/bisr/Papapanou/Sept_2018/PreprocessedData/450k/450k_m_values_gene.txt") {
  
  ### load library
  if(!require(IlluminaHumanMethylation450kanno.ilmn12.hg19)) {
    source("https://bioconductor.org/biocLite.R")
    biocLite("IlluminaHumanMethylation450kanno.ilmn12.hg19")
    library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
  }
  
  
  ### load data
  load(inputPath)
  load(rCntPath)
  
  
  ### get Illumina 450k hg19 annotation data
  anno <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
  anno <- anno[match(rownames(mvalues_450k), anno$Name), "UCSC_RefGene_Name"]
  names(anno) <- rownames(mvalues_450k)
  
  
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
  unique_gene_set <- intersect(unique_gene_set, unique(affy_norm_ge$Gene_Symbol))
  common_affy_probe_set <- affy_norm_ge$Gene_Symbol[which(affy_norm_ge$Gene_Symbol %in% unique_gene_set)]
  names(common_affy_probe_set) <- rownames(affy_norm_ge)[which(affy_norm_ge$Gene_Symbol %in% unique_gene_set)]
  
  
  ### make an empty result matrix
  methyl_gene <- matrix(NA, length(common_affy_probe_set), ncol(mvalues_450k))
  rownames(methyl_gene) <- names(common_affy_probe_set)
  colnames(methyl_gene) <- colnames(mvalues_450k)
  
  
  ### get a CpG with the most negative cor value with gene expression
  common_samples <- intersect(colnames(affy_norm_ge), colnames(mvalues_450k))
  for(i in 1:length(common_affy_probe_set)) {
    ### get indicies of associated CpGs
    idx <- grep(common_affy_probe_set[i], anno)
    
    ### calculate correlations
    corr <- NULL
    for(j in 1:length(idx)) {
      corr <- c(corr, cor(as.numeric(mvalues_450k[idx[j],common_samples]),
                          as.numeric(affy_norm_ge[names(common_affy_probe_set)[i],common_samples]),
                          use = "pairwise.complete.obs",
                          method = "spearman"))
    }
    
    ### index of the most negative correlation
    idx <- idx[which(corr == min(corr))]
    
    ### tie breaker with pearson correlation
    if(length(idx) > 1) {
      corr <- NULL
      for(j in 1:length(idx)) {
        corr <- c(corr, cor(as.numeric(mvalues_450k[idx[j],common_samples]),
                            as.numeric(affy_norm_ge[names(common_affy_probe_set)[i],common_samples]),
                            use = "pairwise.complete.obs",
                            method = "pearson"))
      }
      idx <- idx[which(corr == min(corr))][1]
    }
    
    ### put methylation levels
    methyl_gene[i,] <- as.numeric(mvalues_450k[idx,])
    
    ### print progress
    if(i %% 500 == 0) {
      writeLines(paste(i, "/", length(common_affy_probe_set)))
    }
  }
  
  
  ### write out the result
  methyl_gene <- data.frame(Probe_ID=names(common_affy_probe_set),
                            Gene_Symbol=common_affy_probe_set,
                            methyl_gene,
                            check.names = FALSE, stringsAsFactors = FALSE)
  write.table(methyl_gene, file = outputPath, sep = "\t", row.names = FALSE)
  
}


