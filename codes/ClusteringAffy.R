###
#   File name : ClusteringAffy.R
#   Author    : Hyunjin Kim
#   Date      : Dec 3, 2018
#   Email     : hk2990@cumc.columbia.edu
#   Purpose   : Because the Affy data has too many DE genes, we want to know if the periodontitis samples
#               carry more features than just periodontal disease condition. Therefore, clustering
#               the Affy samples with RNA-Seq samples to see if they are clustered in a group.
#
#   Instruction
#               1. Source("ClusteringAffy.R")
#               2. Run the function "clusteringAffy" - specify the input files (gene expression of the Affy and the RNA-Seq) and output directory
#               3. The clustering results will be generated in the output directory
#
#   Example
#               > source("The_directory_of_ClusteringAffy.R/ClusteringAffy.R")
#               > clusteringAffy()
###

clusteringAffy <- function(affyPath="./data/Affy/panos_affy_data.rda",
                           rnaseqPath="./data/RNA_Seq/panos_rna_seq_data.rda",
                           outputDir="./results/Test/") {
  
  ### load library
  if(!require(sva, quietly = TRUE)) {
    source("https://bioconductor.org/biocLite.R")
    biocLite("sva")
    require(sva, quietly = TRUE)
  }
  
  
  ### load datasets
  load(affyPath)
  load(rnaseqPath)
  
  
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
  
  
  ### normalize the RNA-Seq raw counts
  normCnt <- normalizeRNASEQwithVST(rawCnt, filtering = TRUE)
  
  
  ### correct batch effect of RNA-Seq data
  normCnt <- ComBat(dat = normCnt, batch = sampleInfo$Batch, mod = model.matrix(~Phenotype, data = sampleInfo))
  
  
  ### get common genes
  common_genes <- intersect(affy_norm_ge$Gene_Symbol, rownames(normCnt))
  
  
  ### one gene - one row. If there are more than two probes mapped to one gene, sum up the values
  new_affy_norm_ge <- NULL
  for(i in 1:length(common_genes)) {
    temp <- affy_norm_ge[which(affy_norm_ge$Gene_Symbol == common_genes[i]), -1]
    new_affy_norm_ge <- rbind(new_affy_norm_ge, apply(temp, 2, sum))
  }
  rownames(new_affy_norm_ge) <- common_genes
  new_affy_norm_ge <- as.data.frame(new_affy_norm_ge)
  
  
  ### RNA-Seq normCnt with the common genes
  new_normCnt <- normCnt[common_genes,]
  
  
  ### combine affy & rna-seq data
  combined_data <- data.frame(new_affy_norm_ge, new_normCnt,
                              stringsAsFactors = FALSE, check.names = FALSE)
  combined_sample_info <- data.frame(Phenotype=c(affy_norm_ge_sample_info$Condition, sampleInfo$Phenotype),
                                     Chip=c(rep("Affy", ncol(new_affy_norm_ge)), rep("RNASeq", ncol(new_normCnt))),
                                     stringsAsFactors = FALSE, check.names = FALSE)
  combined_sample_info[which(combined_sample_info$Phenotype == "Affected"), "Phenotype"] <- "Periodontitis"
  combined_sample_info[which(combined_sample_info$Phenotype == "Control"), "Phenotype"] <- "Healthy"
  rownames(combined_sample_info) <- colnames(combined_data)
  
  
  ### remove gingivitis samples since we are not interseted in those samples here
  gIdx <- which(combined_sample_info$Phenotype == "Gingivitis")
  combined_data <- combined_data[,-gIdx]
  combined_sample_info <- combined_sample_info[-gIdx,]
  
  
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
  
  
  ### PCA plot - lebeling with the phenotype info
  pca_plot(normalizedMat = combined_data,
           grp = combined_sample_info$Phenotype,
           title = "pca_combined_affy_rnaseq_phenotype",
           filePath = paste0(outputDir, "pca_combined_affy_rnaseq_phenotype.png"))
  
  
  ### PCA plot - lebeling with the chip info
  pca_plot(normalizedMat = combined_data,
           grp = combined_sample_info$Chip,
           title = "pca_combined_affy_rnaseq_chip",
           filePath = paste0(outputDir, "pca_combined_affy_rnaseq_chip.png"))
  
  
  ### chip effect correction
  combined_data_cc <- ComBat(dat = combined_data, batch = combined_sample_info$Chip,
                             mod = model.matrix(~combined_sample_info$Phenotype))
  
  
  ### PCA plot - lebeling with the phenotype info - After chip effect correction
  pca_plot(normalizedMat = combined_data_cc,
           grp = combined_sample_info$Phenotype,
           title = "pca_combined_affy_rnaseq_phenotype_chip_effect_corrected",
           filePath = paste0(outputDir, "pca_combined_affy_rnaseq_phenotype_chip_effect_corrected.png"))
  
  
  ### PCA plot - lebeling with the chip info - After chip effect correction
  pca_plot(normalizedMat = combined_data_cc,
           grp = combined_sample_info$Chip,
           title = "pca_combined_affy_rnaseq_chip_chip_effect_corrected",
           filePath = paste0(outputDir, "pca_combined_affy_rnaseq_chip_chip_effect_corrected.png"))
  
}




