###
#   File name : JCP_Revision.R
#   Author    : Hyunjin Kim
#   Date      : Mar 24, 2021
#   Email     : hyunjin.kim@stjude.org
#   Purpose   : 1. The present study although is done very carefully, it would be interesting
#                  to get the author’s perspective on the influence of age, gender, race, and
#                  Mendelian inheritance of the study population. These variables may greatly
#                  influence differential DNA methylation/RNA expression datasets.
#                  Discussion regarding this would benefit the readers.
#               2. A pie chart and Venn diagram showing number of genes overlapping
#                  (methylated/unmethylated/differentially expressed) among health, gingivitis,
#                  and periodontitis group would be beneficial and easy to understand
#                  in addition to the tables provided.
#               3. Similarly, a bar graph representing functional/regulatory regions
#                  (promoter, enhancer, non-coding etc.) scanned would be beneficial.
#               4. Differentially methylated/expression regions on genes can be represented in
#                  genome browser view (Integrated genomics viewer). This is important not only for
#                  visualization purpose but also represents functional elements/enrichment
#                  (For example: Hypermethylated region vs. RNA expression level overlap)
#
#   Instruction
#               1. Source("JCP_Revision.R")
#               2. Run the function "jcr_revision" - specify the input files (Differential) of methylation and RNA-Seq and output directory
#               3. The revision results will be generated in the output directory
#
#   Example
#               > source("The_directory_of_JCP_Revision.R/JCP_Revision.R")
#               > jcr_revision(methylPath="./results/DMA/DMP_Periodontitis_Healthy.txt",
#                              methylRegPath="./results/DMA/DMR_Periodontitis_Healthy.txt",
#                              methylLevelPath="./data/Methylation/beta_values_gene.txt",
#                              gexpLevelPath="./data/RNA_Seq/raw_counts.txt",
#                              rnaseqPath="./results/DEA/DE_Results_Periodontitis_vs_Healthy.xlsx",
#                              sampleInfoPath="./data/RNA_Seq/sampleInfo.txt",
#                              outputDir="./results/Revision/")
###

jcr_revision <- function(methylPath="./results/DMA/DMP_Periodontitis_Healthy.txt",
                         methylRegPath="./results/DMA/DMR_Periodontitis_Healthy.txt",
                         rnaseqPath="./results/DEA/DE_Results_Periodontitis_vs_Healthy.xlsx",
                         methylLevelPath="./data/Methylation/beta_values_gene.txt",
                         gexpLevelPath="./data/RNA_Seq/raw_counts.txt",
                         sampleInfoPath="./data/RNA_Seq/sampleInfo.txt",
                         outputDir="./results/Revision/") {
  
  ### load libraries
  if(!require(xlsx)) {
    install.packages("xlsx")
    library(xlsx)
  }
  if(!require(ggplot2)) {
    install.packages("ggplot2")
    library(ggplot2)
  }
  if(!require(ggforce)) {
    install.packages("ggforce")
    library(ggforce)
  }
  if(!require(gridExtra, quietly = TRUE)) {
    install.packages("gridExtra")
    require(gridExtra, quietly = TRUE)
  }
  if(!require(VennDiagram)) {
    install.packages("VennDiagram")
    library(VennDiagram)
  }
  if(!require(RColorBrewer, quietly = TRUE)) {
    install.packages("RColorBrewer")
    require(RColorBrewer, quietly = TRUE)
  }
  
  ### load data
  dmp <- read.table(file = methylPath, header = TRUE, sep = "\t", check.names = FALSE)
  dmr <- read.table(file = methylRegPath, header = TRUE, sep = "\t", check.names = FALSE)
  de <- read.xlsx2(file = rnaseqPath, sheetIndex = 1, stringsAsFactors = FALSE)
  methylLev <- read.table(file = methylLevelPath, header = TRUE, sep = "\t", check.names = FALSE, row.names = 1)
  gexpLev <- read.table(file = gexpLevelPath, header = TRUE, sep = "\t", check.names = FALSE)
  sampleInfo <- read.table(file = sampleInfoPath, header = TRUE, sep = "\t", check.names = FALSE)
  
  ### set row names
  rownames(dmp) <- dmp$Name
  rownames(de) <- de$Gene_Symbol
  rownames(sampleInfo) <- sampleInfo$RNA.seq.ID
  
  ### change character columns to numeric
  dmp[,c(28, 30, 32, 36:41)] <- sapply(dmp[,c(28, 30, 32, 36:41)], as.numeric)
  dmr[,c(2:4, 6:10)] <- sapply(dmr[,c(2:4, 6:10)], as.numeric)
  de[,2:9] <- sapply(de[,2:9], as.numeric)
  
  
  ### check the sample names are the same
  print(identical(rownames(sampleInfo), colnames(methylLev)))
  print(identical(rownames(sampleInfo), colnames(gexpLev)))
  
  #'******************************************************************************
  #' A function to transform RNA-Seq data with VST in DESeq2 package
  #' readCount: RNA-Seq rawcounts in a matrix or in a data frame form
  #'            Rows are genes and columns are samples
  #' filter_thresh: The function filters out genes that have at least one sample
  #'                with counts larger than the 'filter_thresh' value
  #'                e.g., if the 'filter_thresh' = 1, then it removes genes
  #'                that have counts <= 1 across all the samples
  #'                if 0, then there will be no filtering
  #'******************************************************************************
  normalizeRNASEQwithVST <- function(readCount, filter_thresh=1) {
    
    ### load library
    if(!require(DESeq2, quietly = TRUE)) {
      if(!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager")
      BiocManager::install("DESeq2")
      require(DESeq2, quietly = TRUE)
    }
    
    ### make a design matrix for DESeq2 data
    condition <- data.frame(factor(rep("OneClass", ncol(readCount))))
    
    ### Data preparation for DESeq2 format
    deSeqData <- DESeqDataSetFromMatrix(countData=readCount, colData=condition, design= ~0)
    
    if(filter_thresh > 0) {
      ### Remove rubbish rows - this will decrease the number of rows
      keep = apply(counts(deSeqData), 1, function(r){
        return(sum(r > filter_thresh) > 0)
      })
      deSeqData <- deSeqData[keep,]
    }
    
    ### VST
    vsd <- varianceStabilizingTransformation(deSeqData)
    transCnt <- data.frame(assay(vsd), check.names = FALSE)
    
    return (transCnt)
    
  }
  
  ### A function to perform 2D PCA and save a plot
  ### normalizedMat: rows are genes and columns are samples
  ### grp: group information of the samples - color
  ### grp2: group2 information of the samples - shape
  ### num: the number of top genes to be used based on variance (-1 [default]: use all the genes)
  ### component: to draw a plot with PC1 & PC2 or PC2 & PC3
  ### title: title of the plot
  ### suppliment: if it is TRUE, this function additionally generates figures and tables
  ###             related to contributions. For example, which genes are highly contributed to
  ###             the PC1, PC2, etc.
  ### save: if it is saved as png or not
  ### outDir: output directory for the plot
  pca_plot <- function(normalizedMat, grp, grp2=NULL,
                       num = -1, component=c("PC1&PC2", "PC2&PC3"),
                       title="PCA_Plot",
                       suppliment=FALSE,
                       save=TRUE,
                       outDir="./") {
    ### load library
    if(!require(ggfortify, quietly = TRUE)) {
      install.packages("ggfortify")
      library(ggfortify, quietly = TRUE)
    }
    if(!require(FactoMineR, quietly = TRUE)) {
      install.packages("FactoMineR")
      library(FactoMineR, quietly = TRUE)
    }
    if(!require(factoextra, quietly = TRUE)) {
      install.packages("factoextra")
      library(factoextra, quietly = TRUE)
    }
    if(!require(xlsx, quietly = TRUE)) {
      install.packages("xlsx")
      require(xlsx, quietly = TRUE)
    }
    
    ### select the top genes based on variance
    if(num >= 0 && num <= nrow(normalizedMat)) {
      v <- apply(normalizedMat, 1, var)
      v <- v[order(-v)]
      top_genes <- names(v)[1:num]
    } else {
      top_genes <- rownames(normalizedMat)
    }
    
    ### PCA
    pca_result <- PCA(t(normalizedMat[top_genes,]), graph = FALSE)
    colnames(pca_result$ind$coord) <- paste0("PC", 1:ncol(pca_result$ind$coord))
    colnames(pca_result$var$contrib) <- paste0("PC", 1:ncol(pca_result$var$contrib))
    if(is.null(grp2)) {
      pca_group <- data.frame(pca_result$ind$coord, group=grp)
    } else {
      pca_group <- data.frame(pca_result$ind$coord, group=grp, group2=grp2)
    }
    
    colors = rainbow(length(unique(grp)))
    names(colors) = unique(grp)
    
    ### save as png
    if(component[1] == "PC1&PC2") {
      if(is.null(grp2)) {
        p <- ggplot(pca_group,aes(x=PC1,y=PC2,col=group)) +
          labs(title=paste0(title, "_PC1-2")) +
          geom_point(size = 5) +
          # geom_text(aes(label=colnames(normalizedMat)),hjust="inward", vjust="inward") +
          scale_color_manual(values = colors) +
          theme_classic(base_size = 24)
      } else {
        p <- ggplot(pca_group,aes(x=PC1,y=PC2,col=group)) +
          labs(title=paste0(title, "_PC1-2")) +
          geom_point(aes(shape=group2), size = 5) +
          # geom_text(aes(label=colnames(normalizedMat)),hjust="inward", vjust="inward") +
          scale_color_manual(values = colors) +
          theme_classic(base_size = 24)  
      }
      if(save) {
        ggsave(filename = paste0(outDir, title, "_PC1-2", ".png"), width = 10, height = 8, dpi = 300)
      }
      
      if(suppliment) {
        fviz_contrib(pca_result, choice = "var", axes = 1, top = 30)
        ggsave(filename = paste0(outDir, title, "_PC1_contribution.png"), width = 12, height = 8)
        
        fviz_contrib(pca_result, choice = "var", axes = 2, top = 30)
        ggsave(filename = paste0(outDir, title, "_PC2_contribution.png"), width = 12, height = 8)
      }
    } else if(component[1] == "PC2&PC3") {
      if(is.null(grp2)) {
        p <- ggplot(pca_group,aes(x=PC2,y=PC3,col=group)) +
          labs(title=paste0(title, "_PC2-3")) +
          geom_point(size = 5) +
          # geom_text(aes(label=colnames(normalizedMat)),hjust="inward", vjust="inward") +
          scale_color_manual(values = colors) +
          theme_classic(base_size = 24)
      } else {
        p <- ggplot(pca_group,aes(x=PC2,y=PC3,col=group)) +
          labs(title=paste0(title, "_PC2-3")) +
          geom_point(aes(shape=group2), size = 5) +
          # geom_text(aes(label=colnames(normalizedMat)),hjust="inward", vjust="inward") +
          scale_color_manual(values = colors) +
          theme_classic(base_size = 24)  
      }
      if(save) {
        ggsave(filename = paste0(outDir, title, "_PC2-3", ".png"), width = 10, height = 8, dpi = 300)
      }
      
      if(suppliment) {
        fviz_contrib(pca_result, choice = "var", axes = 2, top = 30)
        ggsave(filename = paste0(outDir, title, "_PC2_contribution.png"), width = 12, height = 8)
        
        fviz_contrib(pca_result, choice = "var", axes = 3, top = 30)
        ggsave(filename = paste0(outDir, title, "_PC3_contribution.png"), width = 12, height = 8)
      }
    } else {
      stop("\"component\" parameter should be \"PC1&PC2\" or \"PC2&PC3\"")
    }
    
    if(suppliment) {
      write.xlsx2(data.frame(Gene_Symbol=rownames(pca_result$var$contrib), pca_result$var$contrib,
                             stringsAsFactors = FALSE, check.names = FALSE),
                  file = paste0(outDir, title, "_PC_contribution.xlsx"),
                  sheetName = "PCA_contribution", row.names = FALSE)
    }
    
    if(!save) {
      return(p)
    }
  }
  
  ### create the output directory
  dir.create(outputDir, showWarnings = FALSE, recursive = TRUE)
  
  ### gender, age, race
  #
  # PCA
  #
  # TOP DE or DM - correlation plot & beeswarm plot
  
  ### normalize the gene expression
  normalized_gexp <- normalizeRNASEQwithVST(gexpLev)
  
  ### make empty plot saver
  p <- vector("list", length = 6)
  
  ### pca plot - gender
  sampleInfo$Gender[which(sampleInfo$Gender == "M")] <- "Male"
  sampleInfo$Gender[which(sampleInfo$Gender == "F")] <- "Female"
  p[[1]] <- pca_plot(normalizedMat = data.frame(normalized_gexp,
                                      stringsAsFactors = FALSE, check.names = FALSE),
           grp = as.character(sampleInfo$Gender),
           title = paste0("PCA_GEXP_GENDER"),
           save = FALSE,
           outDir = outputDir)
  p[[2]] <- pca_plot(normalizedMat = data.frame(methylLev,
                                      stringsAsFactors = FALSE, check.names = FALSE),
           grp = as.character(sampleInfo$Gender),
           title = paste0("PCA_Methylation_GENDER"),
           save = FALSE,
           outDir = outputDir)
  
  ### pca plot - race
  p[[3]] <- pca_plot(normalizedMat = data.frame(normalized_gexp,
                                      stringsAsFactors = FALSE, check.names = FALSE),
           grp = as.character(sampleInfo$Race.Ethnicity),
           title = paste0("PCA_GEXP_RACE"),
           save = FALSE,
           outDir = outputDir)
  p[[4]] <- pca_plot(normalizedMat = data.frame(methylLev,
                                      stringsAsFactors = FALSE, check.names = FALSE),
           grp = as.character(sampleInfo$Race.Ethnicity),
           title = paste0("PCA_Methylation_RACE"),
           save = FALSE,
           outDir = outputDir)
  
  ### pca plot - age
  sampleInfo$Age2 <- "Age >= 50"
  sampleInfo$Age2[which(sampleInfo$Age < 50)] <- "Age < 50"
  p[[5]] <- pca_plot(normalizedMat = data.frame(normalized_gexp,
                                      stringsAsFactors = FALSE, check.names = FALSE),
           grp = as.character(sampleInfo$Age2),
           title = paste0("PCA_GEXP_AGE"),
           save = FALSE,
           outDir = outputDir)
  p[[6]] <- pca_plot(normalizedMat = data.frame(methylLev,
                                      stringsAsFactors = FALSE, check.names = FALSE),
           grp = as.character(sampleInfo$Age2),
           title = paste0("PCA_Methylation_AGE"),
           save = FALSE,
           outDir = outputDir)
  
  ### arrange the plots and print out
  g <- grid.arrange(grobs = p, ncol = 2, nrow = 3)
  ggsave(file = paste0(outputDir, "PCA_AGE_GENDER_RACE.png"), g, width = 20, height = 20, dpi = 300)
  
  ### A pie chart and Venn diagram showing number of genes overlapping
  #   (differentially methylated/differentially expressed) among health, gingivitis,
  #   and periodontitis group would be beneficial and easy to understand
  #   in addition to the tables provided.
  
  ### load other EPIC results
  de_ging_heal <- read.xlsx2(file = "./results/DEA/DE_Results_Gingivitis_vs_Healthy.xlsx", sheetIndex = 1, stringsAsFactors = FALSE)
  de_peri_ging <- read.xlsx2(file = "./results/DEA/DE_Results_Periodontitis_vs_Gingivitis.xlsx", sheetIndex = 1, stringsAsFactors = FALSE)
  rownames(de_ging_heal) <- de_ging_heal$Gene_Symbol
  rownames(de_peri_ging) <- de_peri_ging$Gene_Symbol
  
  dmp_ging_heal <- read.table(file = "./results/DMA/DMP_Gingivitis_Healthy.txt", header = TRUE, sep = "\t", check.names = FALSE)
  dmp_peri_ging <- read.table(file = "./results/DMA/DMP_Periodontitis_Gingivitis.txt", header = TRUE, sep = "\t", check.names = FALSE)
  rownames(dmp_ging_heal) <- dmp_ging_heal$Name
  rownames(dmp_peri_ging) <- dmp_peri_ging$Name
  
  ### change character columns to numeric
  de_ging_heal[,2:9] <- sapply(de_ging_heal[,2:9], as.numeric)
  de_peri_ging[,2:9] <- sapply(de_peri_ging[,2:9], as.numeric)
  dmp_ging_heal[,c(28, 30, 32, 36:41)] <- sapply(dmp_ging_heal[,c(28, 30, 32, 36:41)], as.numeric)
  dmp_peri_ging[,c(28, 30, 32, 36:41)] <- sapply(dmp_peri_ging[,c(28, 30, 32, 36:41)], as.numeric)
  
  ### get DE genes
  de_genes1 <- de$Gene_Symbol[which(de$padj < 0.05)]
  de_genes2 <- de_ging_heal$Gene_Symbol[which(de_ging_heal$padj < 0.05)]
  de_genes3 <- de_peri_ging$Gene_Symbol[which(de_peri_ging$padj < 0.05)]
  
  ### get DMPs
  dm_cpgs1 <- dmp$Name[which(dmp$adj.P.Val < 0.05)]
  dm_cpgs2 <- dmp_ging_heal$Name[which(dmp_ging_heal$adj.P.Val < 0.05)]
  dm_cpgs3 <- dmp_peri_ging$Name[which(dmp_peri_ging$adj.P.Val < 0.05)]
  
  ### prepare a color
  colors <- brewer.pal(3, "Pastel2")
  
  ### draw a venn diagram for EPIC DE
  venn.diagram(
    x = list(de_genes1, de_genes2, de_genes3),
    category.names = c("Periodontitis_vs_Health", "Gingivitis_vs_Health", "Periodontitis_vs_Gingivitis"),
    filename = paste0(outputDir, "EPIC_Overlapping_DEGs.png"),
    output=TRUE,
    
    # Output features
    imagetype="png",
    width = 1080,
    height = 960,
    resolution = 300,
    compression = "lzw",
    
    # Circles
    lwd = 2,
    col = "black",
    fill = colors,
    
    # Numbers
    cex = 0.6,
    fontface = "bold",
    fontfamily = "sans",

    # Set names
    cat.cex = 0.6,
    cat.fontface = "bold",
    cat.default.pos = "outer",
    cat.pos = c(-27, 27, 135),
    cat.dist = c(0.055, 0.055, 0.085),
    cat.fontfamily = "sans",
    rotation = 1
  )
  
  ### draw a venn diagram for EPIC DMP
  venn.diagram(
    x = list(dm_cpgs1, dm_cpgs2, dm_cpgs3),
    category.names = c("Periodontitis_vs_Health", "Gingivitis_vs_Health", "Periodontitis_vs_Gingivitis"),
    filename = paste0(outputDir, "EPIC_Overlapping_DMPs.png"),
    output=TRUE,
    
    # Output features
    imagetype="png",
    width = 1080,
    height = 960,
    resolution = 300,
    compression = "lzw",
    
    # Circles
    lwd = 2,
    col = "black",
    fill = colors,
    
    # Numbers
    cex = 0.6,
    fontface = "bold",
    fontfamily = "sans",
    
    # Set names
    cat.cex = 0.6,
    cat.fontface = "bold",
    cat.default.pos = "outer",
    cat.pos = c(-27, 27, 135),
    cat.dist = c(0.055, 0.055, 0.085),
    cat.fontfamily = "sans",
    rotation = 1
  )

  ### 3. Similarly, a bar graph representing functional/regulatory regions
  #      (promoter, enhancer, non-coding etc.) scanned would be beneficial.
  
  ### CpG position relative to regulatory elements
  sig_dmp <- dmp[which(dmp$adj.P.Val < 0.05),]
  
  ### change the column
  sig_dmp$Regulatory_Feature_Group2 <- sig_dmp$Regulatory_Feature_Group
  sig_dmp$Regulatory_Feature_Group2[which(sig_dmp$Regulatory_Feature_Group == "Unclassified_Cell_type_specific")] <- "Unclassified"
  sig_dmp$Regulatory_Feature_Group2[which(sig_dmp$Regulatory_Feature_Group == "Gene_Associated_Cell_type_specific")] <- "Gene_Associated"
  sig_dmp$Regulatory_Feature_Group2[which(sig_dmp$Regulatory_Feature_Group == "Promoter_Associated_Cell_type_specific")] <- "Promoter_Associated"
  sig_dmp$Regulatory_Feature_Group2[which(sig_dmp$Regulatory_Feature_Group == "NonGene_Associated_Cell_type_specific")] <- "NonGene_Associated"
  sig_dmp$Regulatory_Feature_Group2[which(sig_dmp$Regulatory_Feature_Group == "")] <- "Unclassified"
  
  ### unique regulatory features
  unique_features <- unique(sig_dmp$Regulatory_Feature_Group2)
  
  ### make an empty vector
  feature_cnt <- rep(0, length(unique_features))
  names(feature_cnt) <- unique_features
  
  ### count the features
  for(feature in unique_features) {
    feature_cnt[feature] <- length(which(sig_dmp$Regulatory_Feature_Group2 == feature))
  }
  
  ### data frame for the plot
  plot_df <- data.frame(Categories=unique_features,
                        Counts=feature_cnt,
                        stringsAsFactors = FALSE, check.names = FALSE)
  
  ### bar graph
  ggplot(data=plot_df, aes(x=Categories, y=Counts)) +
    geom_bar(aes(fill=Categories), stat="identity") +
    geom_text(aes(label=Counts), vjust=-0.5, color="black", size=10) +
    ggtitle("Regulatory Features of the DMPs") +
    ylim(c(0,41000)) +
    # facet_zoom(ylim = c(0, 1000)) +
    theme_classic(base_size = 50) +
    theme(axis.title.x = element_blank(),
          axis.text = element_text(size = 50),
          axis.text.x = element_text(angle = -45, hjust = 0))
  ggsave(file = paste0(outputDir, "Regulatory_Features_DMPs.png"), width = 20, height = 15)
  
}
