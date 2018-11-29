###
#   File name : IdentifyingConfoundingFactors.R
#   Author    : Hyunjin Kim
#   Date      : Oct 5, 2018
#   Email     : hk2990@cumc.columbia.edu
#   Purpose   : Run PCA and t-SNE to find possible confounding factors on raw counts
#
#   Instruction
#               1. Source("IdentifyingConfoundingFactors.R")
#               2. Run the function "identifyCF" - specify an input file path (raw counts) and output directory
#               3. The PCA and t-SNE plots will be generated under the output directory
#
#   Example
#               > source("The_directory_of_IdentifyingConfoundingFactors.R/IdentifyingConfoundingFactors.R")
#               > identifyCF(rawCntPath="./data/RNA_Seq/panos_rna_seq_data.rda",
#                            outputDir="./data/RNA_Seq/")
###

identifyCF <- function(rawCntPath="//isilon.c2b2.columbia.edu/ifs/archive/shares/bisr/Papapanou/Sept_2018/PreprocessedData/RNA_Seq/panos_rna_seq_data.rda",
                       outputDir="//isilon.c2b2.columbia.edu/ifs/archive/shares/bisr/Papapanou/Sept_2018/PreprocessedData/RNA_Seq/") {
  
  ### load library
  if(!require(sva)) {
    source("https://bioconductor.org/biocLite.R")
    biocLite("sva")
    library(sva)
  }
  
  
  ### load dataset
  load(rawCntPath)
  
  
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
  
  
  ### A function to perform 2D PCA and save a plot
  pca_plot <- function(normalizedMat, grp, title, component=c("PC1&PC2", "PC2&PC3"), filePath) {
    
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
    
    ### save as png
    png(filename=filePath, width = 1000, height = 800)
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
      dev.off()
      stop("\"component\" parameter should be \"PC1&PC2\" or \"PC2&PC3\"")
    }
    dev.off()
    
  }
  
  
  ### A function to perform t-SNE and save a plot
  tsne_plot <- function(normalizedMat, grp, title, filePath) {
    
    ### load library
    if(!require(Rtsne)) {
      install.packages("Rtsne")
      library(Rtsne)
    }
    
    ### TSNE
    set.seed(1234)
    t <- tryCatch({
      writeLines("Perplexity = 30")
      return(Rtsne(t(normalizedMat), perplexity = 30))
    }, error = function(err) {
      t <- tryCatch({
        writeLines("Perplexity = 10")
        return(Rtsne(t(normalizedMat), perplexity = 10))
      }, error = function(err) {
        t <- tryCatch({
          writeLines("Perplexity = 5")
          return(Rtsne(t(normalizedMat), perplexity = 5))
        }, error = function(err) {
          t <- tryCatch({
            writeLines("Perplexity = 3")
            return(Rtsne(t(normalizedMat), perplexity = 3))
          }, error = function(err) {
            writeLines("Perplexity = 2")
            return(Rtsne(t(normalizedMat), perplexity = 2))
          })
        })
      })
    })
    
    colors = rainbow(length(unique(grp)))
    names(colors) = unique(grp)
    
    ### save a plot
    png(filename=filePath, width = 1000, height = 800)
    plot(t$Y, col="white", xlab="tsne1", ylab="tsne2", main = title)
    text(t$Y, labels=colnames(normalizedMat), col=colors[grp])
    legend("topright", legend = unique(grp), col = colors[unique(grp)], pch = 15)
    dev.off()
    
  }
  
  
  ### make a pca plot
  
  ### phenotype
  pca_plot(normalizeRNASEQwithVST(rawCnt), grp = sampleInfo$Phenotype,
           title = "RNA_Seq_PCA", filePath = paste0(outputDir, "pca_rna_seq_pheno.png"))
  
  ### batch
  pca_plot(normalizeRNASEQwithVST(rawCnt), grp = sampleInfo$Batch,
           title = "RNA_Seq_PCA", filePath = paste0(outputDir, "pca_rna_seq_batch.png"))
  
  ### make a t-SNE plot
  ### phenotype
  tsne_plot(normalizeRNASEQwithVST(rawCnt), grp = sampleInfo$Phenotype,
            title = "RNA_Seq_TSNE", filePath = paste0(outputDir, "tsne_rna_seq_pheno.png"))
  
  ### batch
  tsne_plot(normalizeRNASEQwithVST(rawCnt), grp = sampleInfo$Batch,
            title = "RNA_Seq_TSNE", filePath = paste0(outputDir, "tsne_rna_seq_batch.png"))
  
}


