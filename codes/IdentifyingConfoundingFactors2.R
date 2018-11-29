###
#   File name : IdentifyingConfoundingFactors2.R
#   Author    : Hyunjin Kim
#   Date      : Nov 27, 2018
#   Email     : hk2990@cumc.columbia.edu
#   Purpose   : Run PCA and t-SNE to find possible confounding factors on raw counts
#
#   * THIS IS DIFFERENT FROM "IdentifyingConfoundingFactors.R" SINCE THIS IS FOR AFFY CHIP DATA
#     THE "IdentifyingConfoundingFactors.R" IS FOR RNA-SEQ DATA
#
#   Instruction
#               1. Source("IdentifyingConfoundingFactors2.R")
#               2. Run the function "identifyCF2" - specify an input file path (normalized microarray) and output directory
#               3. The PCA and t-SNE plots will be generated under the output directory
#
#   Example
#               > source("The_directory_of_IdentifyingConfoundingFactors2.R/IdentifyingConfoundingFactors2.R")
#               > identifyCF2(normCntPath="./data/Affy/panos_affy_data.rda",
#                             outputDir="./data/Affy/")
###

identifyCF2 <- function(normCntPath="./data/Affy/panos_affy_data.rda",
                        outputDir="./data/Affy/") {
  
  ### load dataset
  load(normCntPath)
  
  
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
    tryCatch({
      writeLines("Perplexity = 30")
      t<- Rtsne(t(normalizedMat), perplexity = 30)
    }, error = function(err) {
      tryCatch({
        writeLines("Perplexity = 10")
        t <- Rtsne(t(normalizedMat), perplexity = 10)
      }, error = function(err) {
        tryCatch({
          writeLines("Perplexity = 5")
          t <- Rtsne(t(normalizedMat), perplexity = 5)
        }, error = function(err) {
          tryCatch({
            writeLines("Perplexity = 3")
            t <- Rtsne(t(normalizedMat), perplexity = 3)
          }, error = function(err) {
            writeLines("Perplexity = 2")
            t <- Rtsne(t(normalizedMat), perplexity = 2)
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
  
  
  ### put "NA" for the empty cells
  affy_norm_ge_sample_info$Diagnosis[which(affy_norm_ge_sample_info$Diagnosis == "")] <- "NA"
  
  
  ### make a pca plot
  
  ### phenotype
  pca_plot(affy_norm_ge[,-1], grp = affy_norm_ge_sample_info$Condition,
           title = "Affy_PCA", filePath = paste0(outputDir, "pca_affy_pheno.png"))
  
  ### diagnosis
  pca_plot(affy_norm_ge[,-1], grp = affy_norm_ge_sample_info$Diagnosis,
           title = "Affy_PCA", filePath = paste0(outputDir, "pca_affy_diagnosis.png"))
  
  
  ### make a t-SNE plot
  
  ### phenotype
  tsne_plot(affy_norm_ge[,-1], grp = affy_norm_ge_sample_info$Condition,
            title = "Affy_TSNE", filePath = paste0(outputDir, "tsne_affy_pheno.png"))
  
  ### diagnosis
  tsne_plot(affy_norm_ge[,-1], grp = affy_norm_ge_sample_info$Diagnosis,
            title = "Affy_TSNE", filePath = paste0(outputDir, "tsne_affy_diagnosis.png"))
  
}


