###
#   File name : PreprocessMethyl.R
#   Author    : Hyunjin Kim
#   Date      : Sep 27, 2018
#   Email     : hk2990@cumc.columbia.edu
#   Purpose   : Convert and normalize the IDAT files of methylation data to readable format
#
#   Workflow
#               https://www.bioconductor.org/packages/devel/workflows/vignettes/methylationArrayAnalysis/inst/doc/methylationArrayAnalysis.html#obtaining-the-data
#
#   Instruction
#               1. Source("PreprocessMethyl.R")
#               2. Run the function "preprocessMethyl" - specify the input directory (.idat), the sample info csv file and output directory
#               3. The converted and normalized methylation data will be generated under the output directory
#
#   Example
#               > source("The_directory_of_PreprocessMethyl.R/PreprocessMethyl.R")
#               > preprocessMethyl(inputDir="E:/Panos/Methylation/rawdata/",
#                                  sampleInfoPath="E:/Panos/Methylation/rawdata/SampleSheet_Hyunjin.csv",
#                                  xReactiveProbePath="E:/Methylation/cross_reactive_probes/xReactiveProbes_EPIC.txt",
#                                  outputDir="./data/Methylation/")
###

preprocessMethyl <- function(inputDir="//isilon.c2b2.columbia.edu/ifs/archive/shares/bisr/Papapanou/Sept_2018/OriginalData/Methylation/rawdata/",
                             sampleInfoPath="//isilon.c2b2.columbia.edu/ifs/archive/shares/bisr/Papapanou/Sept_2018/OriginalData/Methylation/SampleSheet_Hyunjin.csv",
                             xReactiveProbePath="//isilon.c2b2.columbia.edu/ifs/archive/shares/af_lab/Hyunjin/cross_reactive_probes/xReactiveProbes_EPIC.txt",
                             outputDir="//isilon.c2b2.columbia.edu/ifs/archive/shares/bisr/Papapanou/Sept_2018/PreprocessedData/Methylation/") {
  
  ### load libraries
  if(!require(minfi)) {
    source("https://bioconductor.org/biocLite.R")
    biocLite("minfi")
    library(minfi)
  }
  if(!require(IlluminaHumanMethylationEPICmanifest)) {
    source("https://bioconductor.org/biocLite.R")
    biocLite("IlluminaHumanMethylationEPICmanifest")
    library(IlluminaHumanMethylationEPICmanifest)
  }
  if(!require(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)) {
    source("https://bioconductor.org/biocLite.R")
    biocLite("IlluminaHumanMethylationEPICanno.ilm10b2.hg19")
    library(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)
  }
  if(!require(limma)) {
    source("https://bioconductor.org/biocLite.R")
    biocLite("limma")
    library(limma)
  }
  if(!require(sva)) {
    source("https://bioconductor.org/biocLite.R")
    biocLite("sva")
    library(sva)
  }
  if(!require(RColorBrewer)) {
    source("https://bioconductor.org/biocLite.R")
    biocLite("RColorBrewer")
    library(RColorBrewer)
  }
  
  
  ### get target directories
  dirs <- list.files(inputDir)
  
  
  ### get paths of the existing IDAT files
  idatFilePaths <- NULL
  for(i in 1:length(dirs)) {
    f <- list.files(paste0(inputDir, dirs[i]))
    f <- f[which(endsWith(f, ".idat"))]
    idatFilePaths <- c(idatFilePaths, paste0(inputDir, dirs[1], "/", f))
  }
  
  
  ### set sample info
  targets <- read.csv(sampleInfoPath, stringsAsFactors = FALSE, skip = 8)
  targets$Basename <- paste0(inputDir, targets$SentrixBarcode_A, "/", targets$SentrixBarcode_A, "_", targets$SentrixPosition_A)
  
  
  ### put pariodental phenotype to targets
  Phenotype <- c(rep("Periodontitis", 6), rep("Healthy", 6), rep("Gingivitis", 6),
                 rep("Periodontitis", 6), rep("Healthy", 6), rep("Gingivitis", 6))
  targets$Phenotype <- Phenotype
  
  
  ### remove all-NA columns from the sample info
  rIdx <- which(apply(apply(sapply(targets, is.na), 2, as.numeric), 2, sum) == nrow(targets))
  targets <- targets[,-rIdx]
  
  
  ### read the IDAT files
  rgSet <- read.metharray.exp(targets = targets)
  sampleNames(rgSet) <- targets$Sample_Name
  
  
  # ### QC: QC report by minfi package
  # qcReport(rgSet, sampNames=targets$Sample_Name, sampGroups=targets$SentrixBarcode_A, 
  #          pdf=paste0(outputDir, "qc_report.pdf"))
  
  
  ### QC1: calculate detection p-values
  detP <- detectionP(rgSet)
  
  
  ### draw a bar plot with the detection p-values
  pal <- brewer.pal(length(unique(targets$SentrixBarcode_A)),"Pastel1")
  png(paste0(outputDir, "qc_detection_pv.png"), width = 1000, height = 1000)
  barplot(colMeans(detP), col=pal[factor(targets$SentrixBarcode_A)], las=2,
          cex.names=1, main="Mean Detection P-Values", cex.main=3,
          ylim = c(0, max(colMeans(detP))*1.5))
  abline(h=0.05,col="red")
  legend("topright", legend=levels(factor(targets$SentrixBarcode_A)), fill=pal,
         bg="white", title = "Sample Groups", cex = 1.5)
  dev.off()
  
  
  ### no normalization
  rawSet <- preprocessRaw(rgSet)
  
  
  ### quantile normalization
  ### since I do not expect global methylation changes among the phenotypes, quantile normalization has been applied
  ### if they were cancer vs normal, using Funnorm normalization would be recommended
  grSet <- preprocessQuantile(rgSet)
  
  
  ### QC2: effect of normalization
  pal <- brewer.pal(length(unique(targets$SentrixBarcode_A)),"Dark2")
  png(paste0(outputDir, "qc_normalization.png"), width = 1800, height = 1000)
  par(mfrow=c(1,2))
  densityPlot(rawSet, sampGroups=targets$SentrixBarcode_A, main="Raw", cex.main=2.5, legend=FALSE)
  legend("topleft", legend=levels(factor(targets$SentrixBarcode_A)), fill=pal,
         bg="white", title = "Sample Groups", cex = 1.5)
  densityPlot(getBeta(grSet), sampGroups=targets$SentrixBarcode_A, main="Normalized", cex.main=2.5, legend=FALSE)
  legend("topleft", legend=levels(factor(targets$SentrixBarcode_A)), fill=pal,
         bg="white", title = "Sample Groups", cex = 1.5)
  dev.off()
  
  
  ### QC3: clustering - x-axis: log2(Meth median intensity), y-axis: log2(Unmeth median intensity)
  png(paste0(outputDir, "qc_log_median_intensity.png"), width = 1000, height = 1000, res = 120)
  plotQC(getQC(rawSet))
  title(main = "log2(Median Intensities)")
  dev.off()
  
  
  ### get Illumina EPIC hg19 annotation data
  anno <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)
  anno <- anno[match(featureNames(grSet), anno$Name),
               c(1:4,12:19,24:ncol(anno))]
  
  
  ### remove probes on the sex chromosomes
  keep <- !(featureNames(grSet) %in% anno$Name[anno$chr %in% c("chrX", "chrY")])
  grSet <- grSet[keep,]
  
  
  ### remove probes with SNPs at CpG site
  grSet <- dropLociWithSnps(grSet)
  
  
  ### remove cross-reactive probes
  xReactiveProbes <- as.character(read.table(xReactiveProbePath, header = FALSE, sep = "\t")[,1])
  keep <- !(featureNames(grSet) %in% xReactiveProbes)
  grSet <- grSet[keep,]
  
  
  ### get beta and m values
  beta <- getBeta(grSet)
  m <- getM(grSet)
  
  
  ### batch effect correction - sample group
  beta_batch_corrected <- ComBat(dat = beta,
                                 batch = targets$SentrixBarcode_A,
                                 mod = model.matrix(~1, data = targets))
  m_batch_corrected <- ComBat(dat = m,
                              batch = targets$SentrixBarcode_A,
                              mod = model.matrix(~1, data = targets))
  
  
  ### add RNA-Seq batch info - just in case to examine
  targets$Batch <- c(rep("A", 18), rep("B", 18))
  
  
  ### QC4: PCA (MDS) plot of beta and m values
  
  ### sample group - before batch effect correction
  png(paste0(outputDir, "qc_pca_group.png"), width = 1800, height = 1000)
  par(mfrow=c(1,2))
  ### beta
  plotMDS(beta, top = 1000, gene.selection = "common", main = "Beta",
          col=pal[factor(targets$SentrixBarcode_A)], cex = 1.5)
  legend("topleft", legend=levels(factor(targets$SentrixBarcode_A)), fill=pal,
         bg="white", title = "Sample Groups", cex = 1)
  ### m
  plotMDS(m, top = 1000, gene.selection = "common", main = "M",
          col=pal[factor(targets$SentrixBarcode_A)], cex = 1.5)
  legend("topleft", legend=levels(factor(targets$SentrixBarcode_A)), fill=pal,
         bg="white", title = "Sample Groups", cex = 1)
  dev.off()
  
  ### sample group - after batch effect correction
  png(paste0(outputDir, "qc_pca_group_bc.png"), width = 1800, height = 1000)
  par(mfrow=c(1,2))
  ### beta
  plotMDS(beta_batch_corrected, top = 1000, gene.selection = "common", main = "Batch Effect Corrected Beta",
          col=pal[factor(targets$SentrixBarcode_A)], cex = 1.5)
  legend("topleft", legend=levels(factor(targets$SentrixBarcode_A)), fill=pal,
         bg="white", title = "Sample Groups", cex = 1)
  ### m
  plotMDS(m_batch_corrected, top = 1000, gene.selection = "common", main = "Batch Effect Corrected M",
          col=pal[factor(targets$SentrixBarcode_A)], cex = 1.5)
  legend("topleft", legend=levels(factor(targets$SentrixBarcode_A)), fill=pal,
         bg="white", title = "Sample Groups", cex = 1)
  dev.off()
  
  ### sample phenotype - before batch effect correction
  pal <- brewer.pal(length(unique(targets$Phenotype)),"Pastel1")
  png(paste0(outputDir, "qc_pca_phenotype.png"), width = 1800, height = 1000)
  par(mfrow=c(1,2))
  ### beta
  plotMDS(beta, top = 1000, gene.selection = "common", main = "Beta",
          col=pal[factor(targets$Phenotype)], cex = 1.5)
  legend("topleft", legend=levels(factor(targets$Phenotype)), fill=pal,
         bg="white", title = "Sample Phenotypes", cex = 1)
  ### m
  plotMDS(m, top = 1000, gene.selection = "common", main = "M",
          col=pal[factor(targets$Phenotype)], cex = 1.5)
  legend("topleft", legend=levels(factor(targets$Phenotype)), fill=pal,
         bg="white", title = "Sample Phenotypes", cex = 1)
  dev.off()
  
  ### sample phenotype - after batch effect correction
  pal <- brewer.pal(length(unique(targets$Phenotype)),"Pastel1")
  png(paste0(outputDir, "qc_pca_phenotype_bc.png"), width = 1800, height = 1000)
  par(mfrow=c(1,2))
  ### beta
  plotMDS(beta_batch_corrected, top = 1000, gene.selection = "common", main = "Batch Effect Corrected Beta",
          col=pal[factor(targets$Phenotype)], cex = 1.5)
  legend("topleft", legend=levels(factor(targets$Phenotype)), fill=pal,
         bg="white", title = "Sample Phenotypes", cex = 1)
  ### m
  plotMDS(m_batch_corrected, top = 1000, gene.selection = "common", main = "Batch Effect Corrected M",
          col=pal[factor(targets$Phenotype)], cex = 1.5)
  legend("topleft", legend=levels(factor(targets$Phenotype)), fill=pal,
         bg="white", title = "Sample Phenotypes", cex = 1)
  dev.off()
  
  ### row number of Sentrix Position - before batch effect correction
  pal <- brewer.pal(length(unique(targets$SentrixPosition_A)),"Pastel1")
  png(paste0(outputDir, "qc_pca_rowNum.png"), width = 1800, height = 1000)
  par(mfrow=c(1,2))
  ### beta
  plotMDS(beta, top = 1000, gene.selection = "common", main = "Beta",
          col=pal[factor(targets$SentrixPosition_A)], cex = 1.5)
  legend("topleft", legend=levels(factor(substring(targets$SentrixPosition_A, 1, 3))), fill=pal,
         bg="white", title = "Row Number", cex = 1)
  ### m
  plotMDS(m, top = 1000, gene.selection = "common", main = "M",
          col=pal[factor(targets$SentrixPosition_A)], cex = 1.5)
  legend("topleft", legend=levels(factor(substring(targets$SentrixPosition_A, 1, 3))), fill=pal,
         bg="white", title = "Row Number", cex = 1)
  dev.off()
  
  ### RNA-Seq batch info - before batch effect correction
  pal <- c("red", "blue")
  names(pal) <- c("A", "B")
  png(paste0(outputDir, "qc_pca_batch.png"), width = 1800, height = 1000)
  par(mfrow=c(1,2))
  ### beta
  plotMDS(beta, top = 1000, gene.selection = "common", main = "Beta",
          col=pal[factor(targets$Batch)], cex = 1.5)
  legend("topleft", legend=levels(factor(targets$Batch)), fill=pal,
         bg="white", title = "RNA-Seq Batch", cex = 1)
  ### m
  plotMDS(m, top = 1000, gene.selection = "common", main = "M",
          col=pal[factor(targets$Batch)], cex = 1.5)
  legend("topleft", legend=levels(factor(targets$Batch)), fill=pal,
         bg="white", title = "RNA-Seq Batch", cex = 1)
  dev.off()
  
  
  ### QC5: TSNE plot of beta and m values
  
  ### a function to perform t-SNE and save a plot
  tsne_plot <- function(normalizedMat, grp, featureNum, title, filePath=NULL) {
    
    ### load library
    if(!require(Rtsne)) {
      install.packages("Rtsne")
      library(Rtsne)
    }
    
    ### A function to select top variance features
    selectTopV <- function(x, selectNum) {
      v <- apply(x, 1, var)
      x <- x[order(-v),]
      if(selectNum != Inf) {
        x <- x[1:selectNum,]
      }
      
      return (x)
    }
    
    ### choose top variance features
    filteredMat <- selectTopV(normalizedMat, featureNum)
    
    ### TSNE
    set.seed(1234)
    t <- tryCatch({
      writeLines("Perplexity = 30")
      return(Rtsne(t(filteredMat), perplexity = 30))
    }, error = function(err) {
      t <- tryCatch({
        writeLines("Perplexity = 10")
        return(Rtsne(t(filteredMat), perplexity = 10))
      }, error = function(err) {
        t <- tryCatch({
          writeLines("Perplexity = 5")
          return(Rtsne(t(filteredMat), perplexity = 5))
        }, error = function(err) {
          t <- tryCatch({
            writeLines("Perplexity = 3")
            return(Rtsne(t(filteredMat), perplexity = 3))
          }, error = function(err) {
            writeLines("Perplexity = 2")
            return(Rtsne(t(filteredMat), perplexity = 2))
          })
        })
      })
    })
    
    colors = rainbow(length(unique(grp)))
    names(colors) = unique(grp)
    
    if(is.null(filePath)) {
      ### print a plot
      plot(t$Y, col="white", xlab="tsne1", ylab="tsne2", main = title)
      text(t$Y, labels=colnames(filteredMat), col=colors[as.character(grp)])
      legend("topleft", legend = unique(grp), col = colors[as.character(unique(grp))], pch = 15)
    } else {
      ### save a plot
      png(filename=filePath, width = 1000, height = 800)
      plot(t$Y, col="white", xlab="tsne1", ylab="tsne2", main = title)
      text(t$Y, labels=colnames(filteredMat), col=colors[as.character(grp)])
      legend("topleft", legend = unique(grp), col = colors[as.character(unique(grp))], pch = 15)
      dev.off()
    }
    
  }
  
  ### sample group - before batch effect correction
  png(paste0(outputDir, "qc_tsne_group.png"), width = 1800, height = 1000)
  par(mfrow=c(1,2))
  ### beta
  tsne_plot(beta, grp = targets$SentrixBarcode_A, featureNum = 1000, title = "Beta")
  ### m
  tsne_plot(m, grp = targets$SentrixBarcode_A, featureNum = 1000, title = "M")
  dev.off()
  
  ### sample group - after batch effect correction
  png(paste0(outputDir, "qc_tsne_group_bc.png"), width = 1800, height = 1000)
  par(mfrow=c(1,2))
  ### beta
  tsne_plot(beta_batch_corrected, grp = targets$SentrixBarcode_A, featureNum = 1000, title = "Batch Effect Corrected Beta")
  ### m
  tsne_plot(m_batch_corrected, grp = targets$SentrixBarcode_A, featureNum = 1000, title = "Batch Effect Corrected M")
  dev.off()
  
  ### sample phenotype - before batch effect correction
  png(paste0(outputDir, "qc_tsne_phenotype.png"), width = 1800, height = 1000)
  par(mfrow=c(1,2))
  ### beta
  tsne_plot(beta, grp = targets$Phenotype, featureNum = 1000, title = "Beta")
  ### m
  tsne_plot(m, grp = targets$Phenotype, featureNum = 1000, title = "M")
  dev.off()
  
  ### sample phenotype - after batch effect correction
  png(paste0(outputDir, "qc_tsne_phenotype_bc.png"), width = 1800, height = 1000)
  par(mfrow=c(1,2))
  ### beta
  tsne_plot(beta_batch_corrected, grp = targets$Phenotype, featureNum = 1000, title = "Batch Effect Corrected Beta")
  ### m
  tsne_plot(m_batch_corrected, grp = targets$Phenotype, featureNum = 1000, title = "Batch Effect Corrected M")
  dev.off()
  
  ### row number of Sentrix Position - before batch effect correction
  png(paste0(outputDir, "qc_tsne_rowNum.png"), width = 1800, height = 1000)
  par(mfrow=c(1,2))
  ### beta
  tsne_plot(beta, grp = substring(targets$SentrixPosition_A, 1, 3), featureNum = 1000, title = "Beta")
  ### m
  tsne_plot(m, grp = substring(targets$SentrixPosition_A, 1, 3), featureNum = 1000, title = "M")
  dev.off()
  
  ### RNA-Seq batch info - before batch effect correction
  png(paste0(outputDir, "qc_tsne_batch.png"), width = 1800, height = 1000)
  par(mfrow=c(1,2))
  ### beta
  tsne_plot(beta, grp = targets$Batch, featureNum = 1000, title = "Beta")
  ### m
  tsne_plot(m, grp = targets$Batch, featureNum = 1000, title = "M")
  dev.off()
  
  
  ### PCA & TSNE specifically for each phenotype
  
  ### Healthy
  ### sample group - before batch effect correction
  ### PCA
  png(paste0(outputDir, "qc_pca_group_healthy.png"), width = 1800, height = 1000)
  par(mfrow=c(1,2))
  ### beta
  plotMDS(beta[,which(targets$Phenotype == "Healthy")], top = 1000, gene.selection = "common", main = "Beta",
          col=pal[factor(targets$SentrixBarcode_A[which(targets$Phenotype == "Healthy")])], cex = 1.5)
  legend("topleft", legend=levels(factor(targets$SentrixBarcode_A[which(targets$Phenotype == "Healthy")])), fill=pal,
         bg="white", title = "Sample Groups", cex = 1)
  ### m
  plotMDS(m[,which(targets$Phenotype == "Healthy")], top = 1000, gene.selection = "common", main = "M",
          col=pal[factor(targets$SentrixBarcode_A[which(targets$Phenotype == "Healthy")])], cex = 1.5)
  legend("topleft", legend=levels(factor(targets$SentrixBarcode_A[which(targets$Phenotype == "Healthy")])), fill=pal,
         bg="white", title = "Sample Groups", cex = 1)
  dev.off()
  
  ### TSNE
  png(paste0(outputDir, "qc_tsne_group_healthy.png"), width = 1800, height = 1000)
  par(mfrow=c(1,2))
  ### beta
  tsne_plot(beta[,which(targets$Phenotype == "Healthy")], grp = targets$SentrixBarcode_A[which(targets$Phenotype == "Healthy")], featureNum = 1000, title = "Beta")
  ### m
  tsne_plot(m[,which(targets$Phenotype == "Healthy")], grp = targets$SentrixBarcode_A[which(targets$Phenotype == "Healthy")], featureNum = 1000, title = "M")
  dev.off()
  
  ### sample group - after batch effect correction
  ### PCA
  png(paste0(outputDir, "qc_pca_group_bc_healthy.png"), width = 1800, height = 1000)
  par(mfrow=c(1,2))
  ### beta
  plotMDS(beta_batch_corrected[,which(targets$Phenotype == "Healthy")], top = 1000, gene.selection = "common", main = "Batch Effect Corrected Beta",
          col=pal[factor(targets$SentrixBarcode_A[which(targets$Phenotype == "Healthy")])], cex = 1.5)
  legend("topleft", legend=levels(factor(targets$SentrixBarcode_A[which(targets$Phenotype == "Healthy")])), fill=pal,
         bg="white", title = "Sample Groups", cex = 1)
  ### m
  plotMDS(m_batch_corrected[,which(targets$Phenotype == "Healthy")], top = 1000, gene.selection = "common", main = "Batch Effect Corrected M",
          col=pal[factor(targets$SentrixBarcode_A[which(targets$Phenotype == "Healthy")])], cex = 1.5)
  legend("topleft", legend=levels(factor(targets$SentrixBarcode_A[which(targets$Phenotype == "Healthy")])), fill=pal,
         bg="white", title = "Sample Groups", cex = 1)
  dev.off()
  
  ### TSNE
  ### sample group - after batch effect correction
  png(paste0(outputDir, "qc_tsne_group_bc_healthy.png"), width = 1800, height = 1000)
  par(mfrow=c(1,2))
  ### beta
  tsne_plot(beta_batch_corrected[,which(targets$Phenotype == "Healthy")], grp = targets$SentrixBarcode_A[which(targets$Phenotype == "Healthy")], featureNum = 1000, title = "Batch Effect Corrected Beta")
  ### m
  tsne_plot(m_batch_corrected[,which(targets$Phenotype == "Healthy")], grp = targets$SentrixBarcode_A[which(targets$Phenotype == "Healthy")], featureNum = 1000, title = "Batch Effect Corrected M")
  dev.off()
  
  ### RNA-Seq batch info - before batch effect correction
  ### PCA
  png(paste0(outputDir, "qc_pca_batch_healthy.png"), width = 1800, height = 1000)
  par(mfrow=c(1,2))
  ### beta
  plotMDS(beta[,which(targets$Phenotype == "Healthy")], top = 1000, gene.selection = "common", main = "Beta",
          col=pal[factor(targets$Batch[which(targets$Phenotype == "Healthy")])], cex = 1.5)
  legend("topleft", legend=levels(factor(targets$Batch[which(targets$Phenotype == "Healthy")])), fill=pal,
         bg="white", title = "Sample Groups", cex = 1)
  ### m
  plotMDS(m[,which(targets$Phenotype == "Healthy")], top = 1000, gene.selection = "common", main = "M",
          col=pal[factor(targets$Batch[which(targets$Phenotype == "Healthy")])], cex = 1.5)
  legend("topleft", legend=levels(factor(targets$Batch[which(targets$Phenotype == "Healthy")])), fill=pal,
         bg="white", title = "Sample Groups", cex = 1)
  dev.off()
  
  ### TSNE
  png(paste0(outputDir, "qc_tsne_batch_healthy.png"), width = 1800, height = 1000)
  par(mfrow=c(1,2))
  ### beta
  tsne_plot(beta[,which(targets$Phenotype == "Healthy")], grp = targets$Batch[which(targets$Phenotype == "Healthy")], featureNum = 1000, title = "Beta")
  ### m
  tsne_plot(m[,which(targets$Phenotype == "Healthy")], grp = targets$Batch[which(targets$Phenotype == "Healthy")], featureNum = 1000, title = "M")
  dev.off()
  
  
  ### Gingivitis
  ### sample group - before batch effect correction
  ### PCA
  png(paste0(outputDir, "qc_pca_group_gingivitis.png"), width = 1800, height = 1000)
  par(mfrow=c(1,2))
  ### beta
  plotMDS(beta[,which(targets$Phenotype == "Gingivitis")], top = 1000, gene.selection = "common", main = "Beta",
          col=pal[factor(targets$SentrixBarcode_A[which(targets$Phenotype == "Gingivitis")])], cex = 1.5)
  legend("topleft", legend=levels(factor(targets$SentrixBarcode_A[which(targets$Phenotype == "Gingivitis")])), fill=pal,
         bg="white", title = "Sample Groups", cex = 1)
  ### m
  plotMDS(m[,which(targets$Phenotype == "Gingivitis")], top = 1000, gene.selection = "common", main = "M",
          col=pal[factor(targets$SentrixBarcode_A[which(targets$Phenotype == "Gingivitis")])], cex = 1.5)
  legend("topleft", legend=levels(factor(targets$SentrixBarcode_A[which(targets$Phenotype == "Gingivitis")])), fill=pal,
         bg="white", title = "Sample Groups", cex = 1)
  dev.off()
  
  ### TSNE
  png(paste0(outputDir, "qc_tsne_group_gingivitis.png"), width = 1800, height = 1000)
  par(mfrow=c(1,2))
  ### beta
  tsne_plot(beta[,which(targets$Phenotype == "Gingivitis")], grp = targets$SentrixBarcode_A[which(targets$Phenotype == "Gingivitis")], featureNum = 1000, title = "Beta")
  ### m
  tsne_plot(m[,which(targets$Phenotype == "Gingivitis")], grp = targets$SentrixBarcode_A[which(targets$Phenotype == "Gingivitis")], featureNum = 1000, title = "M")
  dev.off()
  
  ### sample group - after batch effect correction
  ### PCA
  png(paste0(outputDir, "qc_pca_group_bc_gingivitis.png"), width = 1800, height = 1000)
  par(mfrow=c(1,2))
  ### beta
  plotMDS(beta_batch_corrected[,which(targets$Phenotype == "Gingivitis")], top = 1000, gene.selection = "common", main = "Batch Effect Corrected Beta",
          col=pal[factor(targets$SentrixBarcode_A[which(targets$Phenotype == "Gingivitis")])], cex = 1.5)
  legend("topleft", legend=levels(factor(targets$SentrixBarcode_A[which(targets$Phenotype == "Gingivitis")])), fill=pal,
         bg="white", title = "Sample Groups", cex = 1)
  ### m
  plotMDS(m_batch_corrected[,which(targets$Phenotype == "Gingivitis")], top = 1000, gene.selection = "common", main = "Batch Effect Corrected M",
          col=pal[factor(targets$SentrixBarcode_A[which(targets$Phenotype == "Gingivitis")])], cex = 1.5)
  legend("topleft", legend=levels(factor(targets$SentrixBarcode_A[which(targets$Phenotype == "Gingivitis")])), fill=pal,
         bg="white", title = "Sample Groups", cex = 1)
  dev.off()
  
  ### TSNE
  ### sample group - after batch effect correction
  png(paste0(outputDir, "qc_tsne_group_bc_gingivitis.png"), width = 1800, height = 1000)
  par(mfrow=c(1,2))
  ### beta
  tsne_plot(beta_batch_corrected[,which(targets$Phenotype == "Gingivitis")], grp = targets$SentrixBarcode_A[which(targets$Phenotype == "Gingivitis")], featureNum = 1000, title = "Batch Effect Corrected Beta")
  ### m
  tsne_plot(m_batch_corrected[,which(targets$Phenotype == "Gingivitis")], grp = targets$SentrixBarcode_A[which(targets$Phenotype == "Gingivitis")], featureNum = 1000, title = "Batch Effect Corrected M")
  dev.off()
  
  ### RNA-Seq batch info - before batch effect correction
  ### PCA
  png(paste0(outputDir, "qc_pca_batch_gingivitis.png"), width = 1800, height = 1000)
  par(mfrow=c(1,2))
  ### beta
  plotMDS(beta[,which(targets$Phenotype == "Gingivitis")], top = 1000, gene.selection = "common", main = "Beta",
          col=pal[factor(targets$Batch[which(targets$Phenotype == "Gingivitis")])], cex = 1.5)
  legend("topleft", legend=levels(factor(targets$Batch[which(targets$Phenotype == "Gingivitis")])), fill=pal,
         bg="white", title = "Sample Groups", cex = 1)
  ### m
  plotMDS(m[,which(targets$Phenotype == "Gingivitis")], top = 1000, gene.selection = "common", main = "M",
          col=pal[factor(targets$Batch[which(targets$Phenotype == "Gingivitis")])], cex = 1.5)
  legend("topleft", legend=levels(factor(targets$Batch[which(targets$Phenotype == "Gingivitis")])), fill=pal,
         bg="white", title = "Sample Groups", cex = 1)
  dev.off()
  
  ### TSNE
  png(paste0(outputDir, "qc_tsne_batch_gingivitis.png"), width = 1800, height = 1000)
  par(mfrow=c(1,2))
  ### beta
  tsne_plot(beta[,which(targets$Phenotype == "Gingivitis")], grp = targets$Batch[which(targets$Phenotype == "Gingivitis")], featureNum = 1000, title = "Beta")
  ### m
  tsne_plot(m[,which(targets$Phenotype == "Gingivitis")], grp = targets$Batch[which(targets$Phenotype == "Gingivitis")], featureNum = 1000, title = "M")
  dev.off()
  
  
  ### Periodontitis
  ### sample group - before batch effect correction
  ### PCA
  png(paste0(outputDir, "qc_pca_group_periodontitis.png"), width = 1800, height = 1000)
  par(mfrow=c(1,2))
  ### beta
  plotMDS(beta[,which(targets$Phenotype == "Periodontitis")], top = 1000, gene.selection = "common", main = "Beta",
          col=pal[factor(targets$SentrixBarcode_A[which(targets$Phenotype == "Periodontitis")])], cex = 1.5)
  legend("topleft", legend=levels(factor(targets$SentrixBarcode_A[which(targets$Phenotype == "Periodontitis")])), fill=pal,
         bg="white", title = "Sample Groups", cex = 1)
  ### m
  plotMDS(m[,which(targets$Phenotype == "Periodontitis")], top = 1000, gene.selection = "common", main = "M",
          col=pal[factor(targets$SentrixBarcode_A[which(targets$Phenotype == "Periodontitis")])], cex = 1.5)
  legend("topleft", legend=levels(factor(targets$SentrixBarcode_A[which(targets$Phenotype == "Periodontitis")])), fill=pal,
         bg="white", title = "Sample Groups", cex = 1)
  dev.off()
  
  ### TSNE
  png(paste0(outputDir, "qc_tsne_group_periodontitis.png"), width = 1800, height = 1000)
  par(mfrow=c(1,2))
  ### beta
  tsne_plot(beta[,which(targets$Phenotype == "Periodontitis")], grp = targets$SentrixBarcode_A[which(targets$Phenotype == "Periodontitis")], featureNum = 1000, title = "Beta")
  ### m
  tsne_plot(m[,which(targets$Phenotype == "Periodontitis")], grp = targets$SentrixBarcode_A[which(targets$Phenotype == "Periodontitis")], featureNum = 1000, title = "M")
  dev.off()
  
  ### sample group - after batch effect correction
  ### PCA
  png(paste0(outputDir, "qc_pca_group_bc_periodontitis.png"), width = 1800, height = 1000)
  par(mfrow=c(1,2))
  ### beta
  plotMDS(beta_batch_corrected[,which(targets$Phenotype == "Periodontitis")], top = 1000, gene.selection = "common", main = "Batch Effect Corrected Beta",
          col=pal[factor(targets$SentrixBarcode_A[which(targets$Phenotype == "Periodontitis")])], cex = 1.5)
  legend("topleft", legend=levels(factor(targets$SentrixBarcode_A[which(targets$Phenotype == "Periodontitis")])), fill=pal,
         bg="white", title = "Sample Groups", cex = 1)
  ### m
  plotMDS(m_batch_corrected[,which(targets$Phenotype == "Periodontitis")], top = 1000, gene.selection = "common", main = "Batch Effect Corrected M",
          col=pal[factor(targets$SentrixBarcode_A[which(targets$Phenotype == "Periodontitis")])], cex = 1.5)
  legend("topleft", legend=levels(factor(targets$SentrixBarcode_A[which(targets$Phenotype == "Periodontitis")])), fill=pal,
         bg="white", title = "Sample Groups", cex = 1)
  dev.off()
  
  ### TSNE
  ### sample group - after batch effect correction
  png(paste0(outputDir, "qc_tsne_group_bc_periodontitis.png"), width = 1800, height = 1000)
  par(mfrow=c(1,2))
  ### beta
  tsne_plot(beta_batch_corrected[,which(targets$Phenotype == "Periodontitis")], grp = targets$SentrixBarcode_A[which(targets$Phenotype == "Periodontitis")], featureNum = 1000, title = "Batch Effect Corrected Beta")
  ### m
  tsne_plot(m_batch_corrected[,which(targets$Phenotype == "Periodontitis")], grp = targets$SentrixBarcode_A[which(targets$Phenotype == "Periodontitis")], featureNum = 1000, title = "Batch Effect Corrected M")
  dev.off()
  
  ### RNA-Seq batch info - before batch effect correction
  ### PCA
  png(paste0(outputDir, "qc_pca_batch_periodontitis.png"), width = 1800, height = 1000)
  par(mfrow=c(1,2))
  ### beta
  plotMDS(beta[,which(targets$Phenotype == "Periodontitis")], top = 1000, gene.selection = "common", main = "Beta",
          col=pal[factor(targets$Batch[which(targets$Phenotype == "Periodontitis")])], cex = 1.5)
  legend("topleft", legend=levels(factor(targets$Batch[which(targets$Phenotype == "Periodontitis")])), fill=pal,
         bg="white", title = "Sample Groups", cex = 1)
  ### m
  plotMDS(m[,which(targets$Phenotype == "Periodontitis")], top = 1000, gene.selection = "common", main = "M",
          col=pal[factor(targets$Batch[which(targets$Phenotype == "Periodontitis")])], cex = 1.5)
  legend("topleft", legend=levels(factor(targets$Batch[which(targets$Phenotype == "Periodontitis")])), fill=pal,
         bg="white", title = "Sample Groups", cex = 1)
  dev.off()
  
  ### TSNE
  png(paste0(outputDir, "qc_tsne_batch_periodontitis.png"), width = 1800, height = 1000)
  par(mfrow=c(1,2))
  ### beta
  tsne_plot(beta[,which(targets$Phenotype == "Periodontitis")], grp = targets$Batch[which(targets$Phenotype == "Periodontitis")], featureNum = 1000, title = "Beta")
  ### m
  tsne_plot(m[,which(targets$Phenotype == "Periodontitis")], grp = targets$Batch[which(targets$Phenotype == "Periodontitis")], featureNum = 1000, title = "M")
  dev.off()
  
  
  ### PCA & TSNE
  ### One clustering the methylation samples in batch A.
  ### The second clustering the methylation sample in batch B.
  ### In both cases, color the samples in the plot according to phenotype. 
  
  ### PCA
  pal <- brewer.pal(length(unique(targets$Phenotype)),"Dark2")
  png(paste0(outputDir, "qc_pca_phenotype_AB.png"), width = 1800, height = 1000)
  par(mfrow=c(1,2))
  ### batch A
  plotMDS(m[,which(targets$Batch == "A")], top = 1000, gene.selection = "common", main = "M - Batch A",
          col=pal[factor(targets$Phenotype[which(targets$Batch == "A")])], cex = 1.5)
  legend("topleft", legend=levels(factor(targets$Phenotype[which(targets$Batch == "A")])), fill=pal,
         bg="white", title = "Sample Phenotypes", cex = 1)
  ### batch B
  plotMDS(m[,which(targets$Batch == "B")], top = 1000, gene.selection = "common", main = "M - Batch B",
          col=pal[factor(targets$Phenotype[which(targets$Batch == "B")])], cex = 1.5)
  legend("topleft", legend=levels(factor(targets$Phenotype[which(targets$Batch == "B")])), fill=pal,
         bg="white", title = "Sample Phenotypes", cex = 1)
  dev.off()
  
  ### TSNE
  png(paste0(outputDir, "qc_tsne_phenotype_AB.png"), width = 1800, height = 1000)
  par(mfrow=c(1,2))
  ### batch A
  tsne_plot(m[,which(targets$Batch == "A")], grp = targets$Phenotype[which(targets$Batch == "A")], featureNum = 1000, title = "M - Batch A")
  ### batch B
  tsne_plot(m[,which(targets$Batch == "B")], grp = targets$Phenotype[which(targets$Batch == "B")], featureNum = 1000, title = "M - Batch B")
  dev.off()
  
  
  # ### gender
  # png(paste0(outputDir, "qc_pca_gender.png"), width = 1800, height = 1000)
  # par(mfrow=c(1,2))
  # ### beta
  # plotMDS(beta, top = 1000, gene.selection = "common", main = "Beta",
  #         col=pal[factor(targets$Gender)], cex = 1.5)
  # legend("topleft", legend=levels(factor(targets$Gender)), fill=pal,
  #        bg="white", title = "Gender", cex = 1)
  # ### m
  # plotMDS(m, top = 1000, gene.selection = "common", main = "M",
  #         col=pal[factor(targets$Gender)], cex = 1.5)
  # legend("topleft", legend=levels(factor(targets$Gender)), fill=pal,
  #        bg="white", title = "Gender", cex = 1)
  # dev.off()
  
  
  # ### QC: gender prediction
  # plotSex(getSex(grSet, cutoff = -2))
  
  
  ### save the beta & m values and sample information
  ### RDA file
  sample_info <- targets
  ### make a README function for the RDA file
  README = function(){
    writeLines(paste(rep("#", 100), collapse = ""))
    writeLines("The methylation dataset is from Dr. Panos Papapanou")
    writeLines("Array: Illumina EPIC Chip")
    writeLines("Reference Used: hg19")
    writeLines("The number of CpG sites: 866091")
    writeLines("Normalization method: Quantile Normalization")
    writeLines(paste(rep("#", 100), collapse = ""))
    writeLines("The \"beta\" object has beta values of the methylation data")
    writeLines("The \"m\" object has m values of the methylation data")
    writeLines("The \"sample_info\" object has all the sample information")
    writeLines("The \"grSet\" object is a \"GenomicRatioSet\" of \"minfi\" package")
    writeLines(paste(rep("#", 100), collapse = ""))
  }
  ### save in .rda
  save(list = c("grSet", "beta", "m", "sample_info", "README"),
       file = paste0(outputDir, "panos_methylation_data.rda"))
  
  ### text files
  beta <- data.frame(CpG_ID=rownames(beta), beta)
  write.table(beta, file = paste0(outputDir, "beta_values.txt"),
              sep = "\t", row.names = FALSE, col.names = TRUE)
  m <- data.frame(CpG_ID=rownames(m), m)
  write.table(m, file =  paste0(outputDir, "m_values.txt"),
              sep = "\t", row.names = FALSE, col.names = TRUE)
  write.table(sample_info, file = paste0(outputDir, "sample_info.txt"),
              sep = "\t", row.names = FALSE)
  
}

