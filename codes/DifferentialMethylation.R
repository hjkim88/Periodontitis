###
#   File name : DifferentialMethylation.R
#   Author    : Hyunjin Kim
#   Date      : Oct 1, 2018
#   Email     : hk2990@cumc.columbia.edu
#   Purpose   : Perform DM analysis on the normalized methylation data
#
#   Workflow
#               https://www.bioconductor.org/packages/devel/workflows/vignettes/methylationArrayAnalysis/inst/doc/methylationArrayAnalysis.html#obtaining-the-data
#
#               1. There are 3 comparisons (Periodontitis vs Healthy, Gingivitis vs Healthy, and Periodontitis vs Gingivitis),
#                  and "limma" is used to get differentially methylated sites.
#                  There is a batch effect among different sample groups, therefore, the effect is considered
#                  as a confoudning factor and included in the design matrix of the linear model.
#                  For the design matrix, [~0+phenotype+samplegroup] is used.
#                  "IlluminaHumanMethylationEPICanno.ilm10b2.hg19" is used to annotate the CpG sites.
#               2. After getting the differentially methylated sites (~850k), "DMRcate" is used to generate
#                  differentially methylated regions. 1000 is used for the lambda parameter in the dmrcate()
#                  function which means gaps >= lambda between significant CpG sites will be in separate DMRs.
#               3. The DMR.plot() draws a plot of the top differentially methylated reguions associated with
#                  chromosome location and nearby genes. It also represents methylation value of the samples.
#               4. Pathway analysis is performed based on genes asoociated with the differentially
#                  methylated CpG sites. Go and Kegg databases are used.
#
#   Instruction
#               1. Source("DifferentialMethylation.R")
#               2. Run the function "dma" - specify the input file path (RDA file) and output directory
#               3. Differentially methylated results will be generated under the output directory
#
#   Example
#               > source("The_directory_of_DifferentialMethylation.R/DifferentialMethylation.R")
#               > dma(inputFilelPath="./data/Methylation/panos_methylation_data.rda",
#                     fdrThreshold=0.05,
#                     dmrPrintNum=5,
#                     outputDir="./results/DMA/")
###

dma <- function(inputFilelPath="//isilon.c2b2.columbia.edu/ifs/archive/shares/bisr/Papapanou/Sept_2018/PreprocessedData/Methylation/panos_methylation_data.rda",
                fdrThreshold=0.05,
                dmrPrintNum=5,
                outputDir="//isilon.c2b2.columbia.edu/ifs/archive/shares/bisr/Papapanou/Sept_2018/DifferentialMethylation/") {
  
  ### load libraries
  if(!require(limma, quietly = TRUE)) {
    if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
    BiocManager::install("limma")
    require(limma, quietly = TRUE)
  }
  if(!require(minfi, quietly = TRUE)) {
    if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
    BiocManager::install("minfi")
    require(minfi, quietly = TRUE)
  }
  if(!require(IlluminaHumanMethylationEPICanno.ilm10b2.hg19, quietly = TRUE)) {
    if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
    BiocManager::install("IlluminaHumanMethylationEPICanno.ilm10b2.hg19")
    require(IlluminaHumanMethylationEPICanno.ilm10b2.hg19, quietly = TRUE)
  }
  if(!require(DMRcate, quietly = TRUE)) {
    if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
    BiocManager::install("DMRcate")
    require(DMRcate, quietly = TRUE)
  }
  if(!require(missMethyl, quietly = TRUE)) {
    if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
    BiocManager::install("missMethyl")
    require(missMethyl, quietly = TRUE)
  }
  
  ### load data
  load(inputFilelPath)
  
  
  ### get Illumina EPIC hg19 annotation data
  anno <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)
  anno <- anno[match(rownames(m), anno$Name),
               c(1:4,12:19,24:ncol(anno))]
  
  
  ### get factor attributes for design matrix
  sample_type <- factor(sample_info$Phenotype, levels = c("Healthy", "Gingivitis", "Periodontitis"))
  sample_group <- factor(paste0("G",sample_info$SentrixBarcode_A))
  sample_batch <- factor(sample_info$Batch)
  
  
  ### make a design matrix
  design <- model.matrix(~0+sample_type+sample_batch)
  colnames(design) <- c(levels(sample_type), levels(sample_batch)[-1])
  
  
  ### fir the linear model with m value
  fit <- lmFit(m, design)
  
  
  ### make a contrast matrix
  contrastMat <- makeContrasts(Gingivitis-Healthy,
                               Periodontitis-Healthy,
                               Periodontitis-Gingivitis,
                               levels = design)
  
  
  ### fit the contrasts
  fit2 <- contrasts.fit(fit, contrastMat)
  fit2 <- eBayes(fit2)
  
  
  ### get the differentially methylated CpGs
  result1 <- topTable(fit2, coef = 1, adjust.method="BH", number=Inf, genelist = anno)
  result2 <- topTable(fit2, coef = 2, adjust.method="BH", number=Inf, genelist = anno)
  result3 <- topTable(fit2, coef = 3, adjust.method="BH", number=Inf, genelist = anno)
  
  
  ### write the DMP results
  write.table(result1, file = paste0(outputDir, "DMP_Gingivitis_Healthy.txt"), sep = "\t", row.names = FALSE)
  write.table(result2, file = paste0(outputDir, "DMP_Periodontitis_Healthy.txt"), sep = "\t", row.names = FALSE)
  write.table(result3, file = paste0(outputDir, "DMP_Periodontitis_Gingivitis.txt"), sep = "\t", row.names = FALSE)
  
  ### write the DMP (FDR < fdrThreshold) results
  if(length(which(result1$adj.P.Val < fdrThreshold)) > 0) {
    write.table(result1[which(result1$adj.P.Val < fdrThreshold),], file = paste0(outputDir, "DMP_Gingivitis_Healthy_", fdrThreshold, ".txt"), sep = "\t", row.names = FALSE)
  }
  if(length(which(result2$adj.P.Val < fdrThreshold)) > 0) {
    write.table(result2[which(result2$adj.P.Val < fdrThreshold),], file = paste0(outputDir, "DMP_Periodontitis_Healthy_", fdrThreshold, ".txt"), sep = "\t", row.names = FALSE)
  }
  if(length(which(result3$adj.P.Val < fdrThreshold)) > 0) {
    write.table(result3[which(result3$adj.P.Val < fdrThreshold),], file = paste0(outputDir, "DMP_Periodontitis_Gingivitis_", fdrThreshold, ".txt"), sep = "\t", row.names = FALSE)
  }
  
  ### annotation for DMR(Differentially Methylated Region)s
  ### Gingivitis - Healthy
  myAnno1 <- cpg.annotate(object = m, datatype = "array", what = "M",
                         arraytype = "EPIC", fdr = fdrThreshold,
                         analysis.type = "differential", design = design,
                         contrasts = TRUE, cont.matrix = contrastMat,
                         coef = "Gingivitis - Healthy")
  ### Periodontitis - Healthy
  myAnno2 <- cpg.annotate(object = m, datatype = "array", what = "M",
                          arraytype = "EPIC", fdr = fdrThreshold,
                          analysis.type = "differential", design = design,
                          contrasts = TRUE, cont.matrix = contrastMat,
                          coef = "Periodontitis - Healthy")
  ### Periodontitis - Gingivitis
  myAnno3 <- cpg.annotate(object = m, datatype = "array", what = "M",
                          arraytype = "EPIC", fdr = fdrThreshold,
                          analysis.type = "differential", design = design,
                          contrasts = TRUE, cont.matrix = contrastMat,
                          coef = "Periodontitis - Gingivitis")
  
  
  ### get DMRs
  ### Gingivitis - Healthy
  if(length(which(myAnno1@ranges$is.sig == TRUE)) > 0) {
    ### get DMRs
    DMRs1 <- dmrcate(myAnno1, lambda=1000, C=2)
    ### extract ranges from DMRs
    results.ranges1 <- extractRanges(DMRs1, genome = "hg19")
    ### filter with Stouffer combined p-values
    if(length(results.ranges1) > 0) {
      results.ranges1 <- results.ranges1[which(results.ranges1$Stouffer < fdrThreshold),]
    }
    ### save the DMR info
    if(length(results.ranges1) > 0) {
      write.table(data.frame(results.ranges1), file = paste0(outputDir, "DMR_Gingivitis_Healthy.txt"),
                  sep = "\t", row.names = FALSE)
    }
  }
  ### Periodontitis - Healthy
  if(length(which(myAnno2@ranges$is.sig == TRUE)) > 0) {
    ### get DMRs
    DMRs2 <- dmrcate(myAnno2, lambda=1000, C=2)
    ### extract ranges from DMRs
    results.ranges2 <- extractRanges(DMRs2, genome = "hg19")
    ### filter with Stouffer combined p-values
    if(length(results.ranges2) > 0) {
      results.ranges2 <- results.ranges2[which(results.ranges2$Stouffer < fdrThreshold),]
    }
    ### save the DMR info
    if(length(results.ranges2) > 0) {
      write.table(data.frame(results.ranges2), file = paste0(outputDir, "DMR_Periodontitis_Healthy.txt"),
                  sep = "\t", row.names = FALSE)
    }
  }
  ### Periodontitis - Gingivitis
  if(length(which(myAnno3@ranges$is.sig == TRUE)) > 0) {
    ### get DMRs
    DMRs3 <- dmrcate(myAnno3, lambda=1000, C=2)
    ### extract ranges from DMRs
    results.ranges3 <- extractRanges(DMRs3, genome = "hg19")
    ### filter with Stouffer combined p-values
    if(length(results.ranges3) > 0) {
      results.ranges3 <- results.ranges3[which(results.ranges3$Stouffer < fdrThreshold),]
    }
    ### save the DMR info
    if(length(results.ranges3) > 0) {
      write.table(data.frame(results.ranges3), file = paste0(outputDir, "DMR_Periodontitis_Gingivitis.txt"),
                  sep = "\t", row.names = FALSE)
    }
  }
  
  
  ### set up parameters for a plot
  pheno <- c("orange", "green", "pink")
  names(pheno) <- levels(factor(sample_info$Phenotype))
  cols <- pheno[as.character(factor(sample_info$Phenotype))]
  
  
  ### draw DMR plots with top DMRs
  for(i in 1:dmrPrintNum) {
    ### Gingivitis - Healthy
    if(length(which(myAnno1@ranges$is.sig == TRUE)) >= i) {
      ### draw DMR plots
      idx <- union(which(sample_info$Phenotype == "Gingivitis"),
                   which(sample_info$Phenotype == "Healthy"))
      png(paste0(outputDir, "DMR", i, "_Gingivitis_Healthy.png"), width = 3600, height = 2400, res = 350)
      DMR.plot(ranges=results.ranges1, dmr=i, CpGs=m[,idx], phen.col=cols[idx], what = "M",
               arraytype = "EPIC", pch=16, toscale=TRUE, plotmedians=TRUE, 
               genome="hg19", samps=1:length(idx))
      dev.off()
    }
    ### Periodontitis - Healthy
    if(length(which(myAnno2@ranges$is.sig == TRUE)) >= i) {
      ### draw DMR plots
      idx <- union(which(sample_info$Phenotype == "Periodontitis"),
                   which(sample_info$Phenotype == "Healthy"))
      png(paste0(outputDir, "DMR", i, "_Periodontitis_Healthy.png"), width = 3600, height = 2400, res = 350)
      DMR.plot(ranges=results.ranges2, dmr=i, CpGs=m[,idx], phen.col=cols[idx], what = "M",
               arraytype = "EPIC", pch=16, toscale=TRUE, plotmedians=TRUE, 
               genome="hg19", samps=1:length(idx))
      dev.off()
    }
    ### Periodontitis - Gingivitis
    if(length(which(myAnno3@ranges$is.sig == TRUE)) >= i) {
      ### draw DMR plots
      idx <- union(which(sample_info$Phenotype == "Periodontitis"),
                   which(sample_info$Phenotype == "Gingivitis"))
      png(paste0(outputDir, "DMR", i, "_Periodontitis_Gingivitis.png"), width = 3600, height = 2400, res = 350)
      DMR.plot(ranges=results.ranges3, dmr=i, CpGs=m[,idx], phen.col=cols[idx], what = "M",
               arraytype = "EPIC", pch=16, toscale=TRUE, plotmedians=TRUE, 
               genome="hg19", samps=1:length(idx))
      dev.off()
    }
  }
  
  
  ### pathway analysis with genes associated with the DMPs
  ### Gingivitis - Healthy
  sigCpGs1 <- result1$Name[result1$adj.P.Val < fdrThreshold]
  try (
    if(length(sigCpGs1) > 0) {
      ### GO
      gst_go1 <- gometh(sig.cpg = sigCpGs1, all.cpg = result1$Name, plot.bias = FALSE,
                    array.type = "EPIC", collection = "GO")
      gst_go1 <- gst_go1[order(gst_go1$FDR),]
      gst_go1 <- gst_go1[which(gst_go1$FDR < fdrThreshold),]
      ### KEGG
      gst_kegg1 <- gometh(sig.cpg = sigCpGs1, all.cpg = result1$Name, plot.bias = FALSE,
                        array.type = "EPIC", collection = "KEGG")
      gst_kegg1 <- gst_kegg1[order(gst_kegg1$FDR),]
      gst_kegg1 <- gst_kegg1[which(gst_kegg1$FDR < fdrThreshold),]
      ### write out the pathway results
      if(nrow(gst_go1) > 0) {
        write.table(cbind(ID=rownames(gst_go1), gst_go1), file = paste0(outputDir, "GO_DMP_Gingivitis_Healthy.txt"),
                    sep = "\t", row.names = FALSE)
      }
      if(nrow(gst_kegg1) > 0) {
        write.table(cbind(ID=rownames(gst_kegg1), gst_kegg1), file = paste0(outputDir, "KEGG_DMP_Gingivitis_Healthy.txt"),
                    sep = "\t", row.names = FALSE)
      }
    }
  , silent = TRUE)
  ### Periodontitis - Healthy
  sigCpGs2 <- result2$Name[result2$adj.P.Val < fdrThreshold]
  try (
    if(length(sigCpGs2) > 0) {
      ### GO
      gst_go2 <- gometh(sig.cpg = sigCpGs2, all.cpg = result2$Name, plot.bias = FALSE,
                    array.type = "EPIC", collection = "GO")
      gst_go2 <- gst_go2[order(gst_go2$FDR),]
      gst_go2 <- gst_go2[which(gst_go2$FDR < fdrThreshold),]
      ### KEGG
      gst_kegg2 <- gometh(sig.cpg = sigCpGs2, all.cpg = result2$Name, plot.bias = FALSE,
                        array.type = "EPIC", collection = "KEGG")
      gst_kegg2 <- gst_kegg2[order(gst_kegg2$FDR),]
      gst_kegg2 <- gst_kegg2[which(gst_kegg2$FDR < fdrThreshold),]
      ### write out the pathway results
      if(nrow(gst_go2) > 0) {
        write.table(cbind(ID=rownames(gst_go2), gst_go2), file = paste0(outputDir, "GO_DMP_Periodontitis_Healthy.txt"),
                    sep = "\t", row.names = FALSE)
      }
      if(nrow(gst_kegg2) > 0) {
        write.table(cbind(ID=rownames(gst_kegg2), gst_kegg2), file = paste0(outputDir, "KEGG_DMP_Periodontitis_Healthy.txt"),
                    sep = "\t", row.names = FALSE)
      }
    }
  , silent = TRUE)
  ### Periodontitis - Gingivitis
  sigCpGs3 <- result3$Name[result3$adj.P.Val < fdrThreshold]
  try (
    if(length(sigCpGs3) > 0) {
      ### GO
      gst_go3 <- gometh(sig.cpg = sigCpGs3, all.cpg = result3$Name, plot.bias = FALSE,
                    array.type = "EPIC", collection = "GO")
      gst_go3 <- gst_go3[order(gst_go3$FDR),]
      gst_go3 <- gst_go3[which(gst_go3$FDR < fdrThreshold),]
      ### KEGG
      gst_kegg3 <- gometh(sig.cpg = sigCpGs3, all.cpg = result3$Name, plot.bias = FALSE,
                        array.type = "EPIC", collection = "KEGG")
      gst_kegg3 <- gst_kegg3[order(gst_kegg3$FDR),]
      gst_kegg3 <- gst_kegg3[which(gst_kegg3$FDR < fdrThreshold),]
      ### write out the pathway results
      if(nrow(gst_go3) > 0) {
        write.table(cbind(ID=rownames(gst_go3), gst_go3), file = paste0(outputDir, "GO_DMP_Periodontitis_Gingivitis.txt"),
                    sep = "\t", row.names = FALSE)
      }
      if(nrow(gst_kegg3) > 0) {
        write.table(cbind(ID=rownames(gst_kegg3), gst_kegg3), file = paste0(outputDir, "KEGG_DMP_Periodontitis_Gingivitis.txt"),
                    sep = "\t", row.names = FALSE)
      }
    }
  , silent = TRUE)
  
}
