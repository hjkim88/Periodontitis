###
#   File name : ClusteringAffy.R
#   Author    : Hyunjin Kim
#   Date      : Dec 3, 2018
#   Email     : hk2990@cumc.columbia.edu
#   Purpose   : Because the Affy data has too many DE genes, we want to know if the periodontitis samples
#               carry more features than just periodontal disease condition. Therefore, clustering
#               the Affy samples with RNA-Seq samples to see if they are clustered in a group.
#               And we also found additional periodontitis samples on GEO database.
#               We also would like to compare them to what we have.
#
#   Instruction
#               1. Source("ClusteringAffy.R")
#               2. Run the function "clusteringAffy" - specify the input files (gene expression of the Affy and the RNA-Seq) and output directory
#               3. The clustering results will be generated in the output directory
#
#   Example
#               > source("The_directory_of_ClusteringAffy.R/ClusteringAffy.R")
#               > clusteringAffy(affyPath="./data/Affy/panos_affy_data.rda",
#                                rnaseqPath="./data/RNA_Seq/panos_rna_seq_data.rda",
#                                gse79705ExpPath="./data/GSE79705/GSE79705_100718_HG18_opt_expr_RMA.calls.txt",
#                                gse79705InfoPath="./data/GSE79705/GSE79705_sample_info.txt",
#                                gse23586ExpPath="./data/GSE23586/GSE23586_series_matrix.txt",
#                                gse23586InfoPath="./data/GSE23586/GSE23586_sample_info.txt",
#                                gse27993ExpPath="./data/GSE27993/GSE27993_series_matrix.txt",
#                                gse27993InfoPath="./data/GSE27993/GSE27993_sample_info.txt",
#                                gse33774ExpPath="./data/GSE33774/GSE33774_series_matrix.txt",
#                                gse33774InfoPath="./data/GSE33774/GSE33774_sample_info.txt",
#                                outputDir="./results/Test/")
###

clusteringAffy <- function(affyPath="//isilon.c2b2.columbia.edu/ifs/archive/shares/bisr/Papapanou/Sept_2018/PreprocessedData/Affy/panos_affy_data.rda",
                           rnaseqPath="//isilon.c2b2.columbia.edu/ifs/archive/shares/bisr/Papapanou/Sept_2018/PreprocessedData/RNA_Seq/panos_rna_seq_data.rda",
                           gse79705ExpPath="//isilon.c2b2.columbia.edu/ifs/archive/shares/bisr/Papapanou/Sept_2018/OriginalData/GSE79705/GSE79705_100718_HG18_opt_expr_RMA.calls.txt",
                           gse79705InfoPath="//isilon.c2b2.columbia.edu/ifs/archive/shares/bisr/Papapanou/Sept_2018/OriginalData/GSE79705/GSE79705_sample_info.txt",
                           gse23586ExpPath="//isilon.c2b2.columbia.edu/ifs/archive/shares/bisr/Papapanou/Sept_2018/OriginalData/GSE23586/GSE23586_series_matrix.txt",
                           gse23586InfoPath="//isilon.c2b2.columbia.edu/ifs/archive/shares/bisr/Papapanou/Sept_2018/OriginalData/GSE23586/GSE23586_sample_info.txt",
                           gse27993ExpPath="//isilon.c2b2.columbia.edu/ifs/archive/shares/bisr/Papapanou/Sept_2018/OriginalData/GSE27993/GSE27993_series_matrix.txt",
                           gse27993InfoPath="//isilon.c2b2.columbia.edu/ifs/archive/shares/bisr/Papapanou/Sept_2018/OriginalData/GSE27993/GSE27993_sample_info.txt",
                           gse33774ExpPath="//isilon.c2b2.columbia.edu/ifs/archive/shares/bisr/Papapanou/Sept_2018/OriginalData/GSE33774/GSE33774_series_matrix.txt",
                           gse33774InfoPath="//isilon.c2b2.columbia.edu/ifs/archive/shares/bisr/Papapanou/Sept_2018/OriginalData/GSE33774/GSE33774_sample_info.txt",
                           outputDir="//isilon.c2b2.columbia.edu/ifs/archive/shares/bisr/Papapanou/Sept_2018/Test/") {
  
  ### load library
  if(!require(sva, quietly = TRUE)) {
    source("https://bioconductor.org/biocLite.R")
    biocLite("sva")
    require(sva, quietly = TRUE)
  }
  if(!require(hgu133plus2.db, quietly = TRUE)) {
    source("https://bioconductor.org/biocLite.R")
    biocLite("hgu133plus2.db")
    require(hgu133plus2.db, quietly = TRUE)
  }
  if(!require(hugene10sttranscriptcluster.db, quietly = TRUE)) {
    source("https://bioconductor.org/biocLite.R")
    biocLite("hugene10sttranscriptcluster.db")
    require(hugene10sttranscriptcluster.db, quietly = TRUE)
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
  normCnt <- ComBat(dat = as.matrix(normCnt), batch = sampleInfo$Batch, mod = model.matrix(~Phenotype, data = sampleInfo))
  
  
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
  combined_data_cc <- ComBat(dat = as.matrix(combined_data), batch = combined_sample_info$Chip,
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
  
  
  ####################################
  ### comparison with GEO datasets ###
  ####################################
  
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
  
  
  ### GSE79705
  
  ### load data
  gse79705_exp <- read.table(file = gse79705ExpPath, header = TRUE, sep = "\t", row.names = 1,
                             stringsAsFactors = FALSE, check.names = FALSE)
  gse79705_sample_info <- read.table(file = gse79705InfoPath, header = TRUE, sep = "\t", row.names = 1,
                                     stringsAsFactors = FALSE, check.names = FALSE)
  
  
  ### pre-process the data
  gse79705_exp <- gse79705_exp[,-1]
  gse79705_exp <- data.frame(Gene_Symbol=sapply(gse79705_exp$GENE_INFO, function(x) {
    temp <- strsplit(x, split = " ", fixed = TRUE)[[1]][8]
    return(substr(temp, 1, nchar(temp)-1))
  }, USE.NAMES = FALSE), gse79705_exp, stringsAsFactors = FALSE, check.names = FALSE)
  gse79705_exp <- gse79705_exp[,-which(colnames(gse79705_exp) == "GENE_INFO")]
  
  
  ### DE analysis - periodontitis vs healthy
  deresult_gse79705 <- limmaWithComparisons(normCnt = gse79705_exp[,-1],
                                            grp = gse79705_sample_info$Phenotype,
                                            exp_class = "Periodontitis",
                                            ctrl_class = "Healthy",
                                            bat_eff = NULL)
  
  
  ### write out the DE result
  write.table(data.frame(Probe=rownames(deresult_gse79705), deresult_gse79705,
                         stringsAsFactors = FALSE, check.names = FALSE),
              file = paste0(outputDir, "limma_gse79705_periodontitis_vs_healthy.txt"),
              sep = "\t", row.names = FALSE)
  
  
  ### sum up for any one gene: probe -> gene mapping, one gene - one row
  unique_genes <- unique(gse79705_exp$Gene_Symbol)
  rIdx <- union(which(is.na(unique_genes)), which(unique_genes == "N/A"))
  if(length(rIdx) > 0) {
    unique_genes <- unique_genes[-rIdx]
  }
  new_gse79705_exp <- NULL
  for(i in 1:length(unique_genes)) {
    temp <- gse79705_exp[which(gse79705_exp$Gene_Symbol == unique_genes[i]), -1]
    new_gse79705_exp <- rbind(new_gse79705_exp, apply(temp, 2, sum))
  }
  rownames(new_gse79705_exp) <- unique_genes
  new_gse79705_exp <- as.data.frame(new_gse79705_exp)
  
  
  ### combine new affy and new gse79705
  common_genes <- intersect(rownames(new_affy_norm_ge), rownames(new_gse79705_exp))
  combined_affy_gse79705 <- data.frame(scale(new_affy_norm_ge[common_genes,]),
                                       scale(new_gse79705_exp[common_genes,]),
                                       stringsAsFactors = FALSE, check.names = FALSE)
  combined_affy_gse79705_sample_info <- data.frame(Phenotype=c(affy_norm_ge_sample_info$Condition,
                                                               gse79705_sample_info$Phenotype),
                                                   Chip=c(rep("Affy", ncol(new_affy_norm_ge)),
                                                          rep("GSE79705", ncol(new_gse79705_exp))),
                                                   stringsAsFactors = FALSE, check.names = FALSE)
  combined_affy_gse79705_sample_info[which(combined_affy_gse79705_sample_info$Phenotype == "Affected"), "Phenotype"] <- "Periodontitis"
  combined_affy_gse79705_sample_info[which(combined_affy_gse79705_sample_info$Phenotype == "Control"), "Phenotype"] <- "Healthy"
  rownames(combined_affy_gse79705_sample_info) <- colnames(combined_affy_gse79705)
  
  
  ### PCA plot - lebeling with the phenotype info
  pca_plot(normalizedMat = combined_affy_gse79705,
           grp = combined_affy_gse79705_sample_info$Phenotype,
           title = "pca_combined_affy_gse79705_phenotype",
           filePath = paste0(outputDir, "pca_combined_affy_gse79705_phenotype.png"))
  
  
  ### PCA plot - lebeling with the chip info
  pca_plot(normalizedMat = combined_affy_gse79705,
           grp = combined_affy_gse79705_sample_info$Chip,
           title = "pca_combined_affy_gse79705_chip",
           filePath = paste0(outputDir, "pca_combined_affy_gse79705_chip.png"))
  
  
  ### chip effect correction
  combined_data_cc <- ComBat(dat = as.matrix(combined_affy_gse79705), batch = combined_affy_gse79705_sample_info$Chip,
                             mod = model.matrix(~combined_affy_gse79705_sample_info$Phenotype))
  
  
  ### PCA plot - lebeling with the phenotype info - After chip effect correction
  pca_plot(normalizedMat = combined_data_cc,
           grp = combined_affy_gse79705_sample_info$Phenotype,
           title = "pca_combined_affy_gse79705_phenotype_chip_effect_corrected",
           filePath = paste0(outputDir, "pca_combined_affy_gse79705_phenotype_chip_effect_corrected.png"))
  
  
  ### PCA plot - lebeling with the chip info - After chip effect correction
  pca_plot(normalizedMat = combined_data_cc,
           grp = combined_affy_gse79705_sample_info$Chip,
           title = "pca_combined_affy_gse79705_chip_chip_effect_corrected",
           filePath = paste0(outputDir, "pca_combined_affy_gse79705_chip_chip_effect_corrected.png"))
  
  
  ### GSE23586
  
  ### load data
  gse23586_exp <- read.table(file = gse23586ExpPath, header = TRUE, sep = "\t", row.names = 1,
                             stringsAsFactors = FALSE, check.names = FALSE)
  gse23586_sample_info <- read.table(file = gse23586InfoPath, header = TRUE, sep = "\t", row.names = 1,
                                     stringsAsFactors = FALSE, check.names = FALSE)
  
  
  ### annotate gene symbols to the data
  ### convert microarray probes to gene symbols
  ### hgu133plus2
  ### get gene symbol mapping info
  x <- hgu133plus2SYMBOL
  mapped_probes <- mappedkeys(x)
  xx <- as.list(x[mapped_probes])
  
  
  ### get common probes
  common_probes <- intersect(rownames(gse23586_exp), names(xx))
  
  
  ### keep common probes only
  gse23586_exp <- gse23586_exp[common_probes,]
  
  
  ### annotate gene symbols
  gse23586_exp <- data.frame(Gene_Symbol=sapply(xx[common_probes], function(y) {
    return(y[1])
  }),
  gse23586_exp, stringsAsFactors = FALSE, check.names = FALSE)
  
  
  ### DE analysis - periodontitis vs healthy
  deresult_gse23586 <- limmaWithComparisons(normCnt = gse23586_exp[,-1],
                                            grp = gse23586_sample_info$Phenotype,
                                            exp_class = "Periodontitis",
                                            ctrl_class = "Healthy",
                                            bat_eff = gse23586_sample_info$Patient)
  
  
  ### write out the DE result
  write.table(data.frame(Probe=rownames(deresult_gse23586), deresult_gse23586,
                         stringsAsFactors = FALSE, check.names = FALSE),
              file = paste0(outputDir, "limma_gse23586_periodontitis_vs_healthy.txt"),
              sep = "\t", row.names = FALSE)
  
  
  ### combine affy and GSE23586 data
  common_probes <- intersect(rownames(affy_norm_ge), rownames(gse23586_exp))
  combined_affy_gse23586 <- data.frame(scale(affy_norm_ge[common_probes,-1]),
                                       scale(gse23586_exp[common_probes,-1]),
                                       stringsAsFactors = FALSE, check.names = FALSE)
  combined_affy_gse23586_sample_info <- data.frame(Phenotype=c(affy_norm_ge_sample_info$Condition,
                                                               gse23586_sample_info$Phenotype),
                                                   Chip=c(rep("Affy", ncol(affy_norm_ge)-1),
                                                          rep("GSE23586", ncol(gse23586_exp)-1)),
                                                   stringsAsFactors = FALSE, check.names = FALSE)
  combined_affy_gse23586_sample_info[which(combined_affy_gse23586_sample_info$Phenotype == "Affected"), "Phenotype"] <- "Periodontitis"
  combined_affy_gse23586_sample_info[which(combined_affy_gse23586_sample_info$Phenotype == "Control"), "Phenotype"] <- "Healthy"
  rownames(combined_affy_gse23586_sample_info) <- colnames(combined_affy_gse23586)
  
  
  ### PCA plot - lebeling with the phenotype info
  pca_plot(normalizedMat = combined_affy_gse23586,
           grp = combined_affy_gse23586_sample_info$Phenotype,
           title = "pca_combined_affy_gse23586_phenotype",
           filePath = paste0(outputDir, "pca_combined_affy_gse23586_phenotype.png"))
  
  
  ### PCA plot - lebeling with the chip info
  pca_plot(normalizedMat = combined_affy_gse23586,
           grp = combined_affy_gse23586_sample_info$Chip,
           title = "pca_combined_affy_gse23586_chip",
           filePath = paste0(outputDir, "pca_combined_affy_gse23586_chip.png"))
  
  
  ### batch effect correction
  combined_data_bc <- ComBat(dat = as.matrix(combined_affy_gse23586), batch = combined_affy_gse23586_sample_info$Chip,
                             mod = model.matrix(~combined_affy_gse23586_sample_info$Phenotype))
  
  
  ### PCA plot - lebeling with the phenotype info - After chip effect correction
  pca_plot(normalizedMat = combined_data_bc,
           grp = combined_affy_gse23586_sample_info$Phenotype,
           title = "pca_combined_affy_gse23586_phenotype_chip_effect_corrected",
           filePath = paste0(outputDir, "pca_combined_affy_gse23586_phenotype_chip_effect_corrected.png"))
  
  
  ### PCA plot - lebeling with the chip info - After chip effect correction
  pca_plot(normalizedMat = combined_data_bc,
           grp = combined_affy_gse23586_sample_info$Chip,
           title = "pca_combined_affy_gse23586_chip_chip_effect_corrected",
           filePath = paste0(outputDir, "pca_combined_affy_gse23586_chip_chip_effect_corrected.png"))
  
  
  ### GSE27993
  
  ### load data
  gse27993_exp <- read.table(file = gse27993ExpPath, header = TRUE, sep = "\t", row.names = 1,
                             stringsAsFactors = FALSE, check.names = FALSE)
  gse27993_sample_info <- read.table(file = gse27993InfoPath, header = TRUE, sep = "\t", row.names = 1,
                                     stringsAsFactors = FALSE, check.names = FALSE)
  
  
  ### annotate gene symbols to the data
  ### convert microarray probes to gene symbols
  ### [HuGene-1_0-st] Affymetrix Human Gene 1.0 ST Array
  ### get gene symbol mapping info
  x <- hugene10sttranscriptclusterSYMBOL
  mapped_probes <- mappedkeys(x)
  xx <- as.list(x[mapped_probes])
  
  
  ### get common probes
  common_probes <- intersect(rownames(gse27993_exp), names(xx))
  
  
  ### keep common probes only
  gse27993_exp <- gse27993_exp[common_probes,]
  
  
  ### annotate gene symbols
  gse27993_exp <- data.frame(Gene_Symbol=sapply(xx[common_probes], function(y) {
    return(y[1])
  }),
  gse27993_exp, stringsAsFactors = FALSE, check.names = FALSE)
  
  
  ### DE analysis - periodontitis vs healthy
  deresult_gse27993 <- limmaWithComparisons(normCnt = gse27993_exp[,-1],
                                            grp = gse27993_sample_info$Phenotype,
                                            exp_class = "Periodontitis",
                                            ctrl_class = "Healthy",
                                            bat_eff = NULL)
  
  
  ### write out the DE result
  write.table(data.frame(Probe=rownames(deresult_gse27993), deresult_gse27993,
                         stringsAsFactors = FALSE, check.names = FALSE),
              file = paste0(outputDir, "limma_gse27993_periodontitis_vs_healthy.txt"),
              sep = "\t", row.names = FALSE)
  
  
  ### sum up for any one gene: probe -> gene mapping, one gene - one row
  unique_genes <- unique(gse27993_exp$Gene_Symbol)
  new_gse27993_exp <- NULL
  for(i in 1:length(unique_genes)) {
    temp <- gse27993_exp[which(gse27993_exp$Gene_Symbol == unique_genes[i]), -1]
    new_gse27993_exp <- rbind(new_gse27993_exp, apply(temp, 2, sum))
  }
  rownames(new_gse27993_exp) <- unique_genes
  new_gse27993_exp <- as.data.frame(new_gse27993_exp)
  
  
  ### combine new affy and new gse27993
  common_genes <- intersect(rownames(new_affy_norm_ge), rownames(new_gse27993_exp))
  combined_affy_gse27993 <- data.frame(scale(new_affy_norm_ge[common_genes,]),
                                       scale(new_gse27993_exp[common_genes,]),
                                       stringsAsFactors = FALSE, check.names = FALSE)
  combined_affy_gse27993_sample_info <- data.frame(Phenotype=c(affy_norm_ge_sample_info$Condition,
                                                               gse27993_sample_info$Phenotype),
                                                   Chip=c(rep("Affy", ncol(new_affy_norm_ge)),
                                                          rep("GSE27993", ncol(new_gse27993_exp))),
                                                   stringsAsFactors = FALSE, check.names = FALSE)
  combined_affy_gse27993_sample_info[which(combined_affy_gse27993_sample_info$Phenotype == "Affected"), "Phenotype"] <- "Periodontitis"
  combined_affy_gse27993_sample_info[which(combined_affy_gse27993_sample_info$Phenotype == "Control"), "Phenotype"] <- "Healthy"
  rownames(combined_affy_gse27993_sample_info) <- colnames(combined_affy_gse27993)
  
  
  ### PCA plot - lebeling with the phenotype info
  pca_plot(normalizedMat = combined_affy_gse27993,
           grp = combined_affy_gse27993_sample_info$Phenotype,
           title = "pca_combined_affy_gse27993_phenotype",
           filePath = paste0(outputDir, "pca_combined_affy_gse27993_phenotype.png"))
  
  
  ### PCA plot - lebeling with the chip info
  pca_plot(normalizedMat = combined_affy_gse27993,
           grp = combined_affy_gse27993_sample_info$Chip,
           title = "pca_combined_affy_gse27993_chip",
           filePath = paste0(outputDir, "pca_combined_affy_gse27993_chip.png"))
  
  
  ### chip effect correction
  combined_data_cc <- ComBat(dat = as.matrix(combined_affy_gse27993), batch = combined_affy_gse27993_sample_info$Chip,
                             mod = model.matrix(~combined_affy_gse27993_sample_info$Phenotype))
  
  
  ### PCA plot - lebeling with the phenotype info - After chip effect correction
  pca_plot(normalizedMat = combined_data_cc,
           grp = combined_affy_gse27993_sample_info$Phenotype,
           title = "pca_combined_affy_gse27993_phenotype_chip_effect_corrected",
           filePath = paste0(outputDir, "pca_combined_affy_gse27993_phenotype_chip_effect_corrected.png"))
  
  
  ### PCA plot - lebeling with the chip info - After chip effect correction
  pca_plot(normalizedMat = combined_data_cc,
           grp = combined_affy_gse27993_sample_info$Chip,
           title = "pca_combined_affy_gse27993_chip_chip_effect_corrected",
           filePath = paste0(outputDir, "pca_combined_affy_gse27993_chip_chip_effect_corrected.png"))
  
  
  ### GSE27993
  
  ### load data
  gse33774_exp <- read.table(file = gse33774ExpPath, header = TRUE, sep = "\t", row.names = 1,
                             stringsAsFactors = FALSE, check.names = FALSE)
  gse33774_sample_info <- read.table(file = gse33774InfoPath, header = TRUE, sep = "\t", row.names = 1,
                                     stringsAsFactors = FALSE, check.names = FALSE)
  
  
  ### annotate gene symbols to the data
  ### convert microarray probes to gene symbols
  ### [HuGene-1_0-st] Affymetrix Human Gene 1.0 ST Array
  ### get gene symbol mapping info
  x <- hugene10sttranscriptclusterSYMBOL
  mapped_probes <- mappedkeys(x)
  xx <- as.list(x[mapped_probes])
  
  
  ### get common probes
  common_probes <- intersect(rownames(gse33774_exp), names(xx))
  
  
  ### keep common probes only
  gse33774_exp <- gse33774_exp[common_probes,]
  
  
  ### annotate gene symbols
  gse33774_exp <- data.frame(Gene_Symbol=sapply(xx[common_probes], function(y) {
    return(y[1])
  }),
  gse33774_exp, stringsAsFactors = FALSE, check.names = FALSE)
  
  
  ### DE analysis - parodontitis vs healthy
  deresult_gse33774 <- limmaWithComparisons(normCnt = gse33774_exp[,-1],
                                            grp = gse33774_sample_info$Phenotype,
                                            exp_class = "Parodontitis",
                                            ctrl_class = "Healthy",
                                            bat_eff = NULL)
  
  
  ### write out the DE result
  write.table(data.frame(Probe=rownames(deresult_gse33774), deresult_gse33774,
                         stringsAsFactors = FALSE, check.names = FALSE),
              file = paste0(outputDir, "limma_gse33774_parodontitis_vs_healthy.txt"),
              sep = "\t", row.names = FALSE)
  
  
  ### sum up for any one gene: probe -> gene mapping, one gene - one row
  unique_genes <- unique(gse33774_exp$Gene_Symbol)
  new_gse33774_exp <- NULL
  for(i in 1:length(unique_genes)) {
    temp <- gse33774_exp[which(gse33774_exp$Gene_Symbol == unique_genes[i]), -1]
    new_gse33774_exp <- rbind(new_gse33774_exp, apply(temp, 2, sum))
  }
  rownames(new_gse33774_exp) <- unique_genes
  new_gse33774_exp <- as.data.frame(new_gse33774_exp)
  
  
  ### combine new affy and new gse33774
  common_genes <- intersect(rownames(new_affy_norm_ge), rownames(new_gse33774_exp))
  combined_affy_gse33774 <- data.frame(scale(new_affy_norm_ge[common_genes,]),
                                       scale(new_gse33774_exp[common_genes,]),
                                       stringsAsFactors = FALSE, check.names = FALSE)
  combined_affy_gse33774_sample_info <- data.frame(Phenotype=c(affy_norm_ge_sample_info$Condition,
                                                               gse33774_sample_info$Phenotype),
                                                   Chip=c(rep("Affy", ncol(new_affy_norm_ge)),
                                                          rep("GSE33774", ncol(new_gse33774_exp))),
                                                   stringsAsFactors = FALSE, check.names = FALSE)
  combined_affy_gse33774_sample_info[which(combined_affy_gse33774_sample_info$Phenotype == "Affected"), "Phenotype"] <- "Periodontitis"
  combined_affy_gse33774_sample_info[which(combined_affy_gse33774_sample_info$Phenotype == "Control"), "Phenotype"] <- "Healthy"
  rownames(combined_affy_gse33774_sample_info) <- colnames(combined_affy_gse33774)
  
  
  ### PCA plot - lebeling with the phenotype info
  pca_plot(normalizedMat = combined_affy_gse33774,
           grp = combined_affy_gse33774_sample_info$Phenotype,
           title = "pca_combined_affy_gse33774_phenotype",
           filePath = paste0(outputDir, "pca_combined_affy_gse33774_phenotype.png"))
  
  
  ### PCA plot - lebeling with the chip info
  pca_plot(normalizedMat = combined_affy_gse33774,
           grp = combined_affy_gse33774_sample_info$Chip,
           title = "pca_combined_affy_gse33774_chip",
           filePath = paste0(outputDir, "pca_combined_affy_gse33774_chip.png"))
  
  
  ### chip effect correction
  combined_data_cc <- ComBat(dat = as.matrix(combined_affy_gse33774), batch = combined_affy_gse33774_sample_info$Chip,
                             mod = model.matrix(~combined_affy_gse33774_sample_info$Phenotype))
  
  
  ### PCA plot - lebeling with the phenotype info - After chip effect correction
  pca_plot(normalizedMat = combined_data_cc,
           grp = combined_affy_gse33774_sample_info$Phenotype,
           title = "pca_combined_affy_gse33774_phenotype_chip_effect_corrected",
           filePath = paste0(outputDir, "pca_combined_affy_gse33774_phenotype_chip_effect_corrected.png"))
  
  
  ### PCA plot - lebeling with the chip info - After chip effect correction
  pca_plot(normalizedMat = combined_data_cc,
           grp = combined_affy_gse33774_sample_info$Chip,
           title = "pca_combined_affy_gse33774_chip_chip_effect_corrected",
           filePath = paste0(outputDir, "pca_combined_affy_gse33774_chip_chip_effect_corrected.png"))
  
}




