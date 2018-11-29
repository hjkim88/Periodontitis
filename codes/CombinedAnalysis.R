###
#   File name : CombinedAnalysis.R
#   Author    : Hyunjin Kim
#   Date      : Oct 14, 2018
#   Email     : hk2990@cumc.columbia.edu
#   Purpose   : Analyze using both Methylation and expression levels of genes
#               It is an analysis using both methylation data and RNA-Seq
#
#   To investigate differential methylation patterns associated with genes that were 
#   differentially expressed in the three periodontal phenotypes.
#
#   Instruction
#               1. Source("CombinedAnalysis.R")
#               2. Run the function "combAnalysis" - specify the input files (Differential) of methylation and RNA-Seq and output directory
#               3. Various results will be generated in the output path
#
#   Example
#               > source("The_directory_of_CombinedAnalysis.R/CombinedAnalysis.R")
#               > combAnalysis(methylPath="./results/DMA/DMP_Periodontitis_Healthy.txt",
#                              methylRegPath="./results/DMA/DMR_Periodontitis_Healthy.txt",
#                              methylRDAPath="./data/Methylation/panos_methylation_data.rda",
#                              rnaseqPath="./results/DEA/DE_Results_Periodontitis_vs_Healthy.xlsx",
#                              rnaseqRDAPath="./data/RNA_Seq/panos_rna_seq_data.rda",
#                              fdrThreshold=0.05,
#                              outputDir="./results/Combined/")
###

combAnalysis <- function(methylPath="//isilon.c2b2.columbia.edu/ifs/archive/shares/bisr/Papapanou/Sept_2018/DifferentialMethylation/DMP_Periodontitis_Healthy.txt",
                         methylRegPath="//isilon.c2b2.columbia.edu/ifs/archive/shares/bisr/Papapanou/Sept_2018/DifferentialMethylation/DMR_Periodontitis_Healthy.txt",
                         methylRDAPath="//isilon.c2b2.columbia.edu/ifs/archive/shares/bisr/Papapanou/Sept_2018/PreprocessedData/Methylation/panos_methylation_data.rda",
                         rnaseqPath="//isilon.c2b2.columbia.edu/ifs/archive/shares/bisr/Papapanou/Sept_2018/DifferentialExpression/DE_Results_Periodontitis_vs_Healthy.xlsx",
                         rnaseqRDAPath="//isilon.c2b2.columbia.edu/ifs/archive/shares/bisr/Papapanou/Sept_2018/PreprocessedData/RNA_Seq/panos_rna_seq_data.rda",
                         fdrThreshold=0.05,
                         outputDir="//isilon.c2b2.columbia.edu/ifs/archive/shares/bisr/Papapanou/Sept_2018/Combined/") {
  
  ### load libraries
  if(!require(xlsx)) {
    install.packages("xlsx")
    library(xlsx)
  }
  if(!require(ggplot2)) {
    install.packages("ggplot2")
    library(ggplot2)
  }
  if(!require(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)) {
    source("https://bioconductor.org/biocLite.R")
    biocLite("IlluminaHumanMethylationEPICanno.ilm10b2.hg19")
    library(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)
  }
  if(!require(FEM)) {
    source("https://bioconductor.org/biocLite.R")
    biocLite("FEM")
    library(FEM)
  }
  
  
  ### load data
  dmp <- read.table(file = methylPath, header = TRUE, sep = "\t", check.names = FALSE)
  dmr <- read.table(file = methylRegPath, header = TRUE, sep = "\t", check.names = FALSE)
  de <- read.xlsx2(file = rnaseqPath, sheetIndex = 1, stringsAsFactors = FALSE)
  load(methylRDAPath)
  load(rnaseqRDAPath)
  
  
  ### set row names
  rownames(dmp) <- dmp$Name
  rownames(de) <- de$Gene_Symbol
  
  
  ### Using all the data - not only using those from differential analyses
  
  ### create new column and put NAs
  de$CpG <- NA
  
  
  ### get DE gene associated CpG sites, put the info on the data
  for(i in 1:nrow(de)) {
    ### find indicies of differentially methylated CpGs associated with the DE genes
    idx <- grep(de$Gene_Symbol[i], dmp$GencodeCompV12_NAME)
    
    if(length(idx) > 0) {
      ### because "grep" can also detct part of full sentence, we need to filter one more
      temp <- as.character(dmp$GencodeCompV12_NAME[idx])
      temp <- sapply(temp, function(x) {
        strsplit(x, ";", fixed = TRUE)
      })
      
      ### get indicies that are real
      keep <- sapply(temp, function(x) {
        return(sum(as.integer(grepl(de$Gene_Symbol[i], x))) > 0)
      })
      idx <- idx[which(keep)]
      
      ### add the CpG info to the gene
      de$CpG[i] <- paste(dmp$Name[idx], collapse = ";")
    }
    
    if(i %% 1000 == 0) {
      writeLines(paste(i, "/", nrow(de)))
    }
  }
  
  
  ### create new column and put NAs
  de$MethylMean <- NA
  
  
  ### calculate mean FDR of CpGs
  for(i in 1:nrow(de)) {
    temp <- strsplit(de$CpG[i], ";")[[1]]
    mean <- 0
    for(j in 1:length(temp)) {
      mean <- mean + dmp[temp[j],"adj.P.Val"]
    }
    de$MethylMean[i] <- (mean / length(temp))
    
    if(i %% 1000 == 0) {
      writeLines(paste(i, "/", nrow(de)))
    }
  }
  
  
  ### create a correlation plot between DM CpGs and DE genes
  cor_data <- data.frame(Expression=de$padj, Methylation=de$MethylMean)
  idx <- sapply(cor_data, is.factor)
  cor_data[idx] <- lapply(cor_data[idx], function(x) as.numeric(as.character(x)))
  fName <- paste0("FDR_Correlation_", substr(basename(methylPath), 1, nchar(basename(methylPath))-4), "_all.png")
  ggplot(data = cor_data, aes(x=Expression, y=Methylation)) +
    geom_point(color = "black", size = 1) +
    labs(title=substr(fName, 1, nchar(fName)-4),
         subtitle=sprintf("P.Cor = %s, p-value = %s",
                          round(cor(cor_data$Expression, cor_data$Methylation, use = "pairwise.complete.obs"), 5),
                          signif(cor.test(cor_data$Expression, cor_data$Methylation)$p.value, 5))) +
    geom_smooth(method = lm, color="gray", se=FALSE) +
    theme_classic(base_size = 16)
  ggsave(filename = paste0(outputDir, fName), width = 10, height = 10)
  
  
  ### DMR
  
  ### a function to return an associated gene list vector of with a given region
  getGeneSet <- function(str) {
    ### split the string
    temp <- strsplit(as.character(str), ", ", fixed = TRUE)[[1]]
    ### remove promoter postfix
    temp <- sapply(temp, function(x) {
      return(substr(x, 1, nchar(x)-4))
    })
    
    return(temp)
  }
  
  
  ### make a data for a plot
  cor_data <- NULL
  for(i in 1:nrow(dmr)) {
    ### get DMR-associated gene list vector
    temp <- getGeneSet(dmr$overlapping.promoters[i])
    
    ### sometimes, there are redundant genes, so get an unique set
    temp <- unique(temp)
    
    ### a region may have more than one genes
    for(j in 1:length(temp)) {
      if(!is.na(de[temp[j],"padj"])) {
        cor_data <- rbind(cor_data, as.numeric(c(de[temp[j],"padj"], dmr$Stouffer[i])))
      }
    }
  }
  
  
  ### make a correlation plot between DM regions and DE genes
  fName <- paste0("FDR_Correlation_", substr(basename(methylRegPath), 1, nchar(basename(methylRegPath))-4), "_all.png")
  ggplot(data = as.data.frame(cor_data), aes(x=V1, y=V2)) +
    geom_point(color = "black", size = 1) +
    labs(title=substr(fName, 1, nchar(fName)-4),
         subtitle=sprintf("P.Cor = %s, p-value = %s",
                          round(cor(cor_data[,1], cor_data[,2], use = "pairwise.complete.obs"), 5),
                          signif(cor.test(cor_data[,1], cor_data[,2])$p.value, 5))) +
    xlab("Expression") +
    ylab("Methylation") +
    geom_smooth(method = lm, color="gray", se=FALSE) +
    theme_classic(base_size = 16)
  ggsave(filename = paste0(outputDir, fName), width = 10, height = 10)
  
  
  
  ### Using results with the cutoff
  
  ### filter the data with the given cutoff
  dmp <- dmp[which(dmp$adj.P.Val < fdrThreshold),]
  de <- de[which(de$padj < fdrThreshold),]
  
  
  ### create new column and put NAs
  de$CpG <- NA
  
  
  ### get DE gene associated CpG sites, put the info on the data
  for(i in 1:nrow(de)) {
    ### find indicies of differentially methylated CpGs associated with the DE genes
    idx <- grep(de$Gene_Symbol[i], dmp$GencodeCompV12_NAME)
    
    if(length(idx) > 0) {
      ### because "grep" can also detct part of full sentence, we need to filter one more
      temp <- as.character(dmp$GencodeCompV12_NAME[idx])
      temp <- sapply(temp, function(x) {
        strsplit(x, ";", fixed = TRUE)
      })
      
      ### get indicies that are real
      keep <- sapply(temp, function(x) {
        return(sum(as.integer(grepl(de$Gene_Symbol[i], x))) > 0)
      })
      idx <- idx[which(keep)]
      
      ### add the CpG info to the gene
      de$CpG[i] <- paste(dmp$Name[idx], collapse = ";")
    }
  }
  
  
  ### create new column and put NAs
  de$MethylMean <- NA
  
  
  ### calculate mean FDR of CpGs
  for(i in 1:nrow(de)) {
    temp <- strsplit(de$CpG[i], ";")[[1]]
    mean <- 0
    for(j in 1:length(temp)) {
      mean <- mean + dmp[temp[j],"adj.P.Val"]
    }
    de$MethylMean[i] <- (mean / length(temp))
  }
  
  
  ### create a correlation plot between DM CpGs and DE genes
  cor_data <- data.frame(Expression=de$padj, Methylation=de$MethylMean)
  idx <- sapply(cor_data, is.factor)
  cor_data[idx] <- lapply(cor_data[idx], function(x) as.numeric(as.character(x)))
  fName <- paste0("FDR_Correlation_", substr(basename(methylPath), 1, nchar(basename(methylPath))-4), "_", fdrThreshold, ".png")
  ggplot(data = cor_data, aes(x=Expression, y=Methylation)) +
    geom_point(color = "black", size = 1) +
    labs(title=substr(fName, 1, nchar(fName)-4),
         subtitle=sprintf("P.Cor = %s, p-value = %s",
                          round(cor(cor_data$Expression, cor_data$Methylation, use = "pairwise.complete.obs"), 5),
                          signif(cor.test(cor_data$Expression, cor_data$Methylation)$p.value, 5))) +
    geom_smooth(method = lm, color="gray", se=FALSE) +
    theme_classic(base_size = 16)
  ggsave(filename = paste0(outputDir, fName), width = 10, height = 10)
  
  
  ### DMR
  
  ### make a data for a plot
  cor_data <- NULL
  for(i in 1:nrow(dmr)) {
    ### get DMR-associated gene list vector
    temp <- getGeneSet(dmr$overlapping.promoters[i])
    
    ### sometimes, there are redundant genes, so get an unique set
    temp <- unique(temp)
    
    ### a region may have more than one genes
    for(j in 1:length(temp)) {
      if(!is.na(de[temp[j],"padj"])) {
        cor_data <- rbind(cor_data, as.numeric(c(de[temp[j],"padj"], dmr$Stouffer[i])))
      }
    }
  }
  
  
  ### make a correlation plot between DM regions and DE genes
  fName <- paste0("FDR_Correlation_", substr(basename(methylRegPath), 1, nchar(basename(methylRegPath))-4), "_", fdrThreshold, ".png")
  ggplot(data = as.data.frame(cor_data), aes(x=V1, y=V2)) +
    geom_point(color = "black", size = 1) +
    labs(title=substr(fName, 1, nchar(fName)-4),
         subtitle=sprintf("P.Cor = %s, p-value = %s",
                          round(cor(cor_data[,1], cor_data[,2], use = "pairwise.complete.obs"), 5),
                          signif(cor.test(cor_data[,1], cor_data[,2])$p.value, 5))) +
    xlab("Expression") +
    ylab("Methylation") +
    geom_smooth(method = lm, color="gray", se=FALSE) +
    theme_classic(base_size = 16)
  ggsave(filename = paste0(outputDir, fName), width = 10, height = 10)
  
  
  
  ### This time, not using FDR. using levels themselves
  
  ### get Illumina EPIC hg19 annotation data
  anno <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)
  anno <- anno[match(rownames(m), anno$Name),
               c(1:4,12:19,24:ncol(anno))]
  
  
  ### make grouped m values
  cor_data <- matrix(NA, nrow(m), 6)
  rownames(cor_data) <- rownames(m)
  colnames(cor_data) <- c("M_Periodontitis", "M_Gingivitis", "M_Healthy",
                          "Exp_Periodontitis", "Exp_Gingivitis", "Exp_Healthy")
  
  
  ### get mean m values of each phenotype
  cor_data[,1] <- apply(m[,which(sample_info$Phenotype == "Periodontitis")], 1, mean)
  cor_data[,2] <- apply(m[,which(sample_info$Phenotype == "Gingivitis")], 1, mean)
  cor_data[,3] <- apply(m[,which(sample_info$Phenotype == "Healthy")], 1, mean)
  
  
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
  normCnt <- normalizeRNASEQwithVST(rawCnt)
  
  
  ### a function to find associated gene list vector of a given CpG site
  getGeneSet2 <- function(cpg_name) {
    ### get genes
    temp <- anno[cpg_name,"GencodeCompV12_NAME"]
    
    ### split
    temp <- strsplit(temp, ";", fixed = TRUE)[[1]]
    
    ### get an unique set
    temp <- unique(temp)
    
    return(temp)
  }
  
  
  ### get mean exp values of each phenotype
  for(i in 1:nrow(cor_data)) {
    ### get a gene list vector associated with the given CpG
    genes <- getGeneSet2(rownames(cor_data)[i])
    
    ### get mean exp values
    temp <- length(which(rownames(normCnt) == genes))
    
    if(temp == 1) {
      ### if there's only one gene associated with the CpG
      cor_data[i,4] <- mean(as.numeric(normCnt[genes, which(sampleInfo$Phenotype == "Periodontitis")]), na.rm = TRUE)
      cor_data[i,5] <- mean(as.numeric(normCnt[genes, which(sampleInfo$Phenotype == "Gingivitis")]), na.rm = TRUE)
      cor_data[i,6] <- mean(as.numeric(normCnt[genes, which(sampleInfo$Phenotype == "Healthy")]), na.rm = TRUE)
    } else if(temp > 1) {
      ### if there are more than one genes associated with the CpG
      cor_data[i,4] <- mean(as.numeric(normCnt[genes[1], which(sampleInfo$Phenotype == "Periodontitis")]), na.rm = TRUE)
      cor_data[i,5] <- mean(as.numeric(normCnt[genes[1], which(sampleInfo$Phenotype == "Gingivitis")]), na.rm = TRUE)
      cor_data[i,6] <- mean(as.numeric(normCnt[genes[1], which(sampleInfo$Phenotype == "Healthy")]), na.rm = TRUE)
      for(j in 2:length(genes)) {
        cor_data <- rbind(cor_data, c(cor_data[i, 1:3],
                                      mean(as.numeric(normCnt[genes[j], which(sampleInfo$Phenotype == "Periodontitis")]), na.rm = TRUE),
                                      mean(as.numeric(normCnt[genes[j], which(sampleInfo$Phenotype == "Gingivitis")]), na.rm = TRUE),
                                      mean(as.numeric(normCnt[genes[j], which(sampleInfo$Phenotype == "Healthy")]), na.rm = TRUE)))
      }
    }
    
    ### print progress
    if(i %% 1000 == 0) {
      writeLines(paste(i, "/", nrow(cor_data)))
    }
  }
  
  
  ### make corrleation plots based on levels (not based on FDR)
  
  ### Periodontitis
  fName <- "Level_Correlation_Periodontitis.png"
  ggplot(data = as.data.frame(cor_data), aes(x=Exp_Periodontitis, y=M_Periodontitis)) +
    geom_point(color = "black", size = 1) +
    labs(title=substr(fName, 1, nchar(fName)-4),
         subtitle=sprintf("P.Cor = %s, p-value = %s",
                          round(cor(cor_data[,"M_Periodontitis"], cor_data[,"Exp_Periodontitis"], use = "pairwise.complete.obs"), 5),
                          signif(cor.test(cor_data[,"M_Periodontitis"], cor_data[,"Exp_Periodontitis"])$p.value, 5))) +
    xlab("Expression") +
    ylab("Methylation") +
    geom_smooth(method = lm, color="gray", se=FALSE) +
    theme_classic(base_size = 16)
  ggsave(filename = paste0(outputDir, fName), width = 10, height = 10)
  
  ### Gingivitis
  fName <- "Level_Correlation_Gingivitis.png"
  ggplot(data = as.data.frame(cor_data), aes(x=Exp_Gingivitis, y=M_Gingivitis)) +
    geom_point(color = "black", size = 1) +
    labs(title=substr(fName, 1, nchar(fName)-4),
         subtitle=sprintf("P.Cor = %s, p-value = %s",
                          round(cor(cor_data[,"M_Gingivitis"], cor_data[,"Exp_Gingivitis"], use = "pairwise.complete.obs"), 5),
                          signif(cor.test(cor_data[,"M_Gingivitis"], cor_data[,"Exp_Gingivitis"])$p.value, 5))) +
    xlab("Expression") +
    ylab("Methylation") +
    geom_smooth(method = lm, color="gray", se=FALSE) +
    theme_classic(base_size = 16)
  ggsave(filename = paste0(outputDir, fName), width = 10, height = 10)
  
  ### Healthy
  fName <- "Level_Correlation_Healthy.png"
  ggplot(data = as.data.frame(cor_data), aes(x=Exp_Healthy, y=M_Healthy)) +
    geom_point(color = "black", size = 1) +
    labs(title=substr(fName, 1, nchar(fName)-4),
         subtitle=sprintf("P.Cor = %s, p-value = %s",
                          round(cor(cor_data[,"M_Healthy"], cor_data[,"Exp_Healthy"], use = "pairwise.complete.obs"), 5),
                          signif(cor.test(cor_data[,"M_Healthy"], cor_data[,"Exp_Healthy"])$p.value, 5))) +
    xlab("Expression") +
    ylab("Methylation") +
    geom_smooth(method = lm, color="gray", se=FALSE) +
    theme_classic(base_size = 16)
  ggsave(filename = paste0(outputDir, fName), width = 10, height = 10)
  
}

