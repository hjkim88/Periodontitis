###
#   File name : LoadRawCounts.R
#   Author    : Hyunjin Kim
#   Date      : Oct 3, 2018
#   Email     : hk2990@cumc.columbia.edu
#   Purpose   : Load raw counts and sample info to make a combined data
#
#   Instruction
#               1. Source("LoadRawCounts.R")
#               2. Run the function "loadCounts" - specify the input directory (.counts.txt & sample info xlsx file) and output directory
#               3. The combined rawcounts data will be generated under the output directory
#
#   Example
#               > source("The_directory_of_LoadRawCounts.R/LoadRawCounts.R")
#               > loadCounts(inputDir="E:/Panos/RNA-Seq/",
#                            outputDir="./data/RNA_Seq/")
###

loadCounts <- function(inputDir="//isilon.c2b2.columbia.edu/ifs/archive/shares/bisr/Papapanou/Sept_2018/OriginalData/RNA-Seq/",
                       outputDir="//isilon.c2b2.columbia.edu/ifs/archive/shares/bisr/Papapanou/Sept_2018/PreprocessedData/RNA_Seq/") {
  
  ### load library
  if(!require(xlsx)) {
    install.packages("xlsx")
    library(xlsx)
  }
  
  
  ### get directories and files
  dirs <- list.files(inputDir)
  
  
  ### load sample info
  sampleInfo <- read.xlsx2(file = paste0(inputDir, dirs[which(endsWith(dirs, "V3.xlsx"))]),
                           sheetIndex = 1, startRow = 1)
  sampleInfo <- sampleInfo[,-ncol(sampleInfo)]
  rownames(sampleInfo) <- sampleInfo$RNA.seq.ID
  
  
  ### get raw counts
  rawCnt <- NULL
  sampleInfo$Phenotype <- NA
  
  ### healthy
  f <- list.files(paste0(inputDir, "healthy"))
  f <- f[which(endsWith(f, "counts.txt"))]
  for(i in 1:length(f)) {
    temp <- read.table(paste0(inputDir, "healthy/", f[i]), row.names = 1, header = TRUE)
    temp <- temp[order(rownames(temp)),,drop=FALSE]
    if(is.null(rawCnt)) {
      rawCnt <- temp
    } else {
      rawCnt <- cbind(rawCnt, temp)
    }
  }
  ### add phenotype to sample info
  sampleInfo[colnames(rawCnt)[(ncol(rawCnt)-length(f)+1):ncol(rawCnt)],"Phenotype"] <- "Healthy"
  
  ### gingivitis
  f <- list.files(paste0(inputDir, "gingivitis"))
  f <- f[which(endsWith(f, "counts.txt"))]
  for(i in 1:length(f)) {
    temp <- read.table(paste0(inputDir, "gingivitis/", f[i]), row.names = 1, header = TRUE)
    temp <- temp[order(rownames(temp)),,drop=FALSE]
    if(is.null(rawCnt)) {
      rawCnt <- temp
    } else {
      rawCnt <- cbind(rawCnt, temp)
    }
  }
  ### add phenotype to sample info
  sampleInfo[colnames(rawCnt)[(ncol(rawCnt)-length(f)+1):ncol(rawCnt)],"Phenotype"] <- "Gingivitis"
  
  ### periodontitis
  f <- list.files(paste0(inputDir, "periodontitis"))
  f <- f[which(endsWith(f, "counts.txt"))]
  for(i in 1:length(f)) {
    temp <- read.table(paste0(inputDir, "periodontitis/", f[i]), row.names = 1, header = TRUE)
    temp <- temp[order(rownames(temp)),,drop=FALSE]
    if(is.null(rawCnt)) {
      rawCnt <- temp
    } else {
      rawCnt <- cbind(rawCnt, temp)
    }
  }
  ### add phenotype to sample info
  sampleInfo[colnames(rawCnt)[(ncol(rawCnt)-length(f)+1):ncol(rawCnt)],"Phenotype"] <- "Periodontitis"
  
  
  ### order based on sample names
  rawCnt <- rawCnt[,order(colnames(rawCnt))]
  sampleInfo <- sampleInfo[order(rownames(sampleInfo)),]
  
  
  ### save in text
  write.table(rawCnt, file = paste0(outputDir, "raw_counts.txt"),
              sep = "\t", row.names = TRUE, col.names = TRUE)
  write.table(sampleInfo, file = paste0(outputDir, "sampleInfo.txt"),
              sep = "\t", row.names = FALSE, col.names = TRUE)
  
  
  ### make a README function for the RDA file
  README = function(){
    writeLines(paste(rep("#", 100), collapse = ""))
    writeLines("The RNA-Seq dataset is from Dr. Panos Papapanou")
    writeLines("The sequencing & alignment were done by Columbia Genome Center")
    writeLines("Coverage depth: 30M")
    writeLines("Sequencing type: single-end sequencing")
    writeLines("The number of genes: 25559")
    writeLines(paste(rep("#", 100), collapse = ""))
    writeLines("The \"rawCnt\" object has raw counts of the RNA-Seq")
    writeLines("The \"sampleInfo\" object has clinical information of the samples")
    writeLines(paste(rep("#", 100), collapse = ""))
  }
  
  
  ### save in RDA
  save(list = c("rawCnt", "sampleInfo", "README"),
       file = paste0(outputDir, "panos_rna_seq_data.rda"))
  
}




