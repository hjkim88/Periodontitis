###
#   File name : PreprocessOlds.R
#   Author    : Hyunjin Kim
#   Date      : Nov 19, 2018
#   Email     : hk2990@cumc.columbia.edu
#   Purpose   : Collect Affy and 450k methylation data and transfrom them into RDA files suitable for the combined analysis
#
#   Instruction
#               1. Source("PreprocessOlds.R")
#               2. Run the function "preprocessOlds" - specify the input files (Affy and 450k) and output file paths
#               3. Two RDA files (one for Affy and the other for 450k) will be generated in the output file paths
#
#   Example
#               > source("The_directory_of_PreprocessOlds.R/PreprocessOlds.R")
#               > preprocessOlds(affyInputFilePath="./data/Affy/PeriDataComBatAdjusted_rf2.RData",
#                                methyl450kInputMFilePath="./data/450k/Mvalues.quantile.RData",
#                                methylSampleInfoPath="./data/450k/sampleSheet.csv",
#                                affyOutputFilePath="./data/Affy/panos_affy_data.rda",
#                                methyl450kOutputFilePath="./data/450k/panos_450k_data.rda")
###

preprocessOlds <- function(affyInputFilePath="//isilon.c2b2.columbia.edu/ifs/archive/shares/bisr/Papapanou/Sept_2018/OriginalData/affy/PeriDataComBatAdjusted_rf2.RData",
                           methyl450kInputMFilePath="//isilon.c2b2.columbia.edu/ifs/archive/shares/bisr/Papapanou/Sept_2018/OriginalData/450k/Mvalues.quantile.RData",
                           methylSampleInfoPath="//isilon.c2b2.columbia.edu/ifs/archive/shares/bisr/Papapanou/Sept_2018/OriginalData/450k/sampleSheet.csv",
                           affyOutputFilePath="//isilon.c2b2.columbia.edu/ifs/archive/shares/bisr/Papapanou/Sept_2018/PreprocessedData/Affy/panos_affy_data.rda",
                           methyl450kOutputFilePath="//isilon.c2b2.columbia.edu/ifs/archive/shares/bisr/Papapanou/Sept_2018/PreprocessedData/450k/panos_450k_data.rda") {
  
  ### load library
  if(!require(hgu133plus2.db)) {
    source("https://bioconductor.org/biocLite.R")
    biocLite("hgu133plus2.db")
    library(hgu133plus2.db)
  }
  
  
  ### load datasets
  load(affyInputFilePath)
  load(methyl450kInputMFilePath)
  mvalues_450k_sample_info <- read.table(file = methylSampleInfoPath, header = TRUE, sep = ",",
                                         check.names = FALSE, stringsAsFactors = FALSE)
  
  
  ### reorganize and rename the structures
  
  ### mvalues_450k_sample_info
  rownames(mvalues_450k_sample_info) <- mvalues_450k_sample_info$Sample_Name
  common_450k_sample_name <- intersect(colnames(Mvalues.quantile), rownames(mvalues_450k_sample_info))
  mvalues_450k_sample_info <- mvalues_450k_sample_info[common_450k_sample_name,]
  Mvalues.quantile <- Mvalues.quantile[,common_450k_sample_name]
  mvalues_450k_sample_info <- mvalues_450k_sample_info[,-c(2,3)]
  
  ### dsn & cmb.dat
  rownames(dsn) <- dsn$mRNA
  common_affy_sample_name <- intersect(colnames(cmb.dat), rownames(dsn))
  dsn <- dsn[common_affy_sample_name,]
  cmb.dat <- cmb.dat[,common_affy_sample_name]
  rownames(dsn) <- dsn$Biopsy
  dsn <- dsn[,-2]
  colnames(cmb.dat) <- rownames(dsn)
  
  ### order
  cmb.dat <- as.data.frame(cmb.dat[,order(as.numeric(colnames(cmb.dat)))])
  dsn <- dsn[colnames(cmb.dat),]
  Mvalues.quantile <- as.data.frame(Mvalues.quantile[,order(as.numeric(colnames(Mvalues.quantile)))])
  mvalues_450k_sample_info <- mvalues_450k_sample_info[colnames(Mvalues.quantile),]
  
  
  # ### keep the same samples between Affy and 450k
  # common_sample <- intersect(rownames(dsn), rownames(mvalues_450k_sample_info))
  # affy_norm_ge_perio <- as.data.frame(cmb.dat[,common_sample])
  # affy_norm_ge_perio_sample_info <- dsn[common_sample,]
  # mvalues_450k_perio <- as.data.frame(Mvalues.quantile[,common_sample])
  # mvalues_450k_perio_sample_info <- mvalues_450k_sample_info[common_sample,]
  
  
  ### Affy
  
  ### rename
  affy_norm_ge <- cmb.dat
  affy_norm_ge_sample_info <- dsn
  
  ### convert microarray probes to gene symbols
  ### hgu133plus2
  ### get gene symbol mapping info
  x <- hgu133plus2SYMBOL
  mapped_probes <- mappedkeys(x)
  xx <- as.list(x[mapped_probes])
  
  ### get common probes
  common_probes <- intersect(rownames(affy_norm_ge), names(xx))
  
  ### keep common probes only
  affy_norm_ge <- affy_norm_ge[common_probes,]
  
  ### annotate gene symbols
  affy_norm_ge <- data.frame(Gene_Symbol=sapply(xx[common_probes], function(y) {
                                                                     return(y[1])
                                                                   }),
                             affy_norm_ge, stringsAsFactors = FALSE, check.names = FALSE)
  
  ### make a README function for Affy RDA
  README = function(){
    writeLines(paste(rep("#", 100), collapse = ""))
    writeLines("Dr. Papapanou's Affy chip data and the sample info were collected and combined into one RDA file.")
    writeLines("The \"affy_norm_ge\" object has normalized gene expression levels of the Affy chip.")
    writeLines("The \"affy_norm_ge_sample_info\" object has sample information of the Affy chip.")
    writeLines(paste(rep("#", 100), collapse = ""))
  }
  
  ### save in .rda
  save(list = c("affy_norm_ge", "affy_norm_ge_sample_info", "README"), file = affyOutputFilePath)
  
  
  ### 450k
  
  ### rename
  mvalues_450k <- Mvalues.quantile
  mvalues_450k_sample_info <- mvalues_450k_sample_info
  
  ### make a README function for 450k RDA
  README = function(){
    writeLines(paste(rep("#", 100), collapse = ""))
    writeLines("Dr. Papapanou's 450k methylation data and the sample info were collected and combined into one RDA file.")
    writeLines("The \"mvalues_450k\" object has methylation levels based on M values of the 450k chip.")
    writeLines("The \"mvalues_450k_sample_info\" object has sample information of the 450k chip.")
    writeLines(paste(rep("#", 100), collapse = ""))
  }
  
  ### save in .rda
  save(list = c("mvalues_450k", "mvalues_450k_sample_info", "README"), file = methyl450kOutputFilePath)
  
}
