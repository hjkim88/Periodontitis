###
#   File name : Split_Matrix_Data.R
#   Author    : Hyunjin Kim
#   Date      : Apr 20, 2021
#   Email     : firadazer@gmail.com
#   Purpose   : split RNA-Seq and methylation data by each sample for GEO submission
#
#   Instruction
#               1. Source("Split_Matrix_Data.R")
#               2. Run the function "split_matrix" - specify the input file paths and the output directory
#               3. The results will be generated under the output directory
#
#   Example
#               > source("The_directory_of_Split_Matrix_Data.R/Split_Matrix_Data.R")
#               > split_matrix(matrix1_path="./etc/GEO_Submission/rnaseq_raw_counts.txt",
#                              matrix2_path="./etc/GEO_Submission/methylation_m_values.txt",
#                              output_dir="./etc/GEO_Submission/")
###

split_matrix <- function(matrix1_path="./etc/GEO_Submission/rnaseq_raw_counts.txt",
                         matrix2_path="./etc/GEO_Submission/methylation_m_values.txt",
                         output_dir="./etc/GEO_Submission/") {
  
  ### read rna-seq raw counts
  mat1 <- read.table(matrix1_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
  
  ### new output dir for rnaseq
  output_dir2 <- paste0(output_dir, "rnaseq_raw_counts/")
  dir.create(output_dir2)
  
  ### write out the splited data
  for(col in colnames(mat1)) {
    result_mat <- data.frame(Gene=rownames(mat1),
                             mat1[,col],
                             stringsAsFactors = FALSE, check.names = FALSE)
    colnames(result_mat) <- c("Gene", col)
    write.table(result_mat,
                file = paste0(output_dir2, "rnaseq_raw_counts_", col, ".txt"),
                sep = "\t", quote = FALSE, row.names = FALSE)
  }
  
  ### read methylation levels
  mat2 <- read.table(matrix2_path, header = TRUE, sep = "\t", row.names = 1,
                     stringsAsFactors = FALSE, check.names = FALSE)
  
  ### new output dir for methylation
  output_dir2 <- paste0(output_dir, "methylation_levels/")
  dir.create(output_dir2)
  
  ### write out the splited data
  for(col in colnames(mat2)) {
    result_mat <- data.frame(Gene=rownames(mat2),
                             mat2[,col],
                             stringsAsFactors = FALSE, check.names = FALSE)
    colnames(result_mat) <- c("Gene", col)
    write.table(result_mat,
                file = paste0(output_dir2, "methylation_levels_", col, ".txt"),
                sep = "\t", row.names = FALSE)
  }
  
}
