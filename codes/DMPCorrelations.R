###
#   File name : DMPCorrelations.R
#   Author    : Hyunjin Kim
#   Date      : Oct 17, 2018
#   Email     : hk2990@cumc.columbia.edu
#   Purpose   : Calculate correlations between before and after the kit effect correction
#
#   Instruction
#               1. Source("DMPCorrelations.R")
#               2. Run the function "getCor" - specify the input files (differential methylation), an attribute you want to compare, and output file path
#               3. A correlation plot between two results will be generated in the output path 
#
#   Example
#               > source("The_directory_of_DMPCorrelations.R/DMPCorrelations.R")
#               > getCor(beforeFilePath="./results/DMA/DMP_Periodontitis_Healthy.txt",
#                        afterFilePath="./results/DMA/DMP_Periodontitis_Healthy_kc.txt",
#                        compare=c("PValue", "logFC"),
#                        outputFilePath="./results/Cor/Cor_Periodontitis_Healthy.png")
###

getCor <- function(beforeFilePath="//isilon.c2b2.columbia.edu/ifs/archive/shares/bisr/Papapanou/Sept_2018/DifferentialMethylation/DMP_Periodontitis_Healthy.txt",
                   afterFilePath="//isilon.c2b2.columbia.edu/ifs/archive/shares/bisr/Papapanou/Sept_2018/DifferentialMethylation/DMP_Periodontitis_Healthy_kc.txt",
                   compare=c("padj", "logFC"),
                   outputFilePath="//isilon.c2b2.columbia.edu/ifs/archive/shares/bisr/Papapanou/Sept_2018/DifferentialMethylation/Cor_Periodontitis_Healthy.png") {
  
  ### load library
  if(!require(ggrepel)) {
    source("https://bioconductor.org/biocLite.R")
    biocLite("ggrepel")
    library(ggrepel)
  }
  
  
  ### load datasets
  bd <- read.table(file = beforeFilePath, header = TRUE, sep = "\t", check.names = FALSE)
  ad <- read.table(file = afterFilePath, header = TRUE, sep = "\t", check.names = FALSE)
  
  
  ### A function to make combined data, to calculate correlations, and to print a plot
  calCor <- function(beforeData, afterData, compare, outputPath, fName) {
    ### get FDR or logFC values
    if(compare == "padj") {
      x <- beforeData$adj.P.Val
      y <- afterData$adj.P.Val
      
      x_lab <- "Before_KC_FDR"
      y_lab <- "After_KC_FDR"
    } else if(compare == "logFC") {
      x <- beforeData$logFC
      y <- afterData$logFC
      
      x_lab <- "Before_KC_logFC"
      y_lab <- "After_KC_logFC"
    }
    
    ### name the vectors
    names(x) <- beforeData$Name
    names(y) <- afterData$Name
    
    ### remove NAs
    if(length(which(is.na(x))) > 0) {
      x <- x[-which(is.na(x))]
    }
    if(length(which(is.na(y))) > 0) {
      y <- y[-which(is.na(y))]
    }
    
    ### only extract commonly shared genes
    common <- intersect(names(x), names(y))
    x <- x[common]
    y <- y[common]
    
    ### make a one data frame
    df <- data.frame(x,y)
    
    ### change factor columns to character columns
    idx <- sapply(df, is.factor)
    df[idx] <- lapply(df[idx], function(x) as.character(x))
    
    ### plot the correlation with ggplot
    ggplot(data = df, aes(x=x, y=y)) +
      geom_point(color = "black", size = 1) +
      labs(title=substr(fName, 1, nchar(fName)-4), subtitle=sprintf("P.Cor = %s, p-value = %s", round(cor(df$x, df$y), 5), signif(cor.test(df$x, df$y)$p.value, 5)), x=x_lab, y=y_lab) +
      geom_smooth(method = lm, color="gray", se=FALSE) +
      theme_classic(base_size = 16)
    ggsave(filename = paste0(outputPath, fName), width = 10, height = 10)
  }
  
  
  ### make a correlation plot
  calCor(beforeData = bd, afterData = ad, compare = compare[1],
         outputPath = paste0(dirname(outputFilePath), "/"), fName = basename(outputFilePath))
  
}
