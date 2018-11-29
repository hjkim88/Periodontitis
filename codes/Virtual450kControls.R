###
#   File name : Virtual450kControls.R
#   Author    : Hyunjin Kim
#   Date      : Nov 21, 2018
#   Email     : hk2990@cumc.columbia.edu
#   Purpose   : Combine Periodontitis samples from the old and Control samples from the new data
#
#   * Because the old 450k methylation data does not have control samples,
#     we employ controls from the new methylation data - EPIC
#     we select 450k CpGs from the EPIC and control samples to make virtual control samples
#
#   Instruction
#               1. Source("Virtual450kControls.R")
#               2. Run the function "virtual450k" - specify the input files (450k & EPIC) and output directory
#               3. The virtual RDA file (Periodontitis - Old, Control - New) will be generated in the output directory
#
#   Example
#               > source("The_directory_of_Virtual450kControls.R/Virtual450kControls.R")
#               > virtual450k(methyl450kPath="./data/450k/panos_450k_data.rda",
#                             methylEpicPath="./data/Methylation/panos_methylation_data.rda",
#                             outputDir="./data/450k/virtual_450k_with_controls/")
###

virtual450k <- function(methyl450kPath="//isilon.c2b2.columbia.edu/ifs/archive/shares/bisr/Papapanou/Sept_2018/PreprocessedData/450k/panos_450k_data.rda",
                        methylEpicPath="//isilon.c2b2.columbia.edu/ifs/archive/shares/bisr/Papapanou/Sept_2018/PreprocessedData/Methylation/panos_methylation_data.rda",
                        outputDir="//isilon.c2b2.columbia.edu/ifs/archive/shares/bisr/Papapanou/Sept_2018/PreprocessedData/450k/virtual_450k_with_controls/") {
  
  ### load datasets
  load(methyl450kPath)
  load(methylEpicPath)
  
  
  ### get common CpG sites
  common_probes <- intersect(rownames(mvalues_450k), rownames(m))
  
  
  ### only keep the commons
  mvalues_450k <- mvalues_450k[common_probes,]
  m <- m[common_probes,]
  
  
  ### combine the old and the new
  m_virtual <- cbind(mvalues_450k, m)
  sample_type <- c(rep("450k_Periodontitis", ncol(mvalues_450k)), paste0(rep("EPIC_", ncol(m)), sample_info$Phenotype))
  
  
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
  
  
  ### generate a PCA plot
  pca_plot(normalizedMat = m_virtual, grp = sample_type, title = "PCA_Virtual_M_Values",
           filePath = paste0(outputDir, "PCA_Virtual_M.png"))
  
  
  ### batch effect correction
  
  ### load library
  if(!require(sva)) {
    source("https://bioconductor.org/biocLite.R")
    biocLite("sva")
    library(sva)
  }
  
  ### correct the batch effect between 450k and EPIC
  mvalues_450k_virtual <- ComBat(dat = m_virtual,
                                 batch = c(rep("450k", ncol(mvalues_450k)), rep("EPIC", ncol(m))),
                                 mod = model.matrix(~Phenotype, data = data.frame(Phenotype=c(rep("Periodontitis", ncol(mvalues_450k)), sample_info$Phenotype))))
  
  ### generate a PCA plot
  pca_plot(normalizedMat = mvalues_450k_virtual, grp = sample_type, title = "PCA_Virtual_M_Values_BC",
           filePath = paste0(outputDir, "PCA_Virtual_M_BC.png"))
  
  
  ### make sample info for the virtual data
  mvalues_450k_virtual_sample_info <- data.frame(Phenotype=c(rep("Periodontitis", ncol(mvalues_450k)), sample_info$Phenotype),
                                                 Platform=c(rep("450k", ncol(mvalues_450k)), rep("EPIC", ncol(m))),
                                                 stringsAsFactors = FALSE, check.names = FALSE)
  rownames(mvalues_450k_virtual_sample_info) <- colnames(mvalues_450k_virtual)
  
  
  ### make a README function for the virtual data
  README <- function(){
    writeLines(paste(rep("#", 100), collapse = ""))
    writeLines("Because the old 450k methylation data does not have control samples,")
    writeLines("we employed controls from the new methylation data - EPIC.")
    writeLines("We selected the same 450k CpGs from the EPIC as 450k to make virtual 450k control samples.")
    writeLines("The \"mvalues_450k_virtual\" has methylation values of common CpG sites between")
    writeLines("450k and EPIC datasets and the batch effect between two platforms was corrected.")
    writeLines("It contains both 450k and EPIC samples.")
    writeLines("The \"mvalues_450k_virtual_sample_info\" has phenotype and platform info of the samples.")
    writeLines(paste(rep("#", 100), collapse = ""))
  }
  
  
  ### save in .rda
  save(list = c("mvalues_450k_virtual", "mvalues_450k_virtual_sample_info", "README"), file = paste0(outputDir, "panos_450k_virtual_fullset.rda"))
  
  
  ### save in txt
  write.table(data.frame(Probe_ID=rownames(mvalues_450k_virtual), mvalues_450k_virtual,
                         stringsAsFactors = FALSE, check.names = FALSE),
              file = paste0(outputDir, "panos_450k_virtual_fullset.txt"),
              row.names = FALSE)
  write.table(data.frame(Sameple_Name=rownames(mvalues_450k_sample_info), mvalues_450k_sample_info,
                         stringsAsFactors = FALSE, check.names = FALSE),
              file = paste0(outputDir, "panos_450k_virtual_fullset_sample_info.txt"),
              row.names = FALSE)
  
}
