###
#   File name : FDRSignificanceQuartile.R
#   Author    : Hyunjin Kim
#   Date      : Oct 26, 2018
#   Email     : hk2990@cumc.columbia.edu
#   Purpose   : -log10(FDR) plot with DMPs and DEGs
#
#   Instruction
#               1. Source("FDRSignificanceQuartile.R")
#               2. Run the function "significancePlot" - specify the input files (Differential) of methylation and RNA-Seq and output directory
#               3. -log10(FDR) plot will be generated in the output directory
#
#   Example
#               > source("The_directory_of_FDRSignificanceQuartile.R/FDRSignificanceQuartile.R")
#               > significancePlot(methylPath="./results/DMA/DMP_Periodontitis_Healthy.txt",
#                                  methylRegPath="./results/DMA/DMR_Periodontitis_Healthy.txt",
#                                  rnaseqPath="./results/DEA/DE_Results_Periodontitis_vs_Healthy.xlsx",
#                                  geneRIFPath="./data/generifs_basic.txt",
#                                  fdrThreshold=0.05,
#                                  outputDir="./results/Combined/")
###

significancePlot <- function(methylPath="//isilon.c2b2.columbia.edu/ifs/archive/shares/bisr/Papapanou/Sept_2018/DifferentialMethylation/DMP_Periodontitis_Healthy.txt",
                             methylRegPath="//isilon.c2b2.columbia.edu/ifs/archive/shares/bisr/Papapanou/Sept_2018/DifferentialMethylation/DMR_Periodontitis_Healthy.txt",
                             rnaseqPath="//isilon.c2b2.columbia.edu/ifs/archive/shares/bisr/Papapanou/Sept_2018/DifferentialExpression/DE_Results_Periodontitis_vs_Healthy.xlsx",
                             geneRIFPath="//isilon.c2b2.columbia.edu/ifs/archive/shares/bisr/Papapanou/Sept_2018/PreprocessedData/generifs_basic.txt",
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
  if(!require(ggrepel)) {
    source("https://bioconductor.org/biocLite.R")
    biocLite("ggrepel")
    library(ggrepel)
  }
  
  
  ### load data
  dmp <- read.table(file = methylPath, header = TRUE, sep = "\t", check.names = FALSE)
  dmr <- read.table(file = methylRegPath, header = TRUE, sep = "\t", check.names = FALSE)
  de <- read.xlsx2(file = rnaseqPath, sheetIndex = 1, stringsAsFactors = FALSE)
  
  
  ### set row names
  rownames(dmp) <- dmp$Name
  rownames(de) <- de$Gene_Symbol
  
  
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
        cor_data <- rbind(cor_data, c(temp[j], paste0("DMR", i), as.numeric(c(de[temp[j],"padj"], dmr$Stouffer[i]))))
      }
    }
  }
  cor_data <- data.frame(cor_data)
  cor_data[,1:2] <- apply(cor_data[,1:2], 2, as.character)
  cor_data[,3:4] <- apply(cor_data[,3:4], 2, as.numeric)
  
  
  ### add CpG info
  cor_data$X5 <- NA
  for(i in 1:nrow(cor_data)) {
    ### get DMR index
    idx <- as.numeric(substring(cor_data[i,2], 4))
    
    ### narrow down with the specific chromosome
    temp <- dmp[which(as.character(dmp$chr) == as.character(dmr$seqnames[idx])), 1:4]
    
    ### get CpG names
    temp <- as.character(temp$Name[intersect(which(temp$pos >= dmr$start[idx]), which(temp$pos <= dmr$end[idx]))])
    
    ### add the info
    cor_data$X5[i] <- paste(temp, collapse = ",")
  }
  
  
  ### write out the text file
  colnames(cor_data) <- c("Gene_Name", "DMR_Name", "FDR_Expression", "FDR_Methylation", "CpG_Names")
  cor_data <- cor_data[order(cor_data[,3]),]
  write.xlsx2(cor_data, file = paste0(outputDir, "List_-log10(FDR)_", substr(basename(methylRegPath), 1, nchar(basename(methylRegPath))-4), "_all.xlsx"), row.names = FALSE)
  
  
  ### mark the significant genes
  cor_data$Label <- ""
  idx <- intersect(which(cor_data$FDR_Expression < fdrThreshold), which(cor_data$FDR_Methylation < fdrThreshold))
  cor_data$Label[idx] <- cor_data$Gene_Name[idx]
  cor_data[,5:6] <- apply(cor_data[,5:6], 2, as.character)
  
  
  ### make a correlation plot between DM regions and DE genes
  df <- data.frame(-log10(cor_data[,3:4]), Label=cor_data$Label)
  df[,3] <- sapply(df[,3], as.character)
  fName <- paste0("Plot_-log10(FDR)_", substr(basename(methylRegPath), 1, nchar(basename(methylRegPath))-4), "_all.png")
  ggplot(data = df, aes(x=FDR_Expression, y=FDR_Methylation)) +
    geom_point(color = "black", size = 1) +
    geom_label_repel(aes(FDR_Expression, FDR_Methylation, label = Label), color = "red", box.padding = unit(0.45, "lines")) +
    labs(title=substr(fName, 1, nchar(fName)-4),
         subtitle=sprintf("P.Cor = %s, p-value = %s",
                          round(cor(df[,1], df[,2], use = "pairwise.complete.obs"), 5),
                          signif(cor.test(df[,1], df[,2])$p.value, 5))) +
    xlab("Expression") +
    ylab("Methylation") +
    # geom_smooth(method = lm, color="gray", se=FALSE) +
    geom_vline(aes(xintercept = -log10(fdrThreshold), linetype = paste0("-log10(", fdrThreshold, ")")), color = "red") +
    geom_hline(aes(yintercept = -log10(fdrThreshold), linetype = paste0("-log10(", fdrThreshold, ")")), color = "red") +
    scale_linetype_manual(name = "FDR Threshold", values = c(2),
                          guide = guide_legend(override.aes = list(color = c("red")))) +
    theme_classic(base_size = 16)
  ggsave(filename = paste0(outputDir, fName), width = 10, height = 10)
  
  
  ### create GeneRIF text
  if(length(idx) > 0) {
    ### load geneRIF
    geneRIF <- read.table(file = geneRIFPath, header = FALSE, sep = "\t", check.names = FALSE)
    
    # *****************************************************************************
    #
    # Map Gene Symbols to Gene IDs
    #
    # geneSymbol:	A single string or a vector of strings representing gene symbol(s)
    # 
    # Returns a vector of the same size as "geneSymbol" where the i-th entry is the 
    # gene ID (as an integer) corresponsing to the i-th gene symbol in "geneSymbol". 
    # The return vector entries are named using the gene symbols. For gene symbols mapped
    # to more than one entrez ids, only the first id is returned. 
    # ATTENTION: 
    # Before used geneSymbol is first stripped of symbols which are not mapped to 
    # at least one  entrez ID. If no gene symbols remain after this stripping, a NULL 
    # value is returned.
    # *****************************************************************************
    geneSymbolToEntrezId <- function(geneSymbol){
      
      if(!require(org.Hs.eg.db)) {
        source("https://bioconductor.org/biocLite.R")
        biocLite("org.Hs.eg.db")
        library(org.Hs.eg.db)
      }
      
      keys = keys(org.Hs.eg.db,  keytype="SYMBOL")
      geneSymbol = geneSymbol[geneSymbol %in% keys]
      if (length(geneSymbol) == 0)
        return(NULL)
      
      tmp = select(org.Hs.eg.db, keys=geneSymbol, keytype="SYMBOL", columns=c("ENTREZID"))
      # Remove duplicate mappings, if any
      t = which(duplicated(tmp[,1]))
      if (length(t) > 0)
        tmp = tmp[-t, ]
      res = as.integer(tmp[,2])
      names(res) = tmp[, 1]
      return(res)
    }
    
    ### get Entrez IDs
    entIDs <- geneSymbolToEntrezId(unique(cor_data$Gene_Name[idx]))
    
    ### extract geneRIF for the specific genes
    add_info <- geneRIF[which(geneRIF$V2 %in% entIDs),c(2,5)]
    colnames(add_info) <- c("Gene_Name", "Gene_RIF")
    
    ### change Entrez IDs to gene symbols
    name_map <- names(entIDs)
    names(name_map) <- entIDs
    add_info[,1] <- name_map[as.character(add_info[,1])]
    
    ### remove duplicates
    rIdx <- duplicated.data.frame(add_info)
    add_info <- add_info[!rIdx,]
    
    ### write out the result
    write.xlsx2(add_info, file = paste0(outputDir, "GeneRIF_", substr(basename(methylRegPath), 1, nchar(basename(methylRegPath))-4), "_all.xlsx"), row.names = FALSE)
  }
  
}
