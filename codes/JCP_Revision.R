###
#   File name : JCP_Revision.R
#   Author    : Hyunjin Kim
#   Date      : Mar 24, 2021
#   Email     : hyunjin.kim@stjude.org
#   Purpose   : 1. The present study although is done very carefully, it would be interesting
#                  to get the author’s perspective on the influence of age, gender, race, and
#                  Mendelian inheritance of the study population. These variables may greatly
#                  influence differential DNA methylation/RNA expression datasets.
#                  Discussion regarding this would benefit the readers.
#               2. A pie chart and Venn diagram showing number of genes overlapping
#                  (methylated/unmethylated/differentially expressed) among health, gingivitis,
#                  and periodontitis group would be beneficial and easy to understand
#                  in addition to the tables provided.
#               3. Similarly, a bar graph representing functional/regulatory regions
#                  (promoter, enhancer, non-coding etc.) scanned would be beneficial.
#               4. Differentially methylated/expression regions on genes can be represented in
#                  genome browser view (Integrated genomics viewer). This is important not only for
#                  visualization purpose but also represents functional elements/enrichment
#                  (For example: Hypermethylated region vs. RNA expression level overlap)
#
#   Instruction
#               1. Source("JCP_Revision.R")
#               2. Run the function "jcr_revision" - specify the input files (Differential) of methylation and RNA-Seq and output directory
#               3. The revision results will be generated in the output directory
#
#   Example
#               > source("The_directory_of_JCP_Revision.R/JCP_Revision.R")
#               > jcr_revision(methylPath="./results/DMA/DMP_Periodontitis_Healthy.txt",
#                                  methylRegPath="./results/DMA/DMR_Periodontitis_Healthy.txt",
#                                  rnaseqPath="./results/DEA/DE_Results_Periodontitis_vs_Healthy.xlsx",
#                                  outputDir="./results/Revision/")
###

jcr_revision <- function(methylPath="./results/DMA/DMP_Periodontitis_Healthy.txt",
                         methylRegPath="./results/DMA/DMR_Periodontitis_Healthy.txt",
                         rnaseqPath="./results/DEA/DE_Results_Periodontitis_vs_Healthy.xlsx",
                         outputDir="./results/Revision/") {
  
  
  
  
  
}
