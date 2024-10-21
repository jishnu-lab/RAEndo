##Process hotnet output##
##Code updated based on discussion with Jishnu##
##Need to process all the delta outputs and take the smallest module size >3 as the significant module size##
library(tidyverse)
library(dplyr)
args = commandArgs(trailingOnly=TRUE)
hotnet_output_directory <- "/ix/djishnu/Priyamvada/Immport/Network_analysis/Network_modules_updated/gene_scores_ldak_500_permute/homosapiens_binary_co_complex_feb2023-ldak_score_ra_prg65" ##directory where the output from hotnet is stored 
ppi_list <- read.table("/ix/djishnu/Priyamvada/Immport/Network_analysis/Networks/HomoSapiens_binary_co_complex_union.txt")  ##table containing information about protein-protein interactions
gene_score_directory <- '/ix/djishnu/Priyamvada/Immport/forJavad' ##Gene scores that were propagated during hotnet
gene_score_file_name <- 'LDAK_score_RA_prg65.csv' 
gene_scores <- read.table(paste0(gene_score_directory, '/', gene_score_file_name))
gene_scores <- rename(gene_scores, Gene = V1, Score = V2)
delta_directories <- list.dirs(hotnet_output_directory, recursive = F)
delta_values <- str_remove(delta_directories, "/ix/djishnu/Priyamvada/Immport/Network_analysis/Network_modules_updated/gene_scores_ldak_500_permute/homosapiens_binary_co_complex_feb2023-ldak_score_ra_prg65/delta_")
threshold = c()
size <- c()
gene_names <- data.frame("Source" = character(), "Target" = character()) #empty data fame
gene_scores_significant_modules <- data.frame("Gene" = character(), "Score" = numeric())
for (dir in 1:length(delta_directories)){
significant_results <- read.table(paste0(delta_directories[dir], "/", "significance.txt"), header = T)
#no of significant results 
  if (length(which(significant_results$p.value <= 0.1)) == 0) {
    print(paste0("No significant results for", " ", delta_directories[dir]))
  } else if (length(which(significant_results$p.value <= 0.1)) == 1) {
    print(paste0("1 significant module for", " ", delta_directories[dir]))
    significant_modules <- significant_results %>% filter(p.value <= 0.1 & Size > 3)
    size <- significant_modules$Size
    threshold <- c(threshold, paste0(size, "for", delta_values[dir]))
    modules <- read.table(paste0(delta_directories[dir], "/", "components.txt"), header = F, fill = TRUE)
    modules <- modules %>% slice(1:size)
    gene_names_list <- as.vector(t(modules))
    gene_names_list <- unique(gene_names_list)
    gene_names_ppi <- ppi_list[ppi_list$V1 %in% gene_names_list,]
    gene_names_ppi <- gene_names_ppi[gene_names_ppi$V2 %in% gene_names_list,] #filter out only genes present in ppi ref
    gene_names_ppi <- rename(gene_names_ppi, "Source" = "V1", "Target" = "V2") 
    gene_scores_ppi <- data.frame("Gene" = unique(c(gene_names_ppi$Source, gene_names_ppi$Target)))
    gene_scores_ppi <- merge(gene_scores_ppi, gene_scores, by.x = "Gene", by.y = "Gene") #get gene scores for both target and source genes 
    colnames(gene_scores_ppi) <- c("Gene", "Score")
    gene_names <- rbind(gene_names, gene_names_ppi)
    gene_scores_significant_modules <- rbind(gene_scores_significant_modules, gene_scores_ppi)
  } else { 
    print(paste0("Multiple significant modules for", " ", delta_directories[dir]))
    significant_modules <- significant_results %>% filter(p.value <= 0.1 & Size > 3)
    size <- min(significant_modules$Size)
    threshold <- c(threshold, paste0(size, "for", delta_values[dir]))
    modules <- read.table(paste0(delta_directories[dir], "/", "components.txt"), header = F, fill = TRUE)
    modules <- modules %>% slice(1:size)
    gene_names_list <- as.vector(t(modules))
    gene_names_list <- unique(gene_names_list)
    gene_names_ppi <- ppi_list[ppi_list$V1 %in% gene_names_list,]
    gene_names_ppi <- gene_names_ppi[gene_names_ppi$V2 %in% gene_names_list,] #filter out only genes present in ppi ref
    gene_names_ppi <- rename(gene_names_ppi, "Source" = "V1", "Target" = "V2") 
    gene_scores_ppi <- data.frame("Gene" = unique(c(gene_names_ppi$Source, gene_names_ppi$Target)))
    gene_scores_ppi <- merge(gene_scores_ppi, gene_scores, by.x = "Gene", by.y = "Gene") #get gene scores for both target and source genes 
    colnames(gene_scores_ppi) <- c("Gene", "Score")
    gene_names <- rbind(gene_names, gene_names_ppi)
    gene_scores_significant_modules <- rbind(gene_scores_significant_modules, gene_scores_ppi)
          }
}
##remove duplicated rows
print(threshold)
gene_names <- gene_names[!duplicated(gene_names),]
gene_scores_significant_modules <- gene_scores_significant_modules[!duplicated(gene_scores_significant_modules),]
#save gene name and scores in separate file 
#gene_score_and_names <- gene_names[,c("Source", "Source_gene_score")]
##save file 
write.table(gene_names, file = paste0(hotnet_output_directory, "/", "significant_modules.txt"), sep = "\t", row.names = FALSE, quote = FALSE)
write.table(gene_scores_significant_modules, file = paste0(hotnet_output_directory, "/", "gene_score_for_significant_genes.txt"), sep = "\t", row.names = FALSE, quote = FALSE)
