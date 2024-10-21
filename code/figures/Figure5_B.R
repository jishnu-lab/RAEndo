## Figure 5B
## The data for this script are generated in IFNA_beta_TNF_union_feat_perm.R
##                                           INFA_beta_feat_perm.R
##                                           TNF_beta_perm.R
################################################################################


## Reading results:
rm(list = ls())
# TNF

TNFRes  <- readRDS(paste0("/ix/djishnu/Javad/RA/Data_Analysis/results/TNFresultsFromExpRes_andPath_FeatPerm",0.9,".rds"))
TNRRand <- readRDS(paste0("/ix/djishnu/Javad/RA/figures/TNFRAndFromExpRes_andPath_FeatPerm0.9.rds"))


# IFNA
IFNARes <- readRDS(paste0("/ix/djishnu/Javad/RA/Data_Analysis/results/IFNA_beta_resultsFromExpRes_andPath_FeatPerm0.9.rds"))
IFNARand <- readRDS(paste0("/ix/djishnu/Javad/RA/figures/IFNARandFromExpRes_andPath_FeatPerm0.9.rds"))


# IFNA + TNF
IFNATNF <- readRDS(paste0("/ix/djishnu/Javad/RA/Data_Analysis/results/IFNA_beta_TNF_resultsFromExpRes_andPath_FeatPerm0.9.rds"))
IFNATNFRand <- readRDS(paste0("/ix/djishnu/Javad/RA/figures/IFNA_TNF_beta_RandPerm_0.9.rds"))



## Plotting

# TNF
cell <- "Fibro"
resList <- TNFRes
RandPer <- TNRRand
plotDf <- rbind(data.frame(perf=unlist(resList[[cell]]$BestPer),type="actual"),data.frame(perf=unlist(RandPer),type="Random"))
library(ggplot2)

source("/ix/djishnu/Javad/RA/R/create_density_plot.R")

# Assuming you have a data frame named 'plotDf'
result <- create_density_plot(plotDf, "Permutation test TNF",
                              p_value_x = 0.95, p_value_y = 1.2, p_value_size = 3)
density_plot <- result$plot
p_value <- result$p_value

print(density_plot)
print(paste("P-value:", p_value))
ggsave(density_plot,file="/ix/djishnu/Javad/RA/figures/PermutationTNF.pdf",width = 10,height = 10)

#IFNA
cell <- "Fibro"
resList <- IFNARes
RandPer <- IFNARand
plotDf <- rbind(data.frame(perf=unlist(resList[[cell]]$BestPer),type="actual"),data.frame(perf=unlist(RandPer),type="Random"))
result <- create_density_plot(plotDf, "Permutation test Interferon Alpha",
                              p_value_x = 0.95, p_value_y = 1.2, p_value_size = 3)
density_plot <- result$plot
p_value <- result$p_value

print(density_plot)
print(paste("P-value:", p_value))
ggsave(density_plot,file="/ix/djishnu/Javad/RA/figures/PermutationIFNA.pdf",width = 10,height = 10)


#IFNA +TNF

cell <- "Fibro"
resList <- IFNATNF
RandPer <- IFNATNFRand
plotDf <- rbind(data.frame(perf=unlist(resList[[cell]]$BestPer),type="actual"),data.frame(perf=unlist(RandPer),type="Random"))

result <- create_density_plot(plotDf, "Permutation test TNF and Interferon Alpha",
                              p_value_x = 0.9, p_value_y = 1.2, p_value_size = 3)
density_plot <- result$plot
p_value <- result$p_value

print(density_plot)
print(paste("P-value:", p_value))
ggsave(density_plot,"/ix/djishnu/Javad/RA/figures/PermutationIFNATNF.pdf",width = 10,height = 10)



