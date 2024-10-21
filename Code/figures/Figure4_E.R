library(dplyr)
library(ggplot2)
library(plotly)


## After doing the real pvalues with glm we got some new modules for each cell type
## Then I selected the modules which has a pvalue less than 0.1 and auc of greater than 0.65
## The final pdfs are saved in panle folder of one drive and it also saved in /ix/djishnu/Javad/RA/figures/figure4B/
## The main output of this file is the report. the goal is to compare the median of the distribution matched with the actual performance.



## Filtering signficant modules for each cell type
cell_types = c("B cell","T cell","Mono","Fibro")
color_pallete <- c("pink","darkgreen","blue","orange")





report <- NULL

for(cell_name in  cell_types){

idx <- which(cell_types %in% cell_name )

color <- color_pallete[idx]

cat(paste0("processing cell type :",cell_name ,"\n"))
dir.create(paste0("/ix/djishnu/Javad/RA/figures/figure4B/",cell_name,"/"))
Pvalres <- read.csv(paste0("/ix/djishnu/Javad/RA/Data_Analysis/results/July22Results/",cell_name,"-pvals.txt"),sep=" ") ## Exact Pvalues

Pvalres <- Pvalres %>% distinct()
ii <- which((Pvalres$adjustedpval<0.1) & Pvalres$auc>0.69)





path <- paste0("/ix/djishnu/Javad/RA/Data_Analysis/results/July22Results/plots/",cell_name)
rdsList    = list.files(path=path,pattern = ".rds")
FileLists <- paste0("plot_k",ii,".rds")


plotresList <- list()

for( f in FileLists){

plotresList  <- c(plotresList,list(readRDS(paste0("/ix/djishnu/Javad/RA/Data_Analysis/results/July22Results/plots/",cell_name,"/",f)))) ## Emperical Pvalues

}


names(plotresList) <- ii



for(f in 1:length(ii)){
# Assuming data is assigned from plotresList
data <- plotresList[[f]]

## Filter empirical Pvalues


actual        <- data %>% filter(typ=="real_module") %>% dplyr::select(auc) %>% unique()
randpermut    <- data %>% filter(typ=="random_module")
empericalPval <- sum(randpermut$auc >actual$auc)/length(randpermut$auc)





# Filter the data
filtered_data <- data %>%
  filter(typ %in% c("real_module", "dist_module"))




# Ensure "real_module" is on the left
filtered_data$typ <- factor(filtered_data$typ, levels = c("real_module", "dist_module"))

# Define your custom color palette
my_palette <- c("real_module" = "black", "dist_module" = "lightblue")  # Adjust colors as needed



## create density plot

# Assume filtered_data is your dataframe
actual_auc <- filtered_data$auc[filtered_data$typ == "real_module"][1]  # Get the actual AUC value

filtered_data_dist <- filtered_data %>% filter(typ=="dist_module")


real_auc  <- filtered_data %>% filter(typ=="real_module") %>% select(auc) %>%unique()
dis_df    <- filtered_data %>% filter(typ=="dist_module") %>% select(auc)
dis_auc   <- median(dis_df$auc)


report <- data.frame(module=ii[f],real_auc=real_auc,dis_auc=dis_auc,celltype=cell_name)

report_long <- report %>%
  pivot_longer(cols = c(auc, dis_auc), names_to = "metric", values_to = "value")

# Create the bar plot with bars close together
pp <- ggplot(report_long, aes(x = metric, y = value)) +
  geom_bar(stat = "identity", fill = color, color = "black", width = 0.2, position = position_dodge(width = 0.1)) +
  labs(title = paste0(cell_name, ": module ", ii[f]), x = "Metric", y = "Value") +
  theme_minimal() +
  coord_cartesian(ylim = c(0.5, 1)) +
  theme(
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank()   # Remove minor grid lines
  )


# Save the plot as a PDF file
ggsave(filename = paste0("/ix/djishnu/Javad/RA/figures/figure4E/", cell_name, "/", ii[f], "_barplot.pdf"),
       plot = pp,
       width = 7,
       height = 7)


}


write.table(report,file="/ix/djishnu/Javad/RA/Data_Analysis/results/July22Results/Figure4median.csv",quote = F,row.names = F,col.names = T,sep=",")

}













