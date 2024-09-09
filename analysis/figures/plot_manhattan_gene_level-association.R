library(qqman)
library(ggrepel)
library(tidyverse)
library(data.table)
library(stringr)
data <- read.table('/ix/djishnu/Javad/ForPriyamvada/LDAK_score.csv')
colnames(data) <- c('Gene', 'log_10_p_val')
gene_pos <- fread('/ix/djishnu/Priyamvada/Covid_Flu/RNA_seq_analysis/ref_files_bedtools/Gene_length.bed')
gene_pos <- gene_pos[, -5]
colnames(gene_pos) <- c('CHR', 'Start', 'End', 'Gene')
gene_pos$CHR <- as.integer(str_remove(gene_pos$CHR, 'chr'))
data <- merge(data, gene_pos)
plot_data <- data[, c('Gene', 'CHR', 'Start', 'End', 'log_10_p_val')]
adjusted_p_val <- -log10(0.05/length(data$Gene))
snp_to_annotate <- data[data$log_10_p_val >= adjusted_p_val, 'Gene']
snps_of_interset <- data[data$log_10_p_val >= adjusted_p_val, 'Gene']
plot_df <- plot_data %>% 
  
  # Compute chromosome size
  group_by(CHR) %>% 
  summarise(chr_len=max(End)) %>% 
  
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(as.numeric(chr_len))-chr_len) %>%
  select(-chr_len) %>%
  
  # Add this info to the initial dataset
  left_join(plot_data, ., by=c("CHR"="CHR")) %>%
  
  # Add a cumulative position of each SNP
  arrange(CHR, End) %>%
  mutate( BPcum=End+tot)
plot_df <- plot_df %>% mutate(is_highlight=ifelse(log_10_p_val >= adjusted_p_val, "yes", "no"))
plot_df <- plot_df %>% mutate(is_annotate=ifelse(log_10_p_val >= adjusted_p_val, "yes", "no"))
axisdf = plot_df %>% group_by(CHR) %>% summarize(center=(max(BPcum) + min(BPcum) ) / 2 )
manhattan_plot <- ggplot(plot_df, aes(x=BPcum, y=log_10_p_val)) +
  
  # Show all points
  geom_point( aes(color=as.factor(CHR)), alpha=0.8, size=1.3) +
  #scale_color_manual(values = rep(c("#C1C0C8", "#CCCCFF"), 22 )) +
  
  # custom X axis:
  scale_x_continuous( label = axisdf$CHR, breaks= axisdf$center ) +
  scale_y_continuous(expand = c(0, 0), breaks = seq(0, 16, 2), limits = c(0, 16)) +  # remove space between plot area and x axis
  #add the geomlines 
  geom_hline(yintercept = adjusted_p_val, linetype = "dashed", color = "#a70000", size = 0.5) +
  #geom_hline(yintercept = 5, linetype = "dashed", color = "#ff5252", size = 0.5) +
  #add label 
  geom_point(data=subset(plot_df, is_annotate=="yes"), size=1) +
  geom_label_repel(data=subset(plot_df, is_annotate=="yes"), aes(label=Gene), size=4, max.overlaps = Inf) +
  
  # To this line
  geom_label_repel(data = subset(plot_df, is_annotate == "yes"), aes(label = Gene), size = 4, max.overlaps = Inf) +
  scale_color_manual(values = rep(c("#11132f", "#7681ce"), 22 ), guide = FALSE) +
  
  # Custom the theme:
  theme_bw() +
  theme( 
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(), 
    axis.title = element_text(face = 'bold', size = 10, color = 'black')
  ) + 
  xlab("CHROMOSOMES") + ylab("NEGATIVE LOG OF P-VALUE")  
ggsave('/ix/djishnu/Priyamvada/Immport/final_plots/manhattan_plot_with_gene_LDAK.pdf',
       manhattan_plot,
       device = 'pdf', width = 12,
       height = 5,
       units = 'in')
#plot only snps with p-val above 
p.adjust(0.05, method = 'BH', n = 14)

