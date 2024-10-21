#plot heritability partition results 
library(tidyverse)
library(ggpubr)
library(ggprism)
results <- readRDS('/ix/djishnu/Javad/ForPriyamvada/H2Partition.rds')
x <- results
x$module <- as.character(x$module)
x <- rbind(x, data.frame('module' = c('1', '4'), 'Heritability' = c(0, 0), 'SD' = c(0, 0), 'Sig0.05' = c(0,0), 'Sig0.1' = c(0, 0)))
x <- x %>% mutate('Stars' = case_when(Sig0.05 == 0 & Sig0.1 == 0 ~ "",
                                      Sig0.05 == 0 & Sig0.1 == 1 ~ "**",
                                      Sig0.05 == 1 & Sig0.1 == 1 ~ "**"))
plot_stars <- data.frame('group1' = x$module, 'group2' = x$module, 'p.adj' = x$Stars, y.position = x$Heritability + x$SD)
png(paste0('/ix/djishnu/Priyamvada/Immport/final_plots', '/', 'heritability_plot.png'), width = 500, height = 400)
plot <- ggplot(x) +
  geom_bar(aes(x = module, y = Heritability), , fill = "#337357", stat = 'identity', color = 'black', linewidth = 0.50, width = 0.75) +
  scale_x_discrete(limits = c("14", "13", "12", "11", "10", "9", "8", "7", "6", "5", "4", "3", "2", "1"),
                   labels = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black", linewidth = 0.75), 
        plot.title = element_text(face = 'bold'), axis.text = element_text(size = 18, color = 'black'), 
        axis.title = element_text(face = 'bold', size = 20, color = 'black'), legend.position = 'none', strip.background = element_blank(),
        strip.text.x = element_text(size = 18, color = 'black')) +
  geom_errorbar(aes(module, ymin=Heritability, ymax=Heritability+SD), width = 0.45,
                position = position_dodge(0.9)) +
  add_pvalue(plot_stars, bracket.size = 0, fontface = 'bold', label.size = 8, tip.length = 0) +
  xlab("MODULES") + ylab("HERITABILITY") 
dev.off()
ggsave(file = 'heritability_plot_v4_one_color.pdf',
       plot = plot,
       device = 'pdf',
       path = '/ix/djishnu/Priyamvada/Immport/final_plots',
       width = 20,
       height = 10,
       units = 'in'
)
