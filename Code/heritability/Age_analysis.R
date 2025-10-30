library(dplyr)
library(ggplot2)
library(readr)

data <- readxl::read_xlsx("/Users/javad/Library/CloudStorage/OneDrive-UniversityofPittsburgh/RA/megacp_998_shared_051418.xlsx",sheet = "megacp_998_051418",skip=1)
data_Rf_Positive <- data %>% filter(`RF +/-`==1)
data_Rf_Positive$Pheno <- as.numeric(data_Rf_Positive$`CCP +/-`) * as.numeric(data_Rf_Positive$`RF +/-`)





df <- data_Rf_Positive %>%
  mutate(
    AGEvis1 = parse_number(AGEvis1),
    Pheno = factor(Pheno, levels = c(0, 1), labels = c("Pheno 0", "Pheno 1"))
  )

stats <- df %>%
  group_by(Pheno) %>%
  summarise(
    mean_age = mean(AGEvis1, na.rm = TRUE),
    sd_age   = sd(AGEvis1, na.rm = TRUE),
    n        = sum(!is.na(AGEvis1)),
    se       = sd_age / sqrt(n),
    ci95     = 1.96 * se,
    .groups = "drop"
  )

ggplot(stats, aes(x = Pheno, y = mean_age)) +
  geom_boxplot() +
  geom_errorbar(aes(ymin = mean_age - ci95, ymax = mean_age + ci95), width = 0.15) +
  labs(x = "Phenotype", y = "Mean age (years)", title = "Mean age by phenotype (±95% CI)") +
  theme_minimal()


library(ggplot2)

# keep only rows with both variables present
d2 <- droplevels(df[!is.na(df$AGEvis1) & !is.na(df$Pheno), c("AGEvis1","Pheno")])

# nonparametric p-value (use t.test(...) if you prefer)
pval  <- wilcox.test(AGEvis1 ~ Pheno, data = d2)$p.value
y_max <- max(d2$AGEvis1, na.rm = TRUE)

d2$Pheno <- factor(d2$Pheno,
                   levels = c("Pheno 0","Pheno 1"),
                   labels = c("CCP-RF+","CCP+RF+"))


AgePlot <-ggplot(d2, aes(x = Pheno, y = AGEvis1)) +
  geom_boxplot(
    width = 0.5,
    outlier.shape = NA,     # drop plotted outliers; remove this line to show them
    fill = "forestgreen",
    color = "forestgreen",
    alpha = 0.6
  ) +
  labs(x = "Pheno", y = "Age", title = "") +
  annotate("text", x = 1.5, y = y_max * 1.03,
           label = sprintf("Wilcoxon p = %.3g", pval), fontface = "bold") +
  expand_limits(y = y_max * 1.06) +
  theme_classic(base_size = 12) +
  theme(legend.position = "none")

ggsave(AgePlot,file="/Users/javad/Library/CloudStorage/OneDrive-UniversityofPittsburgh/RA/Figures/Final_presubmission_draft/Revision/Revision_Panels/Age.pdf")


