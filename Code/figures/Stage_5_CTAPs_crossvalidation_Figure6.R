library(ggplot2)
library(dplyr)
setwd("/ix/djishnu/Javad/RA/figures/figure6/")

auc_results_tr <- readRDS("/ix/djishnu/Javad/RA/auc_results_treatment.rds")

# Filter data to include only AUC values greater than 0.65
filtered_data <- auc_results_tr %>%
  filter(auc > 0.65) %>%
  mutate(module = factor(module))  # Ensure module is a factor for proper x-axis treatment

# Create the hierarchical bar plot without grid lines

p_tr<- ggplot(filtered_data, aes(x = interaction(module, celltype), y = auc, fill = treatment)) +
  geom_bar(stat = "identity", position = position_dodge2(preserve = "single", padding = 0), width = 0.7) +
  scale_x_discrete(labels = function(x) {
    # Add a line break after the module number for a hierarchical appearance
    sapply(strsplit(x, "\\."), function(x) paste(x, collapse = "\n"))
  }) +
  labs(
    title = "AUC Values > 0.65 by Cell Type, Module, and Treatment",
    x = "Module and Cell Type",
    y = "AUC",
    fill = "Treatment"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.spacing = unit(1, "lines"),
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines

    # Add axis lines
    axis.line = element_line(color = "black", linewidth = 0.5, linetype = "solid")
  )+ scale_fill_manual(
    name = "Treatment",  # Title of the legend
    values = c("DMARD.naive" = "blue", "MTX.inadequate" = "red", "TNFi.inadequate" = "green"),  # Example colors
    labels = c(
      "DMARD.naive" = "No treatment",
      "MTX.inadequate" = "Methotrexate",
      "TNFi.inadequate" = "TNFi inadequate"
    )
  )+coord_cartesian(ylim = c(0.5, 0.9))

ggsave("figure6_treatment.pdf", p_tr, width = 8, height = 6, dpi = 300)

library(ggplot2)
library(dplyr)

# Apply Cartesian filter (ylim) to focus on AUC > 0.5
pctap <- ggplot(filtered_data, aes(x = interaction(module, celltype), y = auc, fill = ctap)) +
  geom_bar(stat = "identity", position = position_dodge2(preserve = "single", padding = 0), width = 0.7) +
  scale_x_discrete(labels = function(x) {
    # Add a line break after the module number for a hierarchical appearance
    sapply(strsplit(x, "\\."), function(x) paste(x, collapse = "\n"))
  }) +
  labs(
    title = "AUC Values > 0.65 by Cell Type, Module, and Treatment",
    x = "Module and Cell Type",
    y = "AUC",
    fill = "Treatment"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.spacing = unit(1, "lines"),
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank()   # Remove minor grid lines
  ) +
  coord_cartesian(ylim = c(0.5, 0.9))  # Apply Cartesian filter to focus on AUC > 0.5

# Print the plot
print(pctap)

ggsave("figure6_CTAPs.pdf", pctap, width = 8, height = 6, dpi = 300)


roc_curve_list <- readRDS("CTAPs/roc_curve_list_CTAP.rds")


############################### CTAP ROC curves#################################
################################################################################
#
# ## With confidence inte
#
# rm(list = ls())
# setwd("/ix/djishnu/Javad/RA/")
#
# library(pROC)
# library(ggplot2)
# library(dplyr)
# library(ggpubr)
# library(cowplot)
#
# # Load predictions
# predictions <- readRDS("CTAPs/predictions_CTAP.rds")
#
# # Function to extract metadata from names
# extract_info <- function(name) {
#   parts <- unlist(strsplit(name, "_"))
#   celltype <- parts[1]
#   module <- as.numeric(parts[2])
#   CTAP <- parts[3]
#   return(data.frame(name, celltype, module, CTAP, stringsAsFactors = FALSE))
# }
#
# # Create metadata dataframe
# metadata <- do.call(rbind, lapply(names(predictions), extract_info))
#
# # Define a vector of line types
# line_types <- c("solid", "dashed", "dotdash", "dotted", "longdash", "twodash")
#
# # Container for the final arranged plots for each CTAP
# plot_list <- list()
#
# # Loop over each CTAP
# for (ctap in unique(metadata$CTAP)) {
#   ctap_data <- metadata %>% filter(CTAP == ctap)
#   cell_plots <- list()
#
#   # Loop over each cell type within the current CTAP
#   for (cell in unique(ctap_data$celltype)) {
#     cat("Processing celltype |", cell, "in CTAP |", ctap, "\n")
#     cell_data <- ctap_data %>% filter(celltype == cell)
#     roc_list <- list()
#     auc_list <- data.frame()
#     legend_labels <- c()
#
#     # Compute ROC curves for each prediction
#     for (i in 1:nrow(cell_data)) {
#       name <- cell_data$name[i]
#       module <- cell_data$module[i]
#
#       pred <- predictions[[name]]
#       roc_curve <- roc(response = pred$pred$pred$obs, predictor = pred$pred$pred[,6])
#       auc_value <- as.numeric(auc(roc_curve))
#
#       # Compute the AUC confidence interval
#       ci_auc <- ci.auc(roc_curve)
#
#       # Check if the ROC curve goes below the y = x line (TPR < FPR)
#       roc_df <- data.frame(
#         TPR = rev(roc_curve$sensitivities),
#         FPR = rev(1 - roc_curve$specificities),
#         module = module
#       )
#
#       # Check if any TPR < FPR
#       if (any(roc_df$TPR < roc_df$FPR)) {
#         cat("Skipping module", module, "because the ROC curve goes below the y = x line.\n")
#         next
#       }
#
#       # Only include modules that satisfy:
#       # (1) AUC > 0.65 and (2) lower bound of AUC CI > 0.5
#       if (auc_value > 0.65 && ci_auc[1] > 0.5) {
#         # Calculate pointwise confidence intervals for the ROC curve (sensitivities)
#         ci <- ci.se(roc_curve, specificities = seq(0, 1, 0.01))
#
#         ci_df <- data.frame(
#           FPR = 1 - attr(ci, "specificities"),
#           TPR.lower = ci[, 1],
#           TPR.upper = ci[, 3],
#           module = module
#         )
#
#         roc_list[[i]] <- list(roc_df = roc_df, ci_df = ci_df)
#
#         auc_list <- rbind(auc_list, data.frame(module = module, AUC = auc_value))
#
#         legend_labels <- c(legend_labels, paste("Module", module, "(AUC:", round(auc_value, 3), ")"))
#       } else {
#         cat("Skipping module", module, "because AUC =", round(auc_value, 3),
#             ", lower CI =", round(ci_auc[1], 3), "\n")
#       }
#     }
#
#     # Skip cell types that did not meet the threshold
#     if (length(roc_list) == 0) next
#
#     # Combine all ROC data for this cell type
#     roc_data <- do.call(rbind, lapply(roc_list, function(x) x$roc_df))
#     ci_data <- do.call(rbind, lapply(roc_list, function(x) x$ci_df))
#
#     # Define colors for each cell type (adjust as needed)
#     colors <- c("pink", "blue", "green", "orange")
#     names(colors) <- unique(ctap_data$celltype)
#     roc_data$color <- colors[cell]
#     ci_data$color <- colors[cell]
#
#     # Ensure colors is a named vector for modules
#     names(colors) <- unique(roc_data$module)
#
#     # Build the ROC plot with confidence intervals
#     p <- ggplot(roc_data, aes(x = FPR)) +
#       geom_line(aes(y = TPR, linetype = factor(module), color = factor(module)), size = 1, show.legend = TRUE) +
#       geom_ribbon(data = ci_data, aes(ymin = TPR.lower, ymax = TPR.upper, fill = factor(module)), alpha = 0.2) +
#       geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red", size = 1) +  # y = x line
#       scale_linetype_manual(values = line_types, name = "Modules") +
#       scale_color_manual(values = colors, name = "Modules") +
#       scale_fill_manual(values = colors, name = "Confidence Interval") +
#       labs(
#         title = paste(cell),
#         x = "False Positive Rate",
#         y = "True Positive Rate"
#       ) +
#       theme(
#         panel.background = element_blank(),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         axis.line = element_line(color = "black"),
#         axis.ticks = element_line(color = "black"),
#         legend.position = "right",
#         axis.text.x = element_text(family = "Times", size = 20),
#         axis.text.y = element_text(family = "Times", size = 20)
#       )
#
#     cell_plots[[cell]] <- p
#   }
#
#   # If no valid plots were created for this CTAP, skip to the next one
#   if (length(cell_plots) == 0) next
#
#   # Arrange the cell type plots in a 2x2 grid using ggpubr
#   arranged_plots <- ggarrange(plotlist = cell_plots, ncol = 2, nrow = 2)
#
#   # Extract a common legend from one of the cell plots.
#   common_legend <- NULL
#   for (plot_name in names(cell_plots)) {
#     tmp_leg <- get_legend(cell_plots[[plot_name]])
#     if (!is.null(tmp_leg) && length(tmp_leg$grobs) > 0) {
#       common_legend <- tmp_leg
#       break
#     }
#   }
#
#   if (is.null(common_legend)) {
#     dummy_df <- data.frame(module = factor(auc_list$module, levels = unique(auc_list$module)))
#     dummy_plot <- ggplot(dummy_df, aes(x = module, y = module, linetype = module)) +
#       geom_line() +
#       scale_linetype_manual(values = line_types, drop = FALSE) +
#       guides(linetype = guide_legend(override.aes = list(color = "black"))) +
#       theme(legend.position = "right")
#     common_legend <- get_legend(dummy_plot)
#   }
#
#   final_plot <- plot_grid(arranged_plots, common_legend, ncol = 2, rel_widths = c(3, 1))
#   plot_list[[ctap]] <- final_plot
#
#   ggsave(final_plot, file = paste0("figures/figure6/", ctap, ".pdf"), width = 20, height = 15)
#   print(final_plot)
# }

################################################################################

rm(list = ls())
setwd("/ix/djishnu/Javad/RA/")

library(pROC)
library(ggplot2)
library(dplyr)
library(ggpubr)
library(cowplot)
library(RColorBrewer)
library(grid)  # for unit()

# Load predictions
predictions <- readRDS("CTAPs/predictions_CTAP.rds")

# Function to extract metadata from prediction names
extract_info <- function(name) {
  parts <- unlist(strsplit(name, "_"))
  celltype <- parts[1]
  module <- as.numeric(parts[2])
  CTAP <- parts[3]
  return(data.frame(name, celltype, module, CTAP, stringsAsFactors = FALSE))
}

# Create a metadata dataframe
metadata <- do.call(rbind, lapply(names(predictions), extract_info))

# Define line types for different modules
line_types <- c("dashed", "dashed", "dashed", "dashed", "dashed", "dashed")

# Container for the final arranged plots for each CTAP
plot_list <- list()

# Loop over each CTAP
for (ctap in unique(metadata$CTAP)) {
  ctap_data <- metadata %>% filter(CTAP == ctap)
  cell_plots <- list()

  # Loop over each cell type within the current CTAP
  for (cell in unique(ctap_data$celltype)) {
    cat("Processing celltype |", cell, "in CTAP |", ctap, "\n")
    cell_data <- ctap_data %>% filter(celltype == cell)
    roc_list <- list()
    auc_list <- data.frame()
    legend_labels <- c()

    # Compute ROC curves for each prediction in this cell type
    for (i in 1:nrow(cell_data)) {
      name <- cell_data$name[i]
      module <- cell_data$module[i]

      pred <- predictions[[name]]
      roc_curve <- roc(response = pred$pred$pred$obs, predictor = pred$pred$pred[,6])
      auc_value <- as.numeric(auc(roc_curve))

      if (auc_value > 0.65) {
        roc_df <- data.frame(
          TPR = rev(roc_curve$sensitivities),
          FPR = rev(1 - roc_curve$specificities),
          module = module
        )
        roc_list[[i]] <- roc_df
        auc_list <- rbind(auc_list, data.frame(module = module, AUC = auc_value))
        legend_labels <- c(legend_labels, paste("Module", module, "(AUC:", round(auc_value, 3), ")"))
      }
    }

    # Skip cell types with no valid ROC curves
    if (length(roc_list) == 0) next

    # Combine all ROC data for this cell type
    roc_data <- do.call(rbind, roc_list)

    # Determine the unique modules (as characters) in this plot
    module_keys <- unique(as.character(roc_data$module))
    num_module_keys <- length(module_keys)

    # Create a color palette for the modules.
    # Use at least 3 colors if num_module_keys < 3.
    if (num_module_keys <= 8) {
      palette_size <- max(num_module_keys, 3)
      raw_colors <- brewer.pal(n = palette_size, name = "Dark2")[1:num_module_keys]
    } else {
      raw_colors <- rainbow(num_module_keys)
    }
    module_colors <- setNames(raw_colors, module_keys)
    # Override Module 1's color if it exists (set to red)
    if ("1" %in% module_keys) {
      module_colors["1"] <- "red"
    }
    # Set color for the perfect classifier to dark green
    module_colors["Perfect classifier"] <- "darkgreen"

    # Map module keys to linetypes using the provided line_types vector.
    module_line_types <- setNames(line_types[seq_len(num_module_keys)], module_keys)
    module_line_types["Perfect classifier"] <- "dashed"

    # Data for the perfect classifier (ideal ROC curve: (0,0) → (0,1) → (1,1))
    perfect_data <- data.frame(
      FPR = c(0, 0, 1),
      TPR = c(0, 1, 1)
    )

    # Ensure the module column is treated as a factor
    roc_data$module <- as.factor(roc_data$module)

    # Build the plot:
    p <- ggplot() +
      # Plot ROC curves for the modules; map color and linetype to module
      geom_line(data = roc_data,
                aes(x = FPR, y = TPR, group = module,
                    color = factor(module), linetype = factor(module)),
                size = 1) +
      # Add the perfect classifier ROC curve (using dark green)
      geom_line(data = perfect_data,
                aes(x = FPR, y = TPR,
                    color = "Perfect classifier", linetype = "Perfect classifier"),
                size = 1) +
      scale_color_manual(
        name = "Module (AUC)",
        values = module_colors,
        breaks = c(module_keys, "Perfect classifier"),
        labels = c(legend_labels, "Perfect classifier"),
        drop = FALSE
      ) +
      scale_linetype_manual(
        name = "Module (AUC)",
        values = module_line_types,
        breaks = c(module_keys, "Perfect classifier"),
        labels = c(legend_labels, "Perfect classifier"),
        drop = FALSE
      ) +
      geom_abline(slope = 1, intercept = 0, linetype = "dotted",
                  color = "gray", size = 1) +
      coord_fixed(ratio = 1) +
      labs(title = paste("ROC Curve for", cell, "in", ctap),
           x = "False Positive Rate",
           y = "True Positive Rate") +
      theme_classic() +
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_rect(fill = "white", color = NA),
            axis.line = element_line(color = "black"),
            axis.ticks = element_line(color = "black"),
            legend.background = element_rect(color = "black", fill = "white", linewidth = 0.5),
            legend.key.size = unit(1, "lines"),
            legend.margin = margin(15, 15, 15, 15),
            legend.text = element_text(size = 12, family = "Times"),
            legend.title = element_text(size = 12, family = "Times"),
            legend.position =c(0.8, 0.2),
            axis.text.x = element_text(family = "Times", size = 20),
            axis.text.y = element_text(family = "Times", size = 20),
            axis.title.x = element_text(family = "Times", size = 25, face = "plain"),
            axis.title.y = element_text(family = "Times", size = 25, face = "plain")
      )


  #print(p)

  ggsave(p,filename = paste0("/ix/djishnu/Javad/RA/figures/figure6/",ctap,"/",cell,".pdf"),width = 8,height = 8)
    # Save the plot for this cell type
    cell_plots[[cell]] <- p
  }


}


rm(list = ls())
setwd("/ix/djishnu/Javad/RA/")

library(pROC)
library(ggplot2)
library(dplyr)
library(ggpubr)
library(cowplot)
library(RColorBrewer)
library(grid) # for unit()

# Load predictions
predictions <- readRDS("CTAPs/predictions_CTAP.rds")

# Function to extract metadata from prediction names
extract_info <- function(name) {
  parts <- unlist(strsplit(name, "_"))
  celltype <- parts[1]
  module <- as.numeric(parts[2])
  CTAP <- parts[3]
  return(data.frame(name, celltype, module, CTAP, stringsAsFactors = FALSE))
}

# Create a metadata dataframe
metadata <- do.call(rbind, lapply(names(predictions), extract_info))

# Define line types for different modules
line_types <- c("dashed", "dashed", "dashed", "dashed", "dashed", "dashed")

# Container for the final arranged plots for each CTAP
plot_list <- list()

# Loop over each CTAP
for (ctap in unique(metadata$CTAP)) {
  ctap_data <- metadata %>% filter(CTAP == ctap)
  cell_plots <- list()

  # Loop over each cell type within the current CTAP
  for (cell in unique(ctap_data$celltype)) {
    cat("Processing celltype |", cell, "in CTAP |", ctap, "\n")
    cell_data <- ctap_data %>% filter(celltype == cell)
    roc_list <- list()
    auc_list <- data.frame()
    legend_labels <- c()

    # Compute ROC curves for each prediction in this cell type
    for (i in 1:nrow(cell_data)) {
      name <- cell_data$name[i]
      module <- cell_data$module[i]

      pred <- predictions[[name]]
      roc_curve <- roc(response = pred$pred$pred$obs, predictor = pred$pred$pred[,6])
      auc_value <- as.numeric(auc(roc_curve))

      if (auc_value > 0.65) {
        roc_df <- data.frame(
          TPR = rev(roc_curve$sensitivities),
          FPR = rev(1 - roc_curve$specificities),
          module = module
        )
        roc_list[[i]] <- roc_df
        auc_list <- rbind(auc_list, data.frame(module = module, AUC = auc_value))
        legend_labels <- c(legend_labels, paste("Module", module, "(AUC:", round(auc_value, 3), ")"))
      }
    }

    # Skip cell types with no valid ROC curves
    if (length(roc_list) == 0) next

    # Combine all ROC data for this cell type
    roc_data <- do.call(rbind, roc_list)

    # Determine the unique modules (as characters) in this plot
    module_keys <- unique(as.character(roc_data$module))
    num_module_keys <- length(module_keys)

    # Create a color palette for the modules using Okabe-Ito palette.
    # Okabe-Ito palette is colorblind-friendly and available in base R (v4.0.0+)
    raw_colors <- palette.colors(n = num_module_keys, palette = "Okabe-Ito")
    module_colors <- setNames(raw_colors, module_keys)

    # Explicitly set colors for specific modules and perfect classifier
    if ("1" %in% module_keys) {
      module_colors["1"] <- "#000000" # Black for Module 1
    }
    if ("2" %in% module_keys) {
      module_colors["2"] <- "#0072B2" # Blue for Module 2 (Okabe-Ito)
    }
    if ("5" %in% module_keys) {
      module_colors["5"] <- "#cc79a7" # Cyan/Green for Module 5 (Okabe-Ito)
    }
    if ("10" %in% module_keys) {
      module_colors["10"] <- "#E69F00" # Orange for Module 10 (Okabe-Ito)
    }
    module_colors["Perfect classifier"] <- "darkgreen" # Dark Green for Perfect classifier


    # Map module keys to linetypes using the provided line_types vector.
    module_line_types <- setNames(line_types[seq_len(num_module_keys)], module_keys)
    module_line_types["Perfect classifier"] <- "dashed"

    # Data for the perfect classifier (ideal ROC curve: (0,0) → (0,1) → (1,1))
    perfect_data <- data.frame(
      FPR = c(0, 0, 1),
      TPR = c(0, 1, 1)
    )

    # Ensure the module column is treated as a factor
    roc_data$module <- as.factor(roc_data$module)

    # Build the plot:
    p <- ggplot() +
      # Plot ROC curves for the modules; map color and linetype to module
      geom_line(data = roc_data,
                aes(x = FPR, y = TPR, group = module,
                    color = factor(module), linetype = factor(module)),
                size = 1.5) + # Increased line size to 1.5
      # Add the perfect classifier ROC curve (using dark green)
      geom_line(data = perfect_data,
                aes(x = FPR, y = TPR,
                    color = "Perfect classifier", linetype = "Perfect classifier"),
                size = 1.5) + # Increased line size to 1.5
      scale_color_manual(
        name = "Module (AUC)",
        values = module_colors,
        breaks = c(module_keys, "Perfect classifier"),
        labels = c(legend_labels, "Perfect classifier"),
        drop = FALSE
      ) +
      scale_linetype_manual(
        name = "Module (AUC)",
        values = module_line_types,
        breaks = c(module_keys, "Perfect classifier"),
        labels = c(legend_labels, "Perfect classifier"),
        drop = FALSE
      ) +
      geom_abline(slope = 1, intercept = 0, linetype = "dotted",
                  color = "gray", size = 1) +
      coord_fixed(ratio = 1) +
      labs(title = paste("ROC Curve for", cell, "in", ctap),
           x = "False Positive Rate",
           y = "True Positive Rate") +
      theme_classic() +
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_rect(fill = "white", color = NA),
            axis.line = element_line(color = "black"),
            axis.ticks = element_line(color = "black"),
            legend.background = element_rect(color = "black", fill = "white", linewidth = 0.5),
            legend.key.size = unit(1, "lines"),
            legend.margin = margin(15, 15, 15, 15),
            legend.text = element_text(size = 12, family = "Times"),
            legend.title = element_text(size = 12, family = "Times"),
            legend.position =c(0.8, 0.2),
            axis.text.x = element_text(family = "Times", size = 20),
            axis.text.y = element_text(family = "Times", size = 20),
            axis.title.x = element_text(family = "Times", size = 25, face = "plain"),
            axis.title.y = element_text(family = "Times", size = 25, face = "plain")
      )


    print(p)

    ggsave(p,filename = paste0("/ix/djishnu/Javad/RA/figures/figure6/",ctap,"/",cell,".pdf"),width = 8,height = 8)
    # Save the plot for this cell type
    cell_plots[[cell]] <- p
  }


}
