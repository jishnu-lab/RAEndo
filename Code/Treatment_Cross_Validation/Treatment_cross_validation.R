
##################################################################
##                      Plot for treatment                      ##
##################################################################

setwd("/ix/djishnu/Javad/RA/")
# Load necessary libraries
library(caret)
library(glmnet)
library(pROC)
library(dplyr)

# Load your data
modules <- readRDS("Data/significant_modules_members.RDS")
RA_meta <- readRDS("CTAPs/RAmetaData.rds")
psudobulks <- readRDS("CTAPs/unmathced_psudobulks.rds")  # Assuming you have a list of pseudobulk matrices for each cell type

# Define cell types (assuming you have a list of cell types)
celltypes <- names(psudobulks)

# Initialize a data frame to store AUC results
auc_results <- data.frame(
  celltype = character(),
  module = integer(),
  treatment = character(),
  auc = numeric(),
  stringsAsFactors = FALSE
)

# Initialize a list to store predictions
predictions <- list()

# Loop through each cell type
for (celltype in celltypes) {
  
  # Prepare the data
  psudobulk <- psudobulks[[celltype]]
  y <- RA_meta[match(rownames(psudobulk), RA_meta$donor), c("donor", "ra_group")]
  y$ra_group <- gsub(" ", "\\.", y$ra_group)
  
  # Create a table of donor vs. treatment
  ytab <- table(y$donor, y$ra_group)
  
  # Loop through each treatment
  for (tr in colnames(ytab)) {
    
    # Create a binary factor for the current treatment
    y_tr <- ifelse(y$ra_group == tr, tr, "rest")
    y_tr <- factor(y_tr, levels = c(tr, "rest"))
    
    # Loop through each module
    for (module_name in c(1, 2, 5, 10)) {
      # Check if module exists
      
      
      # Extract features for the current module
      module_genes <- modules[[module_name]]
      X <- psudobulk[, colnames(psudobulk) %in% module_genes, drop = FALSE]
      
      # Ensure X and y are matched
      if (nrow(X) != length(y_tr)) next
      
      # Train the model using glmnet
      set.seed(123)
      ctrl <- trainControl(method = "cv", number = 10, classProbs = TRUE, summaryFunction = twoClassSummary)
      model <- train(x = X, y = y_tr, method = "glmnet", trControl = ctrl, metric = "ROC", tuneGrid = expand.grid(alpha = 0.5, lambda = 0.01))
      
      # Get the cross-validated AUC
      auc_value <- max(model$results$ROC)
      
      # Save results in the data frame
      auc_results <- rbind(auc_results, data.frame(
        celltype = celltype,
        module = module_name,
        treatment = tr,
        auc = auc_value
      ))
      
      # Save predictions for ROC curve
      pred <- predict(model, X, type = "prob")
      predictions[[paste0(celltype, "_", module_name, "_", tr)]] <- pred
      
      # Print progress
      cat("Processed cell type:", celltype, "| Treatment:", tr, "| Module:", module_name, "| AUC:", auc_value, "\n")
    }
  }
}

# Save AUC results
saveRDS(auc_results, "auc_results_treatment.rds")

# Save predictions for ROC curves
saveRDS(predictions, "predictions_treatment.rds")





# Example of how to plot ROC curves using the saved predictions
for (name in names(predictions)) {
  pred <- predictions[[name]]
  roc_curve <- roc(response = y_tr, predictor = pred$case)  # Assuming 'case' is the positive class
  plot(roc_curve, main = paste("ROC Curve for", name))
}
