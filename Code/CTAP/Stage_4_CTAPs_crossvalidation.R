##################################################################
##                      Plot for CTAPs.                         ##
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


# Initialize a list to store predictions
predictions <- list()
auc_results <- list()
selected_variables <- list()

# Loop through each cell type

set.seed(12345)
roc_curve_list <- list()
for (celltype in celltypes) {

  # Prepare the data
  psudobulk <- psudobulks[[celltype]]
  y <- RA_meta[match(rownames(psudobulk), RA_meta$donor), c("donor", "CTAP")]
  #y$ra_group <- gsub(" ", "\\.", y$ra_group)

  # Create a table of donor vs. treatment
  ytab <- table(y$donor, y$CTAP)

  # Loop through each treatment
  for (ct in colnames(ytab)) {

    # Create a binary factor for the current treatment
    y_ct <- ifelse(y$CTAP == ct, ct, "rest")
    y_ct <- factor(y_ct, levels = c(ct, "rest"))

    # Loop through each module
    for (module_name in c(1, 2, 5, 10)) {
      # Check if module exists


      # Extract features for the current module
      module_genes <- modules[[module_name]]
      X <- psudobulk[, colnames(psudobulk) %in% module_genes, drop = FALSE]

      # Ensure X and y are matched
      if (nrow(X) != length(y_ct)) next

      # Train the model using glmnet
      set.seed(123)
      ctrl <- trainControl(method = "repeatedcv", number = 5, repeats = 1, classProbs = TRUE, summaryFunction = twoClassSummary,savePredictions = "final")
      model <- train(
        x = X,
        y = y_ct,
        method = "glmnet",
        trControl = ctrl,
        metric = "ROC",
        tuneGrid = expand.grid(alpha = 1, lambda = seq(0.01, 0.05, length = 10)),
        glmnet.control = list(maxit = 5000000)
      )
      # Get the cross-validated AUC
      #auc_value <- max(model$results$ROC)

      roc_curve<- roc(model$pred$obs,model$pred[,6])
      roc_curve_list[[ct]]<- list()
      roc_curve_list[[ct]][[celltype]] <- list()
      roc_curve_list[[ct]][[celltype]][[as.character(module_name)]] <- roc_curve
      auc_value <- auc(roc_curve)


      ## selected variables
      # Extract the best model's coefficients
      best_model <- model$finalModel
      best_lambda <- model$bestTune$lambda
      coefficients <- coef(best_model, s = best_lambda)

      # Get the selected variables (non-zero coefficients)
      selected_vars <- rownames(coefficients)[which(coefficients != 0)]
      selected_vars <- selected_vars[selected_vars != "(Intercept)"]  # Remove intercept if present
      # Save selected variables to the list
      selected_variables[[paste0(celltype, "_", module_name, "_", ct)]] <- selected_vars


      # Save results in the data frame
      auc_results <- rbind(auc_results, data.frame(
        celltype = celltype,
        module = module_name,
        ctap = ct,
        auc = auc_value
      ))


      # Save both predictions and actual labels
      predictions[[paste0(celltype, "_", module_name, "_", ct)]] <- list(pred = model, y_ct = y_ct)
      # Print progress
      cat("Processed cell type:", celltype, "| CTAP:", ct, "| Module:", module_name, "| AUC:", auc_value, "\n")
    }
  }
}


# Save ROC curve results
saveRDS(roc_curve_list,"CTAPs/roc_curve_list_CTAP.rds")

# Save AUC results
saveRDS(auc_results, "CTAPs/auc_results_CTAP.rds")

# Save predictions for ROC curves
saveRDS(predictions, "CTAPs/predictions_CTAP.rds")


# # Example of how to plot ROC curves using the saved predictions
# for (name in names(predictions)) {
#   pred <- predictions[[name]]
#   roc_curve <- roc(response = predictions[[name]]$y_ct, predictor = predictions[[name]]$pred[,1])  # Assuming 'case' is the positive class
#   plot(roc_curve, main = paste("ROC Curve for", name))
# }












