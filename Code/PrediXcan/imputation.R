library(caret)
library(glmnet)
library(gbm)
library(matrixStats)
library(cvAUC)
library(pROC)
library(rstatix)
library(ggpubr)
library(ggplot2)
library(data.table)
library(tidyverse)
library(coefplot)
library(randomForest)
library(crossV)
library(ROCR)

library(doParallel)


cores <- 10
registerDoParallel(cores = cores)


get_module_genes_lasso = function(fn_data, fn_endotype, module_genes_fn, lambda_val) {
  
  fn_folds <- createFolds(y = fn_endotype, k=10, list = FALSE, returnTrain = FALSE)
  X = data.matrix(fn_data[,names(fn_data) %in% module_genes_fn ])
  
  if (ncol(X)==0) {
    return(c())
  }
  
  if (ncol(X)==1) {
    # return(c(names(data)[which(names(data) %in% module_genes_fn)]))
    return(c())
  }
  
  selected_features_list = c()
  
  # Perform 10-fold cross-validation with feature selection for classification
  for (i in 1:10) {
    
    fold_fn = which(fn_folds == i)
    
    x_train_fn = fn_data[-fold_fn,(names(fn_data) %in% module_genes_fn)]
    # y_train_fn = as.factor(fn_endotype[-fold_fn])
    y_train_fn = fn_endotype[-fold_fn]
    
    # x_test_fn = data[fold_fn,(names(fn_data) %in% module_genes_fn)]
    # y_test_fn = as.factor(fn_endotype[fold])
    
    
    # Fit the Lasso model with the optimal lambda
    lasso_model_fn <- glmnet(x = x_train_fn, y = y_train_fn, alpha = 1, lambda = lambda_val, family = "binomial")
    # lasso_model <- glmnet(x = X, y = as.factor(endotype), alpha = 1, lambda = optimal_lambda, family = "binomial")
    
    # Get the selected features for this run
    coefs <- coef(lasso_model_fn)
    inds<-which(coefs!=0)
    
    #remove intercept from the list
    selected_features <- tail(row.names(coefs)[inds], -1)
    
    
    # Append the selected features to the list
    selected_features_list <- c(selected_features_list,selected_features)
  }
  
  # Use the table function to find the counts of each feature
  feature_counts <- table(selected_features_list)
  
  # Set a threshold (e.g., at least 8 out of 10 runs)
  threshold <- 1
  
  # Filter features based on the threshold
  selected_features_above_threshold <- names(feature_counts[feature_counts >= threshold])

  # return
  if(!is.null(selected_features_above_threshold)){
    selected_features_above_threshold
  }
  else {
    list()
  }
    
  
  
}




# Split dataset into 300/600



celltype_filenames = list.files(".",pattern = "*_predicted_expression_gene_mapped.tsv")
famFile = read.delim("racer-mega-imputed-f-matched-v3-QC.fam",sep = " ",header = FALSE)
endotype = famFile[,6]-1 # -1 to convert from {1,2} to {0,1} for glm

ldak_modules = readRDS("significant_modules_members.RDS")

idx_300 = sample.int(length(endotype),size = 313)

write.csv(idx_300,"idx_300.csv",quote = FALSE,row.names = FALSE, col.names = FALSE)

table(endotype[idx_300])
table(endotype[-idx_300])

master_feature_list = list()
for (lambda in 1:10) {
  master_feature_list[[lambda]] = list()
  for (celltype_filename in celltype_filenames) {
    master_feature_list[[lambda]][[celltype_filename]] = list()
    for (module_idx in 1:length(ldak_modules)) {
      master_feature_list[[lambda]][[celltype_filename]][[module_idx]] = list()
    }
  }  
}


for (celltype_filename in celltype_filenames) {
  loop_data = read.delim(celltype_filename)[idx_300,]
  loop_endotype = endotype[idx_300]
  
  for (lambda in 1:10) {
    
    for (module_idx in 1:length(ldak_modules)) {
      
      
      # for each tissue type, for each lambda, for each module,
      # select features using lasso with a specific lambda value
      
      module_idx_genes = ldak_modules[[module_idx]]
      master_feature_list[[lambda]][[celltype_filename]][[module_idx]] = get_module_genes_lasso(loop_data,
                                                                                                loop_endotype,
                                                                                                module_idx_genes,
                                                                                                lambda*0.01)
    }
  }
}



master_preds_list = list()
master_auc_list = list()
master_y_list = list()
folds <- createFolds(y = endotype[-idx_300], k=10, list = FALSE, returnTrain = FALSE)

for (lambda in 1:10) {
  master_preds_list[[lambda]] = list()
  master_auc_list[[lambda]] = list()
  master_y_list[[lambda]] = list()
  for (celltype_filename in celltype_filenames) {
    master_preds_list[[lambda]][[celltype_filename]] = list()
    master_auc_list[[lambda]][[celltype_filename]] = list()
    master_y_list[[lambda]][[celltype_filename]] = list()
    
    
    for (module_idx in 1:length(ldak_modules)) {
      master_auc_list[[lambda]][[celltype_filename]][[module_idx]] = 0.5
      master_preds_list[[lambda]][[celltype_filename]][[module_idx]] = list()
      master_y_list[[lambda]][[celltype_filename]][[module_idx]] = list()
    }
  }
}



for (lambda in 1:10) {
  
  for (celltype_filename in celltype_filenames) {
    
    loop_data = read.delim(celltype_filename)[-idx_300,]
    loop_endotype = endotype[-idx_300]
    
    for (module_idx in 1:length(ldak_modules)) {
      
      
      module_genes = master_feature_list[[lambda]][[celltype_filename]][[module_idx]]
      # master_preds_list[[lambda]][[celltype_filename]][[module_idx]] = 0
      # master_y_list[[lambda]][[celltype_filename]][[module_idx]] = 0
      # master_preds_list[[lambda]][[celltype_filename]][[module_idx]] = c()
      # master_y_list[[lambda]][[celltype_filename]][[module_idx]] = c()
      
      if(ncol(loop_data[,(names(loop_data) %in% module_genes)])==0) {
        next
      }
      
      for (i in 1:10) {
        
        fold = which(folds == i)
        
        x_train = loop_data[-fold,(names(loop_data) %in% module_genes)]
        y_train = as.factor(loop_endotype[-fold])
        
        x_test = loop_data[fold,(names(loop_data) %in% module_genes)]
        y_test = as.factor(loop_endotype[fold])
        
        # lasso_model <- glmnet(x = x_train, y = y_train, lambda = 0, family = "binomial")
        # 
        # preds = predict(lasso_model, data.matrix(x_test))
        # preds = predict(lasso_model, data.matrix(x_train))
        
        # yhat.RF = predict(RFfit, x_test)
        
        RFfit <- randomForest(x = x_train, y = y_train, importance=TRUE)
        preds = predict(RFfit, x_test, type = "prob")[,2]
        
        master_preds_list[[lambda]][[celltype_filename]][[module_idx]] = c(master_preds_list[[lambda]][[celltype_filename]][[module_idx]],
                                                                           preds)
        master_y_list[[lambda]][[celltype_filename]][[module_idx]] = c(master_y_list[[lambda]][[celltype_filename]][[module_idx]],
                                                                       as.numeric(as.character(y_test)))
        
        
      }
      
      roc_pred = prediction(unlist(master_preds_list[[lambda]][[celltype_filename]][[module_idx]]) , 
                            unlist(master_y_list[[lambda]][[celltype_filename]][[module_idx]]))
      perf = performance(roc_pred, "tpr", "fpr")
      master_auc_list[[lambda]][[celltype_filename]][[module_idx]] = performance(roc_pred,measure = "auc")@y.values[[1]]
      pdf(paste("./roc_plots_named/", paste(lambda*0.01, celltype_filename, module_idx,"roc.pdf",sep = "_"),sep = ""))
      plot(perf)
      title(performance(roc_pred,measure = "auc")@y.values[[1]])
      dev.off()
      
    }
  }  
}





for (celltype_filename in celltype_filenames) {
  loop_data = read.delim(celltype_filename)[-idx_300,]
  loop_endotype = endotype[-idx_300]
  
  for (module_idx in 1:length(ldak_modules)) {
    
    for (lambda in 1:10) {
      
      module_genes = master_feature_list[[lambda]][[celltype_filename]][[module_idx]]
      
      for (i in 1:10) {
        
        fold = which(folds == i)
        
        x_train = loop_data[-fold,(names(loop_data) %in% module_genes)]
        y_train = as.factor(loop_endotype[-fold])
        
        x_test = loop_data[fold,(names(loop_data) %in% module_genes)]
        y_test = as.factor(loop_endotype[fold])
        
        
        
        
      }
      
      module_idx_genes = ldak_modules[[module_idx]]
      master_feature_list[[lambda]][[celltype_filename]] = get_module_genes_lasso(loop_data,
                                                                                  loop_endotype,
                                                                                  module_idx_genes,
                                                                                  lambda)
    }
  }
}


