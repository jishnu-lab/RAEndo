# This file get the emperical pvalues from a loss model on the genes of the modules.
# The procedure is using caret for the analysis





library(biomaRt)
library(tidyverse)
library(dplyr)
library(doParallel)
library(foreach)
library(caret)
library(tidyverse)





################################################################################
#                             Reading Data
################################################################################

setwd("/ix/djishnu/Javad/RA/")
rm(list = ls())
bashInp <- commandArgs(trailingOnly = T)

source("Data_Analysis/Comparison_scripts/Funcs.R")





print(paste0("args is: \n"))
print(paste0("The input is ",bashInp))

expression_mt <- read.table('./Data/low_input_gene_sample_tpm_matrix.725714.tsv')

metadata_df <- read.table('./Data/metadata_for_bulk_RNA_seq.tsv')



##################################################################
##                    Maping ensmbl to genes                    ##
##################################################################



cell_specific_gene_dir <- './Data/Cell_specific_gene_list/'

module_dir <- './Data'

module_file <- 'significant_modules_members.RDS'

module_set <- readRDS(paste0(module_dir, '/', module_file))

len <- sapply(module_set, length)

output_directory <- '/ix/djishnu/Javad/RA/Data_Analysis/results/July22Results'

################################################################################
#                              Forming LDAK genes
################################################################################
module_set_genes <- unlist(module_set)

genes_all        <- read.table('./Data/LDAK_score.csv')

noModule         <- genes_all[!genes_all$V1 %in% module_set_genes,]

ldak_genes_top   <- noModule[noModule$V2 >= 2, ]

moduleandtopLDAK <- c(module_set_genes,ldak_genes_top$V1)

randomgeneall    <- genes_all[!genes_all$V1 %in% moduleandtopLDAK,]

################################################################################
#                              Running Lasso
################################################################################


#cell_types          <- c('B cell','Mono', 'T cell',"Fibro")

cell_types          <- c("Fibro")

#module_types <- c("real_module","random_module","dist_module")
module_types <- c("real_module","dist_module","random_module")

cat(paste0("Processing cell",cell_types[as.numeric(bashInp)],"... /n"))





pvalReport <- NULL

for(cell in cell_types[as.numeric(bashInp)]){

for( k in 1:length(module_set)){

report <- NULL

pval_df <- NULL

for(typ in module_types){

for(rep in 1:100){


     cell_gene_list <- read.table(paste0(cell_specific_gene_dir, '/', cell, '_gene_list_top_75.txt'), fill = T, header = T)

     cell_data  <- create_cell_data(cell = cell,
                                    metadata=metadata_df,
                                    expression_matrix=expression_mt,
                                    cell_specific_genes=cell_gene_list)




     if(typ=="real_module"){



       cat(paste0("processing module ",k,"\n\n"))
       cat(paste0("processing cell type ",cell,"\n"))



    cellmoData      <- create_module_specific_dataset(module = module_set[[k]],
                                                   expression_data = cell_data$cell_expression_data,
                                                cell_type_metadata = cell_data$cell_metadata,
                                                              cell = cell) ## This is the data for a specific cell and specific module
    if(class(cellmoData)=="numeric"){

      print("Yes")

      next

      }

     COLS <- colnames(cellmoData$data)
     X    <- cellmoData$data[,COLS!="Phenotype"]

     modulesgene <- colnames(X)

     y    <- cellmoData$data[,COLS=="Phenotype",drop=F]

     data <- data.frame(y=factor(y$Phenotype,levels=c(0,1),labels = c("X0","X1")),X)


     # Custom summary function






     myControl <- caret::trainControl(
       method = "LOOCV",
       summaryFunction = customSummary,
       classProbs = TRUE, # IMPORTANT!
       selectionFunction="best"
       )

     # After creating a custom trainControl object, it's easy to fit caret models that use AUC rather than accuracy to tune and evaluate the model. You can just pass your custom trainControl object to the train() function via the trControl argument.

     tuneGrid <- expand.grid(
       alpha = 1,
       lambda = seq(0.001, 1, length = 100) # Adjust the lambda sequence as needed
     )



     # Train glm with custom trainControl: model
     model <- train(y ~., data, method = "glmnet", trControl = myControl,tuneGrid=tuneGrid)








     # Extract the best tuning parameters
     best_lambda <- model$bestTune$lambda
     best_alpha <- model$bestTune$alpha


     x <- as.matrix(data[,-which(names(data) == "y")])
     y <- data$y

     final_model <- glmnet::glmnet(x, y, alpha = best_alpha, lambda = best_lambda,family = "binomial")

     # Extract non-zero coefficients
     final_coefficients <- coef(final_model)
     non_zero_coefficients <- rownames(final_coefficients)[which(final_coefficients != 0)]

     # Create a formula with the non-zero coefficients, excluding the intercept
     formula <- as.formula(paste("y ~", paste(non_zero_coefficients[-1], collapse = " + ")))

     # Refit the model using glm to get p-values

     levels(data$y) <- c(0,1)

     glm_model <- glm(formula, data = data,family = 'binomial')

     # Fit the null model
     null_model <- glm(y ~ 1, data = data,family = 'binomial')

     # Perform likelihood ratio test to compare the final model with the null model
     lrt <- anova(null_model, glm_model, test = "Chisq")

     # Extract the p-value of the entire model
     model_p_value <- lrt$`Pr(>Chi)`[2]




     print(max(model$results$ROC))

     cat(paste0(typ,"\n"))

     report <-  rbind(report,data.frame(typ=typ,k=k,cell=cell,auc=ifelse(max(model$results$ROC)<0.5,0.5,max(model$results$ROC))))


     pvalReport <- rbind(pvalReport,data.frame(typ=typ,k=k,cell=cell,pval=model_p_value,auc=ifelse(max(model$results$ROC)<0.5,0.5,max(model$results$ROC)),sigVars=c(paste0(non_zero_coefficients,collapse = ","))))

     }


       if(typ=="random_module"){


         if(class(cellmoData)=="numeric"){

           print("Yes")

           next

         }

         ss <- sample(intersect(colnames(cell_data$cell_expression_data),noModule$V1),size=cellmoData$overalp)


         RandomcellmoData     <-      create_module_specific_dataset(module = ss,
                                                           expression_data = cell_data$cell_expression_data,
                                                           cell_type_metadata = cell_data$cell_metadata,
                                                           cell = cell)



         COLS <- colnames(RandomcellmoData$data)

         X    <- RandomcellmoData$data[,COLS!="Phenotype"]

         y    <- RandomcellmoData$data[,COLS=="Phenotype",drop=F]


         cat(paste0(typ,"\n\n"))

         data <-data.frame(y=factor(y$Phenotype,levels=c(0,1),labels = c("X0","X1")),X)
            myControl <- trainControl(
              method = "LOOCV",
           number =10,
           summaryFunction = customSummary,
           classProbs = TRUE # IMPORTANT!
           # verboseIter = TRUE
         )

         # After creating a custom trainControl object, it's easy to fit caret models that use AUC rather than accuracy to tune and evaluate the model. You can just pass your custom trainControl object to the train() function via the trControl argument.

         model <- train(y ~., data, method = "glmnet", trControl = myControl)


         report <-  rbind(report,data.frame(typ=typ,k=k,cell=cell,auc=ifelse(max(model$results$ROC)<0.5,1-max(model$results$ROC),max(model$results$ROC))))



       }


       if(typ=="dist_module"){


         if(class(cellmoData)=="numeric"){

           print("Yes")

           next

         }

         ss <- summary(noModule$V2)

         high_threshold    <- ss[5]
         medium_threshold  <- ss[3]
         # Categorize the vector elements into High and Low
         high <- noModule[noModule$V2 >= high_threshold,]
         medium <- noModule[noModule$V2 < high_threshold & noModule$V2 >= medium_threshold,]
         low <- noModule[ noModule$V2 < medium_threshold,]

         actual_mod_score <- genes_all[genes_all$V1 %in% modulesgene,]

           highgenesLen <-   sum(actual_mod_score$V2 > high_threshold )
           mediumgenesLen <- sum(actual_mod_score$V2 < high_threshold & actual_mod_score$V2 >= medium_threshold )
           lowLen <-  sum( actual_mod_score$V2 < medium_threshold )


           ss <- NULL
           if(highgenesLen!=0){ss1 <- sample(intersect(colnames(cell_data$cell_expression_data),high$V1),size=highgenesLen)
           ss <- c(ss,ss1)
           }
           if(mediumgenesLen!=0){ss2 <- sample(intersect(colnames(cell_data$cell_expression_data),medium$V1),size=mediumgenesLen)
           ss <- c(ss,ss2)
           }
           if(lowLen!=0){ss3 <- sample(intersect(colnames(cell_data$cell_expression_data),medium$V1),size=lowLen)
           ss <- c(ss,ss3)
           }





         #ss <- sample(intersect(colnames(cell_data$cell_expression_data),noModule$V1),size=cellmoData$overalp)




         DistcellmoData     <-      create_module_specific_dataset(module = ss,
                                                                     expression_data = cell_data$cell_expression_data,
                                                                     cell_type_metadata = cell_data$cell_metadata,
                                                                     cell = cell)



         COLS <- colnames(DistcellmoData$data)

         X    <- DistcellmoData$data[,COLS!="Phenotype"]

         y    <- DistcellmoData$data[,COLS=="Phenotype",drop=F]


         cat(paste0(typ,"\n\n"))

         data <-data.frame(y=factor(y$Phenotype,levels=c(0,1),labels = c("X0","X1")),X)
         myControl <- trainControl(
           method = "LOOCV",
           number =10,
           summaryFunction = customSummary,
           classProbs = TRUE # IMPORTANT!
           # verboseIter = TRUE
         )

         # After creating a custom trainControl object, it's easy to fit caret models that use AUC rather than accuracy to tune and evaluate the model. You can just pass your custom trainControl object to the train() function via the trControl argument.

         model <- train(y ~., data, method = "glmnet", trControl = myControl)


         report <-  rbind(report,data.frame(typ=typ,k=k,cell=cell,auc=ifelse(max(model$results$ROC)<0.5,1-max(model$results$ROC),max(model$results$ROC))))



       }



if(typ=="top_mod"){


         if(class(cellmoData)=="numeric"){

           print("Yes")

           next

         }


        topNoMod <- intersect(noModule$V1,ldak_genes_top$V1)
        ss <- sample(topNoMod,size=cellmoData$overalp)


         topcellmoData     <-      create_module_specific_dataset(module = ss,
                                                                     expression_data = cell_data$cell_expression_data,
                                                                     cell_type_metadata = cell_data$cell_metadata,
                                                                     cell = cell)



         COLS <- colnames(topcellmoData$data)

         X    <- topcellmoData$data[,COLS!="Phenotype"]

         y    <- topcellmoData$data[,COLS=="Phenotype",drop=F]


         cat(paste0(typ,"\n\n"))

         data <-data.frame(y=factor(y$Phenotype,levels=c(0,1),labels = c("X0","X1")),X)
         myControl <- trainControl(
           method = "LOOCV",
           number =10,
           summaryFunction = customSummary,
           classProbs = TRUE # IMPORTANT!
           # verboseIter = TRUE
         )

         # After creating a custom trainControl object, it's easy to fit caret models that use AUC rather than accuracy to tune and evaluate the model. You can just pass your custom trainControl object to the train() function via the trControl argument.

         model <- train(y ~., data, method = "glmnet", trControl = myControl)


         report <-  rbind(report,data.frame(typ=typ,k=k,cell=cell,auc=ifelse(max(model$results$ROC)<0.5,1-max(model$results$ROC),max(model$results$ROC))))



       }











       }

}


pval <- sum(report[report$cell==cell & report$typ=="random_module","auc"]>unique(report[report$cell==cell & report$typ=="real_module","auc"]) )/length(report[report$cell==cell& report$typ=="random_module","auc"])
pval_df <- data.frame(cell=cell,k=k,auc=unique(report[report$cell==cell & report$typ=="real_module","auc"]),pval=pval)



file_path <- paste0(output_directory, "pvalreport3.txt")
pvalExact_path <- paste0(output_directory, "pvalexact.txt")
file_exists <- file.exists(file_path)

# Write the data to the file
write.table(pval_df, file=file_path, append=TRUE, row.names=FALSE, col.names=!file_exists)

write.table(pvalReport,file=pvalExact_path,append=TRUE, row.names=FALSE, col.names=!file_exists)




if(class(cellmoData)!="numeric"){

  # Here I am saving all of the plots
  dir.create(paste0(output_directory,"/plots/",cell))

  file_name <- paste0(paste0(output_directory,"/plots/",cell),"/plot_k", k, sep = "")

  # Add p-values to your plot
  p <- ggplot(report, aes(x = typ, y = auc, fill = typ)) +
    geom_boxplot() +
    facet_wrap(~ cell, scales = "free_y", ncol = 2) +
    labs(x = "Module Type", y = "AUC", title = paste0("Comparison of different settings by Cell Type for Module",k),
         fill = "Type") +
    theme_minimal()

  ggsave(file = paste0(file_name,".pdf"), plot = p, width = 6, height = 4, dpi = 300)
  saveRDS(report,file=paste0(file_name,".rds"))





}else{


  next
}



}


  pvalReport$adjustedpval <- p.adjust(pvalReport$pval,method="BH")
  pvalReport <- pvalReport %>% distinct()
  write.table(pvalReport,file=paste0("/ix/djishnu/Javad/RA/Data_Analysis/results/July22Results/",cell_types[as.numeric(bashInp)],"-pvals.txt"))



}

# ## Show this for the power of the test
# n <- dim(X)[1]
# predictions <- c()
#
# for (i in 1:n) {
#   model <- glmnet(X[-i,], y$Phenotype[-i], family="binomial",alpha = 1,lambda = 0.01430420)
#   predictions <- c(predictions, predict(model, newx=X[i,]))
# }
#
# library(pROC)
#
# roc(y$Phenotype, predictions)






# # Load required libraries
# library(ggplot2)
# library(gridExtra)
#
# # Initialize an empty list to store plots
# plotList <- list()
# dataList <- list()
# res <- list.files("/ix/djishnu/Javad/RA/Data_Analysis/results/June18_results/",pattern=".rds")
#
# # Loop through .res
# for (k in 1:14) {
#   # Read data from RDS file
#   dataList[[k]] <- readRDS(paste0("/ix/djishnu/Javad/RA/Data_Analysis/results/June18_results/plot_k",k,".pdf.rds"))
#
# }







## Filtering Significant Modules



#
#
# # Define the desired order of levels for the typ variable
# order_levels <- c("real_module", "random_module", "dist_module")
#
# # Convert typ into a factor with the desired order of levels
# report$typ <- factor(report$typ, levels = order_levels)
#
# my_palette <- c("real_module" = "black", "dist_module" = "grey", "random_module" = "lightblue")
#
# p <- ggplot(report, aes(x = typ, y = auc, fill = typ)) +
#   geom_boxplot() +
#   facet_wrap(~ cell, scales = "free_y", ncol = 2) +
#   labs(x = NULL, y = "AUC", title = paste0("Moduel",k),
#        fill = "Type") + scale_fill_manual(values = my_palette) +  # Use custom color palette
#   theme(
#     panel.grid.major = element_blank(),
#     panel.grid.minor = element_blank(),
#     panel.background = element_rect(fill = "white"),
#     panel.border = element_rect(color = "black", fill = NA),
#     strip.text = element_text(color = "black", size = 14, face = "bold", hjust = 0.5),
#     strip.background = element_rect(fill = "lightgreen", color = "black"),
#     axis.title.x = element_blank(),
#     axis.text.x = element_blank(),
#     axis.ticks.x = element_blank()
#   )
#
# p
#
# plotList <- list()
# M <- 0
#
# for(k in c(1,3,6,9,10,12,13)){
#
# report <- dataList[[k]]
#
# M <- M+1
# report$typ <- fct_recode(report$typ,
#                          "same distribution" = "dist_module",
#                          "random" = "random_module",
#                          "actual" = "real_module")
#
# p <- ggplot(report, aes(x = typ, y = auc, fill = typ)) +
#   geom_boxplot() +
#   facet_wrap(~ cell, scales = "free_y", ncol = 2) +
#   labs(x = NULL, y = "AUC", title = paste0("Moduel",k),
#        fill = "Type") +
#   theme(axis.title.x=element_blank(),
#         axis.text.x=element_blank(),
#         axis.ticks.x=element_blank(),legend.position="none")+theme_bw()
#
# plotList[[M]] <- p
#
#
#
# }
#
#
#
#
# # # Combine plots into one figure
# # library(ggpubr)
# #
# # allplots <- grid.arrange(plotList[[1]],
# #                          plotList[[3]],
# #                          plotList[[5]],nrow=1)
# #
# #
# #
# #
# # ggsave(file = "/ix/djishnu/Javad/RA/PlotsAllfigs.pdf",allplots, width = 6, height = 4, dpi = 300)
# #
# # # Print the combined figure
# # print(allplots)
# #
#
#
#
#
#
# library(ggplot2)
#
# plotList <- list()
# M <- 0
#
#
# p <- ggplot(report, aes(x = auc, fill = typ)) +
#   geom_density(alpha = 0.7) +
#   theme_light() +
#   ylab("Density") +
#   xlab("AUC") +
#   facet_wrap(~ cell, scales = "free_y", ncol = 2) +
#   labs(title = paste0("Density Plot of AUC for module ",k))+
#   ggplot2::geom_vline(xintercept = realval, linetype = "longdash", colour = "red",size=2)
# #+theme(legend.position = "none")
#
# # Print the plot
# plotList[[1]] <- p
#
# for(k in c(3,6,9,10,12,13)){
#
# M <- M+1
# report <- dataList[[k]]
# realval <- report[report$typ=="real_module","auc"][1]
# report  <-report[report$typ!="real_module",]
#
#
#
#
# # Create the density plot
# p <- ggplot(report, aes(x = auc, fill = typ)) +
#   geom_density(alpha = 0.7) +
#   theme_light() +
#   ylab("Density") +
#   xlab("AUC") +
#   facet_wrap(~ cell, scales = "free_y", ncol = 2) +
#   labs(title = paste0("Density Plot of AUC for module ",k))+
#   ggplot2::geom_vline(xintercept = realval, linetype = "longdash", colour = "red",size=2)+theme(legend.position = "none")
#
# # Print the plot
# plotList[[M]] <- p
#
# }
#
#
#
#
#
#
#
#
# library(gridExtra)
# library(lemon)
# legend <- lemon::g_legend(p1 + theme(legend.position='bottom'))
#
#
#
#
# allplots <- grid.arrange(plotList[[1]],
#              plotList[[2]],
#              plotList[[3]],
#              plotList[[4]],
#              plotList[[5]],
#              plotList[[6]],ncol=2)+legend
#
#
# ggsave(allplots,filename = "/ix/djishnu/Javad/RA/allplots.pdf")
#
#
#
# ################################################################################
# #
# #                     Getting Pvalues for each cell type and module
# ################################################################################
#
# for (k in 1:14) {
#   # Read data from RDS file
#   dataList[[k]] <- readRDS(paste0("/ix/djishnu/Javad/RA/Data_Analysis/results/plot_k",k,".pdf.rds"))
#
#
# }
#
# cell_types          <- c('B cell', 'Fibro', 'Mono', 'T cell')
#
# report <- data.frame(k=numeric,cell=character,pval=numeric)
#
# for (k in 1:14) {
#
# dat <- dataList[[k]]
# for (c  in 1:length(cell_types)){
#
#
#   randomval <- dat %>% filter(typ=="random_module",cell==cell_types[c])
#   realval <- dat[dat$typ == "real_module" & dat$cell == cell_types[c], ]
#   actual <- unique(realval$auc)
#   pval <- 1-sum(actual>randomval$auc)/length(randomval$auc)
#   if(is.na(pval)){pval <- 1}
#
#
#   report <- rbind(report,data.frame(k=k,cell=cell_types[c],pval=pval))
#
#
# }}
#
# saveRDS(report,"Data_Analysis/results/ExpressionPval.rds")
# report$Sig01 <- ifelse(report$pval < 0.01,1,0)
# report$Sig0.1 <- ifelse(report$pval < 0.1,1,0)
# report$Sig0.25 <- ifelse(report$pval < 0.25,1,0)
# idx <- apply(report,1,function(x){any(x[4:6]!=0)})
# filtredReprot <- report[idx,]
# saveRDS(filtredReprot,"Data_Analysis/results/ExpressionPvalFinal.rds")
#
