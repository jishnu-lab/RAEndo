##################################################################
##                             Goal                             ##
##################################################################

# The main goal is to check performance of signficant IFNA genes
# compared to random using priyamvad's reactome list of genes on May 8
# This code is based on feature permutation

##################################################################
##                         Reading data                         ##
##################################################################


rm(list = ls()) # Clear existing workspace to ensure a fresh environment
cat("\014")  # Clears the console output
setwd("/ix/djishnu/Javad/RA/")
source("/ix/djishnu/Javad/RA/R/LOOV.R")
# Load necessary libraries
library(tidyverse)
library(biomaRt)


# Source required functions for analysis
source("Data_Analysis/Comparison_scripts/Funcs.R")
source("R/ZeroFiltering.R")


# Reading data files
expression_mt <- read.table('./Data/low_input_gene_sample_tpm_matrix.725714.tsv')
metadata_df <- read.table('./Data/metadata_for_bulk_RNA_seq.tsv')
#ensembltab_filtered <- read.csv("Data/EnsemblGene.map")

hsapiens_genes_all <- getBM(attributes = c("ensembl_gene_id",
                                           "hgnc_symbol", "gene_biotype"),
                            mart = useMart("ensembl", dataset = "hsapiens_gene_ensembl"))


toRemove <- c("RNA", "pseudogene", "artifact", "ribozyme") #pattern to match
position_to_remove <- grep(paste(toRemove,collapse="|"),
                           unique(hsapiens_genes_all$gene_biotype), value=TRUE)

hsapiens_genes <- subset(hsapiens_genes_all, subset = !(hsapiens_genes_all$gene_biotype %in% position_to_remove))
expression_matrix <- subset(expression_mt, subset = (rownames(expression_mt) %in% hsapiens_genes$ensembl_gene_id))
celltypeFilter <- metadata_df$type=="Fibro"
filtered_expression_matrix <- expression_matrix[, which(celltypeFilter)]
y <- metadata_df$Phenotype[celltypeFilter]

cell_expression_data <- t(filtered_expression_matrix)
mapping <- setNames(hsapiens_genes_all$hgnc_symbol, hsapiens_genes_all$ensembl_gene_id)
colnames(cell_expression_data) <- mapping[colnames(cell_expression_data)]
IFNApway <- read.table("/ix/djishnu/Javad/RA/Data/TNF_V3.txt",sep=",")
cell_types    <- c('Fibro', 'Mono', 'T cell') # Define cell types to analyze

sVarList <- list()
cellList <- list()
resList <- list()
zeroThre      <- 0.9

for(cell in cell_types[1]){

  ## Reading the data
  cat(paste0("Processing ",cell," ... \n"))

  # Read X
  cell_expression_data <- read.csv(file=paste0("./Data/",cell,"Exp.csv"),row.names = 1)
  column_thresh        <- ceiling(nrow(cell_expression_data)*(1-zeroThre))
  cell_expression_data <- ZeroFiltering(cell_expression_data,row_thresh =28 ,column_thresh =column_thresh )

  # Two ways
  ## Extracting TNF from the gene expression itself
  selectedGenes  <- union(colnames(cell_expression_data)[grep("^TNF",colnames(cell_expression_data))],IFNApway$V1)
  ii <- intersect(colnames(cell_expression_data),selectedGenes)

  ## Using reactome pathway
  #ii <- intersect(colnames(cell_expression_data),TNF$V1)
  sdzero<- which(apply(cell_expression_data[,ii],2,sd)==0)

  x <- scale(cell_expression_data[,ii],T,T)
  if(length(sdzero)!=0) {x <- x[,-sdzero]}
  # Read y

  #y <- read.csv(file=paste0("Data/y_",cell,".csv"),row.names = 1)

  train_data   <- data.frame(x,y=factor(y))
  levels(train_data$y) <- c("X1","X0")

  ## Doing Lasso


  lambda_value <- seq(0.0,0.12,0.01)
  resList[[cell]] <- LOOV(x,y,lambdas=lambda_value,family="binomial")
  fit  <- glmnet(x=x,y=y,alpha=1,lambda = names(resList[[cell]]$bestLambda))






  # Saving the significant variables into sVarList

  sVarList[[cell]] <- resList[[cell]]$sigVars[[resList[[cell]]$bestLambda]]
  RandPer <- list()

  ## For niter replicate run cross validation
  niter <- 1000
  similarityList  <-list()

  for(rep in 1:niter){
    # Generate the random Xs
    xR <- x[,sample(size=length(sVarList[[cell]]),1:ncol(x)),drop=F]
    resR  <- LOOV(x[,sample(size=length(sVarList[[cell]]),1:ncol(x))],y = y,lambdas = 0,family = "binomial")

    ## Permutation on sample space

    # Now you can call the function with a custom similarity threshold
    original_vector <- y
    # permuted_vector <- guided_permutation(original_vector,chnumber=2)
    # similarityList[[rep]] <- calculate_similarity(original_vector, permuted_vector)
    # resR  <- LOOV(x[,sVarList[[cell]]],y = permuted_vector,lambdas = 0)
    RandPer[[rep]]<- resR$perList

  }



  # run glmnet
  plotDf <- rbind(data.frame(perf=unlist(resList[[cell]]$BestPer),type="actual"),data.frame(perf=unlist(RandPer),type="Random"))
  cellList[[cell]]$pvalue<-sum(unlist(RandPer)>unlist(resList[[cell]]$BestPer))/(niter)
  cellList[[cell]]$plot <- ggplot(data=plotDf)+aes(x=type,y=perf,fill=type)+geom_boxplot()+scale_fill_manual(values=c("black","gray"))+theme_minimal() +  # Use a minimal theme that has less visual noise
    theme(panel.grid.major = element_blank(),  # Remove major grid lines
          panel.grid.minor = element_blank(),  # Remove minor grid lines
          panel.background = element_rect(fill="white", colour="white"),
          legend.position = "none") +ylab("auc")+xlab("type") # Set background to white





}

saveRDS(cellList,file=paste0("figures/IFNA_TNF_FromExpRes_andPath_FeatPerm_",zeroThre,".rds"))
saveRDS(resList,file = paste0("Data_Analysis/results/IFNA_TNFresultsFromExpRes_andPath_FeatPerm",zeroThre,".rds"))
saveRDS(RandPer,file=paste0("figures/IFNA_TNFandFromExpRes_andPath_FeatPerm",zeroThre,".rds"))
ggsave(cellList[[cell]]$plot,file=paste0("figures/IFNA_TNF_FromExpRes_andPath_FeatPerm_",zeroThre,".pdf"),height = 10,width = 10)
