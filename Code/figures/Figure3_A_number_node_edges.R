library(qgraph)

#################################################################
##               Making only sig modules network               ##
#################################################################
rm(list = ls())
library(viridis)
setwd("/ix/djishnu/Javad/RA/")

PPISig <- read.table("Data/significant_modules.txt",header = T)
ref <- read.table("Data/HomoSapiens_binary_co_complex_union.txt")
ref <- ref[!((ref$V1%in%PPISig$Source)&(ref$V1%in%PPISig$Target)),]
refS <- ref[sample(size=2000,dim(ref)[1]),]

colnames(refS) <- colnames(PPISig)
allPPI <- rbind(PPISig,refS)
colnames(allPPI) <- colnames(PPISig)

Modules <- readRDS("Data/significant_modules_members.RDS")
Modules[[15]] <- unique(c(refS$Source,refS$Target))


#################################################################
##                     Decode it to numbers                     ##
#################################################################
# The qgplot doesn't perfrom well with groups when the node are in character

unique_genes <- unlist(unique(c(allPPI$Source,allPPI$Target)))
unique_genes <- unique(unlist(Modules))

gene_mapping <- data.frame(gene = unique_genes, number = seq_along(unique_genes))

## Making groups
grps <- lapply(Modules, function(module) {
  # Mapping genes to numbers using gene_mapping
  numbers <- gene_mapping$number[match(module, gene_mapping$gene)]
  return(numbers)
})



## Delete this
PPISig_num <- cbind(gene_mapping$number[match(PPISig$Source,gene_mapping$gene)],
                    gene_mapping$number[match(PPISig$Target,gene_mapping$gene)])



allPPI <- cbind(gene_mapping$number[match(allPPI$Source,gene_mapping$gene)],
                gene_mapping$number[match(allPPI$Target,gene_mapping$gene)])


adj_module1 <- allPPI[(allPPI[,1] %in% grps[[15]]) & (allPPI[,2] %in% grps[[15]]) ,]


my_palette <- c(rainbow(14),"#FFFFFF")

qgraph(allPPI,
       directed = FALSE,
       layout="spring",
       groups = grps,
       vsize=1,
       labels =FALSE,
       repulsion=1,
       color=my_palette)


legend("topright",legend = c(1:15), fill = c(rainbow(14),"#FFFFFF"), title = "Modules")


#################################################

report <- data.frame(k=numeric(),number_nodes = numeric(),
                     number_edges = numeric())

for(i in c(seq_along(Modules))){

  number_nodes  <- length(Modules[[i]])
  number_edges <- nrow(allPPI[(allPPI[,1] %in% grps[[i]]) & (allPPI[,2] %in% grps[[i]]) ,])
  k <- i
  report <- rbind(report,data.frame(k=i,number_nodes=number_nodes,number_edges=number_edges))
}



