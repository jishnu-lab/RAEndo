## Making barplot for each module of the cell type
library(dplyr)
setwd("/ix/djishnu/Javad/RA/")
cell_types =c("B cell","Mono","T cell","Fibro")
color_pallete <- c("pink","blue","green","orange")
plotData <- list()
plotList <- list()
library(orca)
library(plotly)


for(f in 1:length(cell_types)){

plotData[[f]] <- read.table(paste0("Data_Analysis/results/July22Results/",cell_types[f],"-pvals.txt"))

annotations   <- list()

# Loop through the data to create annotations for significant bars
for(i in 1:nrow(plotData[[f]])) {
  if(as.numeric(plotData[[f]]$auc[i]) > 0.69 && plotData[[f]]$adjustedpval[i] < 0.05) {
    annotations[[length(annotations) + 1]] <- list(
      x = plotData[[f]]$k[i],
      y = as.numeric(plotData[[f]]$auc[i])+0.03,
      text = "*",
      showarrow = FALSE,
      font = list(size = 30,color = 'red')
    )
  }
}



# Create the bar plot with a black border
p <-  plotList[[f]]<- plot_ly(plotData[[f]], x = ~k, y = ~as.numeric(auc), type = "bar",
        marker = list(
          color = color_pallete[f],  # Fill color for the bars
          line = list(
            color = 'black',  # Border color
            width = 3  # Border width
          )
        )) %>%
  layout(title = list(
    text = cell_types[f],  # Set the title text
    font = list(size = 24),  # Set the title font size
    x = 0.5,  # Center the title horizontally
    y = 0.9  # Lower the title vertically
  ),
  yaxis = list(
    range = c(0.49, 0.85),
    showgrid = FALSE,
    tickfont = list(size = 18,weight = "bold"),
    titlefont = list(size = 40,family = "Times New Roman",weight = "bold"),
    title = "AUC",
    zerolinecolor = 'black',
    zerolinewidth = 2,
    linecolor = 'black',
    linewidth = 5
  ),
  xaxis = list(
    tickmode = "array",
    tickvals = plotData[[1]]$k,
    ticktext = plotData[[1]]$k,
    zerolinecolor = 'black',
    zerolinewidth = 2,
    linecolor = 'black',
    linewidth = 5,
    showgrid = FALSE,
    tickfont = list(size = 18,weight = "bold"),
    titlefont = list(size = 40,family = "Times New Roman",weight = "bold"),
    title = "Modules"
  ),
  plot_bgcolor = "rgba(0,0,0,0)",
  paper_bgcolor = "rgba(0,0,0,0)",
  annotations = annotations
  )

save_image(p, paste0("/ix/djishnu/Javad/RA/figures/figure4A/plot_", cell_types[f], ".pdf"))
}









