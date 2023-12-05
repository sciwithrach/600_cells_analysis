library(Seurat)
library(SingleCellExperiment)
library(clustree)
library(readr)

#source('00_constants.R')

ClustreeLoop <- function(srat) { # Loop adding resolution to object and plotting clustree
  
  answer <- 'n'
  
  while (answer == 'n') {
    # Take user input for additional resolutions
    test.res <- c(
      readline('Enter a resolution to plot e.g. 0.01:')
    )
    
    test.res <- parse_number(test.res)
    
    # Plot new resolutions
    for (i in 1:length(test.res)) {
      srat <- FindClusters(srat, resolution = test.res[i])
    }
    
    plot <- clustree(srat)
    
    print(plot)
    
    # Accept?
    answer <- readline('Enter y to accept and continue OR enter n to add new resolutions:')
  }

}

Clustering <- function(object.name) {
  
  #Read in object
  #srat <- ReadObject(object.name)
  
  #Find neighbours
  srat <- Seurat:::FindNeighbors(object.name, dims=1:npcs)
  
  #Cluster stability analysis
  srat <- FindClusters(srat, resolution = 0.01)
  srat <- FindClusters(srat, resolution = 0.1)
  srat <- FindClusters(srat, resolution = 0.3)

  #Plot clustree
  plot <- clustree(srat)
  print(plot)
  
  #Run loop to add resolutions
  ClustreeLoop(srat)
  
  #Save final resolution for other functions
  final.res <- as.numeric(readline('Enter chosen clustering resolution (e.g. 0.3):'))
  #final.res <<- parse_number(final.res)
  
  #Save clustree
  #print('Saving Clustree Plot...')
  #SaveFigure(plot, paste0(object.name, "_Clustree"), width = 5, height = 8)
  
  #Run Louvain clustering
  srat <- FindClusters(srat, resolution = final.res)
  
  #Save Seurat object
  #print('Saving Seurat Object after clustering...')
  #SaveObject(srat, paste0(object.name, '_AfterClustering'))
  
  return(srat)
}