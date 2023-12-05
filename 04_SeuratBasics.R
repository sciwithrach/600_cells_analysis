# Standard boring Seurat stuff

#source('00_constants.R')

library(Seurat)
library(ggplot2)
library(sctransform)
library(glmGamPoi)
library(patchwork)
library(dplyr)

SeuratBasics <- function(object.name, transform = NULL) {
  
  #Run all the Seurat stuff
#  srat <- ReadObject(object.name)
  
  # If transform = T, run SCTransform (normalise, scale, find HVFs)
  if (!is.null(transform)) {
    srat <- SCTransform(object.name, vst.flavor = "v2", verbose = FALSE)
    srat <- RunPCA(srat, npcs = 50, verbose=F)
  } else {
    srat <- RunPCA(object.name, npcs = 50, verbose=F)
  }
  
  p <- ElbowPlot(srat, ndims=50)
  
  #Show the elbow plot
  print(p)
  
  #Save the elbow plot
  #print('Saving Elbow Plot...')
  #SaveFigure(p, paste0(object.name, '_ElbowPlot', sep = ""), 
  #           width = 10, 
  #           height = 5)
  
  #Save Seurat object after PCA
  #print('Saving Seurat object after PCA...')
  #SaveObject(srat, paste0(object.name, '_after_PCA'))
  
  #Let user define number of Principal Components
  npcs <- readline('Enter number of PCs to use for UMAP (e.g. 30):')
  npcs <<- as.numeric(gsub("\\D", "", npcs)) # Note - this operator creates a global variable named 'npcs' to be used in the next step
  
  #Run UMAP using user-defined PCs
  srat <- RunUMAP(srat, reduction = 'pca', dims = 1:npcs)
  
  #Show UMAP
  p <- DimPlot(srat, reduction='umap')
  print(p)
  
  #Save UMAP
  #print('Saving UMAP...')
  #SaveFigure(p, paste0(object.name, '_UMAP', sep = ""),
  #           width = 10,
  #           height = 10)
  
  #Save Seurat object after UMAP
  #print('Saving Seurat object after UMAP...')
  #SaveObject(srat, paste0(object.name, '_after_UMAP'))
  
  return(srat)
}

#Run function!
#SeuratBasics(object.name)

