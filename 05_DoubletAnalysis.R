library(Seurat)
library(SingleCellExperiment)
library(scDblFinder)

#source('00_constants.R')

DoubletRemoval <- function(object.name, npcs) {
  
  #Read in objects
#  srat <- ReadObject(object.name)
  sce <- as.SingleCellExperiment(object.name)
  
  #Run doublet analysis
  sce <- scDblFinder(sce, dims = npcs)
  
  #Plot nice scDoubletFinder tutorial plot
  p <- plot_grid(
    plotUMAP(sce, 
             colour_by = "scDblFinder.score"), 
    plotUMAP(sce,
             colour_by = "scDblFinder.class"), 
    plotUMAP(sce, 
             colour_by = "sample"), ncol = 3)
  
  #Show and save plot
  print(p)
 # print('Saving Doublet Analysis Plots...')
 # SaveFigure(p, paste0(object.name, 'DoubletAnalysisPlot', sep = ""))
  
  #Filter SCE object and convert back to Seurat
  sce <- sce[, sce$scDblFinder.class == 'singlet']
  srat_filtered <- as.Seurat(sce)
  
  
  #Save Seurat object
#  print('Saving Seurat Object after doublet removal...')
#  SaveObject(srat_filtered, paste0(object.name, '_AfterDoubletRemoval'))
  
  return(srat_filtered)
}

#Run function!
#DoubletRemoval(object.name, npcs)

