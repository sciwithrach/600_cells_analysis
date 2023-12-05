# Create a Seurat object from raw Parse data and save the output as an RDS

#Load libraries
library(Seurat)

#ReadParseBio with correct count_matrix name
ReadParseBio_count <- function(data.dir, ...) {
  mtx <- file.path(data.dir, "count_matrix.mtx")
  cells <- file.path(data.dir, "cell_metadata.csv")
  features <- file.path(data.dir, "all_genes.csv")
  return(ReadMtx(
    mtx = mtx,
    cells = cells,
    features = features,
    cell.column = 1,
    feature.column = 2,
    cell.sep = ",",
    feature.sep = ",",
    skip.cell = 1,
    skip.feature = 1,
    mtx.transpose = TRUE
  ))
}

#Define function
CreateObj <- function(data.dir, object.name) {
  
  parse <- ReadParseBio_count(here(data.dir, 
                          'all-sample', 
                          'DGE_filtered'))
  
  cell.meta <- read.csv(here(data.dir, 
                             'all-sample', 
                             'DGE_filtered', 
                             "cell_metadata.csv"), 
                              row.names = 1)
  
  srat <- CreateSeuratObject(parse, 
                     min.features = 100, 
                     min.cells = 100, 
                     names.field = 0, 
                     meta.data = cell.meta)
  
  SaveObject(srat, object.name)
  
  return(srat)
}

#Run function
#CreateObj(data.dir, object.name)
