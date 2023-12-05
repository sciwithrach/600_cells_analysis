# Quick per cell QC and filtering with mitochondrial DNA

#source('00_constants.R')

library(Seurat)
library(SingleCellExperiment)
library(scuttle)
library(stringr)

MitoAndQC <- function(srat) {
  
  # Load Seurat object
#  srat <- ReadObject(obj.name)

  #Convert to SingleCellExperiment object
  sce <- as.SingleCellExperiment(srat)

  #Run quick per cell QC
  sce <- scuttle::quickPerCellQC(sce)

  #Add mitoDNA proportion
  location <- rowRanges(sce)
  is.mito <- sapply(names(seqnames(location)), str_detect, '^mt*')
  sce <- addPerCellQCMetrics(sce, subsets = list(Mito = is.mito))
  qc_df <- perCellQCMetrics(sce, subsets = list(Mito = is.mito))

  #Plot mitoDNA proportion
  p <- plotColData(sce, x = 'sum', y = 'subsets_Mito_percent')
  print(p)

  #Prompt %mito to exclude
  mito_pc <- readline(
    prompt = 'Enter maximum % of mitochondrial DNA per cell to INCLUDE (e.g. 80):')
  mito_pc <- as.numeric(gsub("\\D", "", mito_pc))
  
  #Define filters
  sce_keep_mito <- sce$subsets_Mito_percent < mito_pc
  sce_keep_qc <- sce$discard == F
  print(sum(sce_keep_mito))
  print(sum(sce_keep_qc))
  
  #Filter object
  sce <- sce[, sce_keep_mito & sce_keep_qc]
  print(dim(sce))
  
  #Convert to Seurat object and save
  srat.filtered <- as.Seurat(sce, counts = 'counts', data = 'logcounts')
  return(srat.filtered)
  
  #SaveObject(srat.filtered, paste0(obj.name, '_filtered', sep=""))

}

#Run function!
#MitoAndQC(obj.name)