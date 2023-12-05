#source('00_constants.R')

library(MAST)
library(dplyr)

MarkerAnalysisPlot <- function(srat, markers, final.res) {
  
#  #Read in object
#  srat <- object.name
#  assign(srat, 'srat', envir = sys.frame())
  
  resolution <- paste0('integrated_snn_res.', final.res, sep = "")
  
#  markers <- lapply(
#    levels(paste0(get(srat), '@meta.data$', resolution, sep = "")),
#    function(x) {
#      FindAllMarkers(object.name, 
#                     test.use = 'MAST', 
#                     group.by = resolution)
#    }
#  )
  
  #Get top 5 markers for mouse
  
  top_markers <- lapply(markers, function(x) top_n(x, n=5, wt=avg_log2FC))
  top_markers <- bind_rows(top_markers, .id = "id")
  
  #Plot top 5 markers for mouse
  # top_plot <- unique(rownames(top_markers))
  top_plot <- sub("\\.\\d+$", "", rownames(top_markers))
  plot <- DotPlot(srat, features = top_plot, group.by = resolution) + coord_flip()
  
  #Show and save plot
  return(plot)
  
  #print('Saving plot of top 5 markers per cluster...')
  #SaveFigure(plot, paste0(object.name, "_Top5MarkersPerCluster", sep = ""),
  #           width = 9, height = 20)
}
