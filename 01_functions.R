# Convenience functions for Parse pipeline

#Save figures
SaveFigure <- function(plots, name, type = "png", width, height, res){
  if(type == "png") {
    png(paste0(fig_path, '/', name, ".", type),
        width = width, height = height, units = "in", res = 200,
        bg = NA)
  } else {
    pdf(paste0(fig_path, '/', name, ".", type),
        width = width, height = height,
        bg = NA)
  }
  print(plots)
  dev.off()
}

#Save Seurat objects
SaveObject <- function(object, name){
  saveRDS(object, here(out_dir, 'r_objects', paste0(name, ".RDS")))
}

#Read Seurat objects
ReadObject <- function(name){
  readRDS(here(out_dir, 'r_objects', paste0(name, ".RDS")))
}