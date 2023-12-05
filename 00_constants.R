# Constants for scRNAseq pipeline
# IMPORTANT:
# Change the names of the top-level directory (line 13) 
# and path for outputs (line 19) before running

#set up renv environment
renv::restore()

#set seed for reproducibility
set.seed(444)

#Set top-level directory using the here package
here::i_am("pipeline/00_constants.R") #CHANGE DIRECTORY BEFORE RUNNING
library(here)
return(here())
setwd(here())

#Set path for outputs
out_dir <- 'poster_analysis' # CHANGE NAME BEFORE RUNNING

ifelse(!dir.exists(here(out_dir)), 
       dir.create(out_dir),
       "Folder exists already")

#Set paths to save and read Seurat objects
obj_path <- here(out_dir, 'r_objects')
ifelse(!dir.exists(obj_path), 
       dir.create(obj_path), 
       "Folder exists already")

#Make and/or set path to save figures
fig_path <- here(out_dir, 'r_figures')
ifelse(!dir.exists(fig_path), 
       dir.create(fig_path), 
       "Folder exists already")
