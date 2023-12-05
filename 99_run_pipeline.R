# SOURCE CONSTANTS AND FUNCTIONS ###########################

source('00_constants.R')
source('01_functions.R')
# source all things

# CREATE SEURAT OBJECT #####################################

source('02_createobject.R')

data.dir <- # Data directory containing cells, features, and metadata from Parse
object.name <- # Name for Seurat object to be saved as RDS

CreateObj(data.dir, object.name)

# RUN PER CELL QC ##########################################

source('03_percellQC.R')

MitoAndQC(object.name) 

# RUN SEURAT BASICS ########################################

source('04_SeuratBasics.R')

SeuratBasics(paste0(object.name, '_filtered', sep=""))

# DOUBLET ANALYSIS #########################################

source('05_DoubletAnalysis.R')

DoubletRemoval(paste0(object.name, '_filtered_AfterUMAP'), npcs = npcs)

# SEURAT BASICS AFTER DOUBLET REMOVAL ######################

SeuratBasics(paste0(object.name, '_filtered_AfterUMAP_AfterDoubletRemoval'))

# CLUSTERING ###############################################

source('06_Clustering.R')

Clustering(paste0(object.name, '_filtered_AfterUMAP_AfterDoubletRemoval_filtered'))

# MARKER ANALYSIS ##########################################

source('07_MarkerAnalysis.R')

MarkerAnalysis(object.name, final.res)

# BIOMARKER ANALYSIS #######################################

BioMarkerAnalysis(object.name)

# GO ANALYSIS ##############################################

GOAnalysis(object.name)

# ORTHOLOGUES ##############################################
library(orthogene)