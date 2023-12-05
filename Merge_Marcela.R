### merge Zeisel 10x data with new dropseq data

## make meta.data columns with clusters obtained for each data set
seurat_OB[["X10x_clusters"]] <- seurat_OB$seurat_clusters
My_new_OB[["res06_clusters"]] <- My_new_OB$seurat_clusters

data.list <- list(seurat_OB, My_new_OB)

for (i in 1:length(data.list)) {
  data.list[[i]] <- NormalizeData(data.list[[i]], verbose = FALSE)
  data.list[[i]] <- FindVariableFeatures(data.list[[i]], selection.method = "vst", 
                                         nfeatures = 3000, verbose = FALSE)
}

### for some reason the function bit sometimes does not work, so I'm running it again for getting a list of 10000 variable features
My_new_OB <- FindVariableFeatures(My_new_OB, selection.method = "vst", nfeatures = 10000)
seurat_OB <- FindVariableFeatures(seurat_OB, selection.method = "vst", nfeatures = 10000)

var_genes_NewOB <- VariableFeatures(My_new_OB)
var_genes_10x <- VariableFeatures(seurat_OB)
var_genes <- intersect(var_genes_NewOB,var_genes_10x)

all.genes_NewOB <- rownames(My_new_OB)
all.genes_10x <- rownames(seurat_OB)
all.genes_int <- intersect(all.genes_NewOB,all.genes_10x)

########## TEST THIS... I'm looking into using more genes for the integration (_0).. the one I run only had 2000
test.anchors <- FindIntegrationAnchors(object.list = data.list, dims = 1:30)
test.anchors_0 <- FindIntegrationAnchors(object.list = data.list, dims = 1:30, 
                                         anchor.features = var_genes)
test.integrated <- IntegrateData(anchorset = test.anchors, dims = 1:30)
test.integrated_0 <- IntegrateData(anchorset = test.anchors_0, dims = 1:30, features.to.integrate = all.genes_int)

test.integrated <- ScaleData(test.integrated, assay = "RNA") ## note 11-05-20.. why did I run this?

DefaultAssay(test.integrated) <- "integrated"
DefaultAssay(test.integrated_0) <- "integrated"

# Run the standard workflow for visualization and clustering
test.integrated <- ScaleData(test.integrated, verbose = FALSE)
test.integrated_0 <- ScaleData(test.integrated_0, verbose = FALSE)
### do not run findVariableGenes again.. it's already some sort of intersection between the two datasets

test.integrated <- RunPCA(test.integrated, npcs = 40, verbose = FALSE)
test.integrated_0 <- RunPCA(test.integrated_0, npcs = 40, verbose = FALSE)

DimPlot(test.integrated, reduction = "pca")
FeaturePlot(test.integrated, "Eno2", reduction = "pca")

DimPlot(test.integrated_0, reduction = "pca")
FeaturePlot(test.integrated_0, "Eno2", reduction = "pca")

## trying to label batches
chips_10x <- as.character(seurat_OB$ChipID)
NewOB_Runs <- as.character(My_new_OB$orig_ident)
batches <- as.character(c(chips_10x, NewOB_Runs))
test.integrated$batches <- (batches)
DimPlot(test.integrated, reduction = "pca", group.by = "batches")
test.integrated_0$batches <- (batches)
DimPlot(test.integrated_0, reduction = "pca", group.by = "batches")

ElbowPlot(test.integrated, ndims = 40, reduction = "pca")
ElbowPlot(test.integrated_0, ndims = 40, reduction = "pca")


test.integrated <- RunUMAP(test.integrated, dims = 1:30)
DimPlot(test.integrated, reduction = "umap", group.by = "batches")
FeaturePlot(test.integrated, "Calb2", reduction = "umap")

test.integrated_0 <- RunUMAP(test.integrated_0, dims = 1:28)
DimPlot(test.integrated_0, reduction = "umap", group.by = "batches")
FeaturePlot(test.integrated_0, "Calb2", reduction = "umap")

test.integrated <- FindNeighbors(test.integrated, dims = 1:30)
test.integrated <- FindClusters(test.integrated, resolution = c(0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4))
table(Idents(test.integrated))
DimPlot(test.integrated, label = TRUE) + NoLegend()

test.integrated_0 <- FindNeighbors(test.integrated_0, dims = 1:28)
test.integrated_0 <- FindClusters(test.integrated_0, resolution = c(0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4))
table(Idents(test.integrated_0))
DimPlot(test.integrated_0, label = TRUE) + NoLegend()


## check cluster stability
head(test.integrated[[]])
clustree(test.integrated, prefix = "integrated_snn_res.", node_colour = "sc3_stability")

clustree(test.integrated_0, prefix = "integrated_snn_res.", node_colour = "sc3_stability")

### clusters seem more stable at resolution = 0.6
test.integrated <- FindClusters(test.integrated, resolution = 0.6)
table(Idents(test.integrated))
DimPlot(test.integrated, label = TRUE) + NoLegend()

test.integrated_0 <- FindClusters(test.integrated_0, resolution = 0.6)
table(Idents(test.integrated_0))
DimPlot(test.integrated_0, label = TRUE) + NoLegend()

### absolutely great!!!






### ploting gene expression... cannot use integrated assay for this
## so change the default assay to RNA that has the original values... should I scale them?
DefaultAssay(test.integrated_0) <- "RNA"
test.integrated_0 <- NormalizeData(test.integrated_0, normalization.method = "LogNormalize", scale.factor = 10000)
all.genes_int <- rownames(test.integrated_0)
test.integrated_0 <- ScaleData(test.integrated_0, features = all.genes_New_OB)

DotPlot(test.integrated_0, features = c("Ace2", "Tmprss2", "Ctsb", "Ctsl", "Dpp4", "Furin",
                                        "Anpep", "Tmprss11d", "St6gal1", "St3gal4", "Ceacam1"),
        scale.by = "size")

DoHeatmap(test.integrated_0, features = c("Ace2", "Tmprss2", "Ctsb", "Ctsl", "Dpp4", "Furin",
                                          "Anpep", "Tmprss11d", "St6gal1", "St3gal4", "Ceacam1"))


## need to label all the cell classes

seurat_clusters <- as.data.frame(test.integrated_0$seurat_clusters)
colnames(seurat_clusters) <- "Cell_class"
seurat_clusters$Cell_class <- as.character(seurat_clusters$Cell_class)

seurat_clusters$Cell_class[seurat_clusters$Cell_class == "0"] <- "Neurons"
seurat_clusters$Cell_class[seurat_clusters$Cell_class == "1"] <- "Neurons"
seurat_clusters$Cell_class[seurat_clusters$Cell_class == "2"] <- "Neurons"
seurat_clusters$Cell_class[seurat_clusters$Cell_class == "3"] <- "Neurons"
seurat_clusters$Cell_class[seurat_clusters$Cell_class == "4"] <- "Neurons"
seurat_clusters$Cell_class[seurat_clusters$Cell_class == "5"] <- "Astrocytes"
seurat_clusters$Cell_class[seurat_clusters$Cell_class == "6"] <- "OECs"
seurat_clusters$Cell_class[seurat_clusters$Cell_class == "7"] <- "Neurons"
seurat_clusters$Cell_class[seurat_clusters$Cell_class == "8"] <- "Neurons"
seurat_clusters$Cell_class[seurat_clusters$Cell_class == "9"] <- "Immune"
seurat_clusters$Cell_class[seurat_clusters$Cell_class == "10"] <- "Oligodendrocytes"
seurat_clusters$Cell_class[seurat_clusters$Cell_class == "11"] <- "Neurons"
seurat_clusters$Cell_class[seurat_clusters$Cell_class == "12"] <- "Neurons"
seurat_clusters$Cell_class[seurat_clusters$Cell_class == "13"] <- "Neurons"
seurat_clusters$Cell_class[seurat_clusters$Cell_class == "14"] <- "Astrocytes"
seurat_clusters$Cell_class[seurat_clusters$Cell_class == "15"] <- "Vascular"
seurat_clusters$Cell_class[seurat_clusters$Cell_class == "16"] <- "OPCs"
seurat_clusters$Cell_class[seurat_clusters$Cell_class == "17"] <- "Vascular"
seurat_clusters$Cell_class[seurat_clusters$Cell_class == "18"] <- "Neurons"
seurat_clusters$Cell_class[seurat_clusters$Cell_class == "19"] <- "Neurons"
seurat_clusters$Cell_class[seurat_clusters$Cell_class == "20"] <- "Immune"
seurat_clusters$Cell_class[seurat_clusters$Cell_class == "21"] <- "Neurons"
seurat_clusters$Cell_class[seurat_clusters$Cell_class == "22"] <- "Vascular"
seurat_clusters$Cell_class[seurat_clusters$Cell_class == "23"] <- "IPCs"
seurat_clusters$Cell_class[seurat_clusters$Cell_class == "24"] <- "Neurons"
seurat_clusters$Cell_class[seurat_clusters$Cell_class == "25"] <- "Blood"

test.integrated_0$Cell_class <- seurat_clusters$Cell_class

### for later... do Cell_type

final.integrated <- test.integrated_0[,!test.integrated_0$Cell_class == "Blood"]
### saved workspace with final.integragred seurat object alone in final_integrated_data.RData
## save as object
saveRDS(final.integrated, file = "integrated.rds")


final.integrated <- readRDS(file = "integrated.rds")


## plots for figure 1
DimPlot(final.integrated, reduction = "umap", group.by = "Cell_class")


DotPlot(final.integrated, 
        features = c( "Ceacam1","St3gal4","St6gal1","Tmprss11d",  
                      "Anpep","Furin","Dpp4","Hspa5","Bsg","Ctsl", "Ctsb","Tmprss2","Ace2",
                      "Foxc1", "Olig2", "Mog", "Frzb", "Snap25", "Top2a", "Itgam", "Mlc1"),
        scale.by = "size", group.by = "Cell_class") + RotatedAxis()


DotPlot(final.integrated, 
        features = c( "Ceacam1","St3gal4","St6gal1","Tmprss11d",  
                      "Anpep","Furin","Dpp4","Hspa5","Bsg","Ctsl", "Ctsb","Tmprss2","Ace2",
                      "Dcn", "Kcnj8","Foxc1",  "Itgam", "Mrc1", "Frzb","Mlc1", 
                      "Olig2", "Mog", "Top2a","Dcx", "Vip", "Th", "Calb2","Htr2c", "Eomes","Gad1","Snap25"),
        scale.by = "size") + RotatedAxis()


DotPlot(final.integrated, 
        features = c( "Dcn", "Kcnj8","Foxc1",  "Itgam", "Mrc1", "Frzb","Mlc1", 
                      "Olig2", "Mog", "Top2a","Dcx", "Vip", "Th", "Calb2","Htr2c", "Eomes","Gad1","Snap25",
                      "Ceacam1","St3gal4","St6gal1","Tmprss11d",  
                      "Anpep","Furin","Dpp4","Hspa5","Bsg","Ctsl", "Ctsb","Tmprss2","Ace2"),
        scale.by = "size", dot.min = 0.01) + RotatedAxis()



VlnPlot(final.integrated, 
        features = c("Foxc1", "Olig2", "Mog", "Frzb", "Snap25", "Top2a", "Itgam", "Mlc1"),
        group.by = "Cell_class")


DotPlot(final.integrated, 
        features = c("Ace2", "Bsg","Tmprss2", "Ctsb", "Ctsl", "Dpp4", "Furin",
                     "Anpep", "Tmprss11d", "St6gal1", "St3gal4", "Ceacam1", "Hspa5"),
        scale.by = "size", group.by = "Cell_class")

FeaturePlot(final.integrated, features = "Kcnj8", label = FALSE)
FeaturePlot(final.integrated, features = "Ace2", label = FALSE)
FeaturePlot(final.integrated, features = "Tmprss2", label = FALSE)
## fos scale colour viridis
FeaturePlot(final.integrated, features = "Th", label = FALSE) + scale_color_viridis(direction = -1)

FeaturePlot(final.integrated, features = "Kcnj8", label = FALSE) + NoLegend() + labs(title = NULL)


VlnPlot(final.integrated[,final.integrated$seurat_clusters == "11"], 
        features = c( "Ceacam1","St3gal4","St6gal1","Tmprss11d",  
                      "Anpep","Furin","Dpp4","Hspa5","Bsg","Ctsl", "Ctsb","Tmprss2","Ace2",
                      "Th", "Gad1", "Snap25"))

DotPlot(final.integrated[,final.integrated$seurat_clusters == c("11", "24")], 
        features = c( "Ceacam1","St3gal4","St6gal1","Tmprss11d",  
                      "Anpep","Furin","Dpp4","Hspa5","Bsg","Ctsl", "Ctsb","Tmprss2","Ace2",
                      "Th", "Gad1", "Snap25"), scale.by = "size") + RotatedAxis()


# make an sce with normalised counts to get the pct cells expressing covid genes per class and cluster
OB_counts <- final.integrated@assays$RNA@counts
OB_normcounts <- final.integrated@assays$RNA@data
sce_OB <- SingleCellExperiment(assays = list(counts = OB_counts, normcounts = OB_normcounts))
sce_OB$Cell_class <- final.integrated$Cell_class
sce_OB$seurat_clusters <- final.integrated$seurat_clusters

sce_OB <- calculateQCMetrics(sce_OB)

sce_OB_Astrocytes <- sce_OB[,sce_OB$Cell_class == "Astrocytes"]
sce_OB_Astrocytes <- calculateQCMetrics(sce_OB_Astrocytes)
sce_OB_Immune <- sce_OB[,sce_OB$Cell_class == "Immune"]
sce_OB_Immune <- calculateQCMetrics(sce_OB_Immune)
sce_OB_IPCs <- sce_OB[,sce_OB$Cell_class == "IPCs"]
sce_OB_IPCs <- calculateQCMetrics(sce_OB_IPCs)
sce_OB_Neurons <- sce_OB[,sce_OB$Cell_class == "Neurons"]
sce_OB_Neurons <- calculateQCMetrics(sce_OB_Neurons)
sce_OB_OECs <- sce_OB[,sce_OB$Cell_class == "OECs"]
sce_OB_OECs <- calculateQCMetrics(sce_OB_OECs)
sce_OB_Oligos <- sce_OB[,sce_OB$Cell_class == "Oligodendrocytes"]
sce_OB_Oligos <- calculateQCMetrics(sce_OB_Oligos)
sce_OB_OPCs <- sce_OB[,sce_OB$Cell_class == "OPCs"]
sce_OB_OPCs <- calculateQCMetrics(sce_OB_OPCs)
sce_OB_Vascular <- sce_OB[,sce_OB$Cell_class == "Vascular"]
sce_OB_Vascular <- calculateQCMetrics(sce_OB_Vascular)

sce_OB_pericytes <- sce_OB[,sce_OB$seurat_clusters == "17"]
sce_OB_pericytes <- calculateQCMetrics(sce_OB_pericytes)

sce_OB_Vasc_nonPericytes <- sce_OB_Vascular[,!sce_OB_Vascular$seurat_clusters == "17"]
sce_OB_Vasc_nonPericytes <- calculateQCMetrics(sce_OB_Vasc_nonPericytes)


pct_Ace2_Astrocytes <- rowData(sce_OB_Astrocytes)["Ace2",]$pct_dropout_by_counts
pct_Ace2_Astrocytes <- 100-(pct_Ace2_Astrocytes)
pct_Ace2_immune <- rowData(sce_OB_Immune)["Ace2",]$pct_dropout_by_counts
pct_Ace2_immune <- 100-(pct_Ace2_immune)
pct_Ace2_IPCs <- rowData(sce_OB_IPCs)["Ace2",]$pct_dropout_by_counts
pct_Ace2_IPCs <- 100 - (pct_Ace2_IPCs)
pct_Ace2_Neurons <- rowData(sce_OB_Neurons)["Ace2",]$pct_dropout_by_counts
pct_Ace2_Neurons <- 100 - (pct_Ace2_Neurons)
pct_Ace2_OECs <- rowData(sce_OB_OECs)["Ace2",]$pct_dropout_by_counts
pct_Ace2_OECs <- 100 - (pct_Ace2_OECs)
pct_Ace2_Oligos <- rowData(sce_OB_Oligos)["Ace2",]$pct_dropout_by_counts
pct_Ace2_Oligos <- 100 - (pct_Ace2_Oligos)
pct_Ace2_OPCs <- rowData(sce_OB_OPCs)["Ace2",]$pct_dropout_by_counts
pct_Ace2_OPCs <- 100 - (pct_Ace2_OPCs)
pct_Ace2_Vascular <- rowData(sce_OB_Vascular)["Ace2",]$pct_dropout_by_counts
pct_Ace2_Vascular <- 100 - (pct_Ace2_Vascular)

pct_Ace2_pericytes <- rowData(sce_OB_pericytes)["Ace2",]$pct_dropout_by_counts
pct_Ace2_pericytes <- 100 - (pct_Ace2_pericytes)

pct_Ace2_NONpericytes <- rowData(sce_OB_Vasc_nonPericytes)["Ace2",]$pct_dropout_by_counts
pct_Ace2_NONpericytes <- 100 - (pct_Ace2_NONpericytes)

pct_ACe2 <- as.numeric(cbind(pct_Ace2_Astrocytes, pct_Ace2_immune, pct_Ace2_IPCs, pct_Ace2_Neurons, pct_Ace2_OECs,
                             pct_Ace2_Oligos, pct_Ace2_OPCs, pct_Ace2_Vascular, pct_Ace2_pericytes, 
                             pct_Ace2_NONpericytes))

barplot(pct_ACe2, xlab = c("Astrocytes", "immune", "IPCs", "Neurons", "OECs",
                           "Oligos", "OPCs", "Vascular", "pericytes", 
                           "NONpericytes"))




### make an sce with DA neurons
DA_counts <- final.integrated@assays$RNA@counts[,final.integrated$seurat_clusters == "11"]

sce_DAcells <- SingleCellExperiment(assays = list(counts = DA_counts))
plotExpression(sce_DAcells,  features = c( "Ceacam1","St3gal4","St6gal1","Tmprss11d",  
                                           "Anpep","Furin","Dpp4","Hspa5","Bsg","Ctsl", "Ctsb","Tmprss2","Ace2",
                                           "Th", "Gad1", "Snap25"), exprs_values = "counts", log2_values = TRUE)



### label Zeisel and Datta datasets
## use batches column
batch_labels <- as.data.frame(final.integrated$batches)
colnames(batch_labels) <- "dataset"
batch_labels$dataset <- as.character(batch_labels$dataset)

batch_labels$dataset[batch_labels$dataset == "10X04"] <- "Zeisel"
batch_labels$dataset[batch_labels$dataset == "10X28"] <- "Zeisel"
batch_labels$dataset[batch_labels$dataset == "10X49"] <- "Zeisel"
batch_labels$dataset[batch_labels$dataset == "121316run1"] <- "New"
batch_labels$dataset[batch_labels$dataset == "121316run2"] <- "New"
batch_labels$dataset[batch_labels$dataset == "121516run1"] <- "New"
batch_labels$dataset[batch_labels$dataset == "121516run2"] <- "New"

final.integrated$dataset <- batch_labels$dataset

DimPlot(final.integrated[,final.integrated$dataset == "Zeisel"], reduction = "umap",
        group.by = "Cell_class") 


### the scale data has negative values... maybe I messed up the normalisation/scaling.. run again#############
DefaultAssay(final.integrated) <- "RNA"
final.integrated <- NormalizeData(final.integrated, normalization.method = "LogNormalize", scale.factor = 10000)
all.genes_final <- rownames(final.integrated)
final.integrated <- ScaleData(final.integrated, features = all.genes_final)
### it's ok... the scale data can have negative values


### to rule out that the negative values on dot plots are ok #############################
### I'll make a new seurat object with just the counts... so the scale.data is not available, to see if negative 
### valuest still show up on the dotPlot
counts_int <- final.integrated@assays$RNA@counts
## Build Seural object
Seurat_counts <- CreateSeuratObject(counts = counts_int, project = "GoodCells", min.cells = 0, 
                                    min.features = 0)

Seurat_counts[["Cell_class"]] <- final.integrated$Cell_class

DotPlot(Seurat_counts, 
        features = c( "Hspa5","Ceacam1","St3gal4","St6gal1","Tmprss11d",  
                      "Anpep","Furin","Dpp4","Ctsl", "Ctsb","Tmprss2","Bsg","Ace2",
                      "Foxc1", "Olig2", "Mog", "Frzb", "Snap25", "Top2a", "Itgam", "Mlc1"),
        scale.by = "size", group.by = "Cell_class") + RotatedAxis()
### negative values are still there... they even show up on the plots they have on the vignettes
## I think that DotPlot takes the raw counts and scales them
### see https://github.com/satijalab/seurat/issues/783
## The values in DotPlot are extracted from the @data slot, averaged, and then passed to scale. 
## These are then Min-Maxed based on the col.min and col.max parameter values.

## from the vignette
## Shifts the expression of each gene, so that the mean expression across cells is 0
## Scales the expression of each gene, so that the variance across cells is 1
## This step gives equal weight in downstream analyses, so that highly-expressed genes do not dominate
####################################

