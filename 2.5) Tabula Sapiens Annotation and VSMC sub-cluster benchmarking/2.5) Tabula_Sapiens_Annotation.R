
install.packages("Seurat")


# Packages # load libraries
library(Seurat)
library(SeuratDisk)
library(tidyverse)

setwd('C:/Users/Home/Documents/Dissertation/Carotid_Annotation/')



Convert("TS_Vasculature.h5ad", ".h5seurat")
# This creates a copy of this .h5ad object reformatted into .h5seurat inside the example_dir directory

# This .d5seurat object can then be read in manually

obj <- LoadH5Seurat(".h5seurat", assays = "RNA")

# specify the assays as RNA or there will be an error

# Visualising the tabula sapiens object

DimPlot(obj, reduction = "umap")

obj@meta.data

## Annotating the harmony dataset with cell labels from Tabula Sapiens

# Reading in the harmony data

setwd('C:/Users/Home/Documents/Dissertation/Carotid_Datasets/Harmony')


Integrated <- readRDS('Carotid_Harmonised.rds')

# Adding features to Tabula Sapiens


Indents(obj)

# Pre-Processing of tabula sapiens: Normalising data, etc.


# perform standard workflow steps to figure out if we see any batch effects --------
obj <- NormalizeData(object = obj)
obj <- FindVariableFeatures(object = obj)
obj <- ScaleData(object = obj)
obj <- RunPCA(object = obj)
ElbowPlot(obj)

# Running UMAP and clustering
obj <- FindNeighbors(object = obj, dims = 1:30)
obj <- FindClusters(object = obj)
obj <- RunUMAP(object = obj, dims = 1:30)

# Visualising the UMAP clusters 

DimPlot(obj, reduction = 'umap')

Idents(obj)


# Labelling the clusters

## Each cell has been asigned to a cell type in the metadata under 'cell ontology class'
view(obj@meta.data)

## Settings idents (clusters) or the tabula reference as the Seurat annotations provided
Idents(obj)
Idents(obj) <- obj@meta.data$cell_ontology_class
Idents(obj)

DimPlot(obj, reduction = 'umap', label = TRUE)

# Saving the tabula sapiens reference as an RDS
setwd('C:/Users/Home/Documents/Dissertation/Carotid_Annotation/')

saveRDS(obj, file = "TabulaSapiensVasculature.rds")
        

