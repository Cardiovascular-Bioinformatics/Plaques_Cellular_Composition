# Harmony Integration
library(harmony)
library(Seurat)
library(SeuratData)
library(tidyverse)
library(ggplot2)

# Loading the data
fernandes <- readRDS('/Users/judepops/Documents/Dissertation_Local/Output/Fernandes/Fernandes_Filtered_Final.rds')
alsaigh <- readRDS('/Users/judepops/Documents/Dissertation_Local/Output/Alsaigh/Alsaigh_Filtered_Final.rds')
pan <- readRDS('/Users/judepops/Documents/Dissertation_Local/Output/Pan/Pan_Filtered_Final.rds')


fernandes_csv <- read.csv('fernandes.csv')
pan_csv <- read.csv('pan.csv')
alsaigh_csv <- read.csv('alsaigh.csv')



# Merging the data into one big seurat object

## Merging into carotid

carotid <- merge(alsaigh, y = c(fernandes,pan),
                 project = 'Merge')

## saving the merged data

saveRDS(carotid, file = 'Pan_Fern_Als.rds')


## Clustering data and running pca
carotid <- NormalizeData(object = carotid)
carotid <- FindVariableFeatures(object = carotid)
carotid <- ScaleData(object = carotid)
carotid <- RunPCA(object = carotid)
ElbowPlot(carotid)
carotid <- FindNeighbors(carotid, dims = 1:10)
carotid <- FindClusters(carotid, resolution = 0.5)
carotid <- RunUMAP(object = carotid, dims = 1:20)

# Visualising umap
DimPlot(carotid, reduction = "umap", group.by = 'Source')


# run Harmony -----------
## returns corrected dimensionality reductions = 'embeddings'

carotid.harmony <- carotid %>%
  RunHarmony(group.by.vars = 'Source', plot_convergence = FALSE)

#Harmony returns corrected dimensionality values as embeddings (not an altered expression matrix)

#Results are saved within the reduction slot of harmony
carotid.harmony@reductions


#How to get the embeddings: provide name of object and the reduction
carotid.harmony.embed <- Embeddings(carotid.harmony, "harmony")
carotid.harmony.embed[1:10,1:10]

#Just like there are multiple harmony components there are multiple principle components. 
#The embeddings can be used to perform non linear dimensionality reductions (UMAP) and
#downstream analysis (clustering)


carotid.harmony@meta.data

# Using Harmony embeddings to normalise

carotid.harmony <- carotid.harmony %>%
  RunUMAP(reduction = 'harmony', dims = 1:20) %>% #the first 20 harmony dimensions 
  #(not PCA dimensions) --> remains the same
  FindNeighbors(reduction = "harmony", dims = 1:20) %>%
  FindClusters(resolution = 0.5)

carotid.harmony <- carotid.harmony %>%
  RunUMAP(reduction = 'harmony', dims = 1:20) %>% #the first 20 harmony dimensions 
  #(not PCA dimensions) --> remains the same
  FindNeighbors(reduction = "harmony", dims = 1:20) %>%
  FindClusters(resolution = 0.8)

# Saving integrated harmony dataset

saveRDS(carotid.harmony, file = '/Users/judepops/Documents/DISSERTATION/Integrating/Harmony_Integration.rds')


# Converting integrated harmony dataset to h5ad file

# Loading the libraries
library(Seurat)
library(SeuratData)
library(SeuratDisk)

# Converting to h5

setwd("/Users/judepops/Documents/DISSERTATION/Integrating/")


SaveH5Seurat(harmony, filename = "/Users/judepops/Documents/DISSERTATION/Integrating/Harmony_Integration.h5Seurat")
Convert("Harmony_Integration.h5Seurat", dest = "h5ad")



SaveH5Seurat(carotid, filename = "/Users/judepops/Documents/DISSERTATION/Integrating/Pan_Fern_Als.h5Seurat")
Convert("Pan_Fern_Als.h5Seurat", dest = "h5ad")


