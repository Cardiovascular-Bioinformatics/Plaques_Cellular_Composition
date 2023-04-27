#KNN/SNN clustering on each of the datasets - 0.8 resolution

# Packages
library(Seurat)
library(tidyverse)
library(ggplot2)
library(tidyverse)
library(gridExtra)

# Loading in the processed datasets (normalised and filtered)
fernandes <- readRDS('/Users/judepops/Documents/DISSERTATION/Processing/FernandesRDS.rds')
alsaigh <- readRDS('/Users/judepops/Documents/DISSERTATION/Processing/AlsaighRDS.rds')
pan <- readRDS('/Users/judepops/Documents/DISSERTATION/Processing/PanRDS.rds')


# Running kNN clustering on each of the experimental datasets 
# PCs calculated and normalisation already done during processing step

fernandes <- FindNeighbors(fernandes, dims = 1:10)
fernandes <- FindClusters(fernandes, resolution = 0.8)
saveRDS(fernandes, '/Users/judepops/Documents/DISSERTATION/Processing/FernandesRDS.rds')
fernandes <- fernandes@meta.data 
fernandes <- fernandes %>% rownames_to_column(var = 'X')



alsaigh <- FindNeighbors(alsaigh, dims = 1:10)
alsaigh <- FindClusters(alsaigh, resolution = 0.8)
saveRDS(alsaigh, '/Users/judepops/Documents/DISSERTATION/Processing/AlsaighRDS.rds')
alsaigh <- alsaigh@meta.data 
alsaigh <- alsaigh %>% rownames_to_column(var = 'X')


pan <- FindNeighbors(pan, dims = 1:10)
pan <- FindClusters(pan, resolution = 0.8)
saveRDS(pan, '/Users/judepops/Documents/DISSERTATION/Processing/PanRDS.rds')

pan <- pan@meta.data 
pan <- pan %>% rownames_to_column(var = 'X')

# Exporting metadata to csv files
setwd('/Users/judepops/Documents/Untitled Folder')

write_csv(fernandes, 'fernandezknn.csv')
write_csv(alsaigh, 'alsaighknn.csv')
write_csv(pan, 'panknn.csv')






source.cluster <- function(obj) {
  obj <- FindNeighbors(obj, dims = 1:10)
  obj <- FindClusters(obj, resolution = 0.8)
  obj
}



fernandes <- source.cluster(fernandes)












## Checking Integrated data set clustering



# Harmony

harmony_integration <- readRDS('/Users/judepops/Documents/DISSERTATION/Integrating/Harmony_Integration.rds')
harmony_integration@reductions

DimPlot(harmony_integration, reduction = 'umap')
DimPlot(harmony_integration, reduction = 'umap', group.by = 'Source')

setwd('/Users/judepops/Documents/Untitled Folder')
harmony_integration_csv <- harmony_integration@meta.data
harmony_integration_csv <- harmony_integration_csv %>% rownames_to_column(var = 'X')
write_csv(harmony_integration_csv, 'harmony_0.8.csv')

# Clustering at 0.8 resolution

harmony_integration <- harmony_integration %>% 
  RunUMAP(reduction = "harmony", dims = 1:10) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:10) %>% 
  FindClusters(resolution = 0.8) 


# rPCA

rPCA_integration <- readRDS('/Users/judepops/Documents/DISSERTATION/Integrated_Leiden_Louvain0.8/rPCA_Integrated_Leiden_Louvain.rds')
DimPlot(rPCA_integration, reduction = 'umap', label=TRUE)
DimPlot(rPCA_integration, reduction = 'umap', group.by = 'Source')
DimPlot(rPCA_integration, reduction = 'umap', group.by = 'leiden')


# Clustering at 0.8 resolution
DefaultAssay(rPCA_integration) <- "integrated"
rPCA_integration <- FindNeighbors(rPCA_integration, reduction = "pca", dims = 1:10)
rPCA_integration <- FindClusters(rPCA_integration, resolution = 0.8)

setwd('/Users/judepops/Documents/Untitled Folder')
rPCA_integration_csv <- rPCA_integration@meta.data
rPCA_integration_csv <- rPCA_integration_csv %>% rownames_to_column(var = 'X')
write_csv(rPCA_integration_csv, 'rPCA_0.8.csv')


saveRDS(rPCA_integration, '/Users/judepops/Documents/DISSERTATION/Integrated_Leiden_Louvain0.8/rPCA_Integrated_Leiden_Louvain.rds')



# Scanorama

Scanorama_integration <- readRDS('/Users/judepops/Documents/DISSERTATION/Integrated_Leiden_Louvain0.8/Panorama_Integrated_Leiden_Louvain.rds')
DimPlot(Scanorama_integration, reduction = 'umap', label=TRUE)
DimPlot(Scanorama_integration, reduction = 'umap', group.by = 'Source')
DimPlot(Scanorama_integration, )


#Clustering at 0.8 resolution
Scanorama_integration <- ScaleData(Scanorama_integration)
Scanorama_integration <- RunUMAP(Scanorama_integration, dims = 1:10)
Scanorama_integration <- FindNeighbors(Scanorama_integration, dims = 1:10) 
Scanorama_integration <- FindClusters(Scanorama_integration, resolution = c(0.5,0.8))

saveRDS(Scanorama_integration, '/Users/judepops/Documents/DISSERTATION/Integrated_Leiden_Louvain0.8/Panorama_Integrated_Leiden_Louvain.rds')

setwd('/Users/judepops/Documents/Untitled Folder')
Scanorama_integration_csv <- Scanorama_integration@meta.data
Scanorama_integration_csv <- Scanorama_integration_csv %>% rownames_to_column(var = 'X')
write_csv(Scanorama_integration_csv, 'panorama_0.8.csv')

