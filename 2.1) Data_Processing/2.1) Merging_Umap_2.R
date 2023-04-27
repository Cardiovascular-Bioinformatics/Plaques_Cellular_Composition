# Libraries
library(harmony)
library(Seurat)
library(SeuratData)
library(tidyverse)
library(ggplot2)

# Loading the data
fernandes <- readRDS('/Users/judepops/Documents/Dissertation_Local/Output/Fernandes/Fernandes_Filtered_Final.rds')
alsaigh <- readRDS('/Users/judepops/Documents/Dissertation_Local/Output/Alsaigh/Alsaigh_Filtered_Final.rds')
pan <- readRDS('/Users/judepops/Documents/Dissertation_Local/Output/Pan/Pan_Filtered_Final.rds')

# Merging the data into one big seurat object

## Merging into carotid

carotid <- merge(alsaigh, y = c(fernandes,pan),
                 add.cell.ids = c('Alsaigh','Fernandes','Pan'),
                 project = 'Merge')



# Qc

## Checking Mitochondrial features

carotid$mito.percent <- PercentageFeatureSet(carotid, pattern = '^MT-')
head(carotid@meta.data)


# Filtering out samples with low gene expression/high mitochondrial gene expression

carotid.filtered <- subset(carotid, subset = nCount_RNA > 800 &
                             nFeature_RNA > 200 & 
                             mito.percent < 5)

# standard workflow steps - to visualise the data
carotid.filtered <- NormalizeData(carotid.filtered)
carotid.filtered <- FindVariableFeatures(carotid.filtered)
carotid.filtered <- ScaleData(carotid.filtered)
carotid.filtered <- RunPCA(carotid.filtered)
ElbowPlot(carotid.filtered)
carotid.filtered <- RunUMAP(carotid.filtered, dims = 1:20, reduction = 'pca')

#Checking before Umap plot

before <- DimPlot(carotid.filtered, reduction = 'umap', group.by = 'Source')

cells_umap <- DimPlot(carotid.filtered, reduction = "umap", label = TRUE,
              repel = TRUE)

# Saving the final object as merged and filtered 

saveRDS(carotid.filtered, file = '/Users/judepops/Documents/Dissertation_Local/Output/Merged/Merged_Carotid_Final.rds')

# Converting the final object into a h5 file for running in python

Carotid_Atlas <- readRDS('/Users/judepops/Documents/Dissertation_Local/Output/Merged/Merged_Carotid_Final.rds')

# Libraries
library(Seurat)
library(SeuratData)
library(SeuratDisk)

# Converting to h5

## Setting wd
setwd("/Users/judepops/Documents/Dissertation_Local/Output/Merged")

## Converting
SaveH5Seurat(Carotid_Atlas, filename = "/Users/judepops/Documents/Dissertation_Local/Output/Merged/Carotid_Atlas.h5Seurat")
Convert("Carotid_Atlas.h5Seurat", dest = "h5ad")


