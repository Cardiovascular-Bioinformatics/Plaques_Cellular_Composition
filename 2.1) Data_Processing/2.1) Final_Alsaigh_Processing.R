# script to integrate scRNA-Seq datasets to correct for batch effects
setwd("/Users/judepops/Documents/Dissertation_Local/Data")


# load libraries
library(Seurat)
library(ggplot2)
library(tidyverse)
library(gridExtra)

# get data location
dirs <- list.dirs(path = 'Alsaigh/', recursive = F, full.names = F)

for(x in dirs){
  #Assigning the name to each file as the SampleX
  name <- gsub('_filtered_feature_bc_matrix','', x)
  
  #Creating our count matrix as cts
  cts <- ReadMtx(mtx = paste0('Alsaigh/',x,'/matrix.mtx.gz'),
                 features = paste0('Alsaigh/',x,'/features.tsv.gz'),
                 cells = paste0('Alsaigh/',x,'/barcodes.tsv.gz'))
  
  # create seurat objects
  assign(name, CreateSeuratObject(counts = cts))
}




# merge datasets into one object for faster processing
# first parameter is the first object (largest) and the rest are provided as vector of objects
merged_alsaigh <- alsaigh


# QC & filtering ----------------------- (on the merged object)

View(merged_alsaigh@meta.data)


# calculate mitochondrial percentage
merged_alsaigh$mitoPercent <- PercentageFeatureSet(merged_alsaigh, pattern='^MT-')

# explore QC


# filtering
merged_alsaigh_filtered <- subset(merged_alsaigh, subset = nCount_RNA > 800 &
                                    nFeature_RNA > 500 &
                                    mitoPercent < 10)

merged_alsaigh_filtered

merged_alsaigh


# perform standard workflow steps to figure out if we see any batch effects --------
merged_alsaigh_filtered <- NormalizeData(object = merged_alsaigh_filtered)
merged_alsaigh_filtered <- FindVariableFeatures(object = merged_alsaigh_filtered)
merged_alsaigh_filtered <- ScaleData(object = merged_alsaigh_filtered)
merged_alsaigh_filtered <- RunPCA(object = merged_alsaigh_filtered)
ElbowPlot(merged_alsaigh_filtered)
merged_alsaigh_filtered <- FindNeighbors(object = merged_alsaigh_filtered, dims = 1:20)
merged_alsaigh_filtered <- FindClusters(object = merged_alsaigh_filtered)
merged_alsaigh_filtered <- RunUMAP(object = merged_alsaigh_filtered, dims = 1:20)

view(merged_alsaigh_filtered)

#merged_alsaigh_filtered <- readRDS("/Users/judepops/Documents/Dissertation_Local/Output/alsaigh/alsaigh_Filtered.rds")
merged_alsaigh_filtered$Source <- 'alsaigh'
saveRDS(merged_alsaigh_filtered, file = "/Users/judepops/Documents/Dissertation_Local/Output/alsaigh/Alsaigh_Filtered_Final.rds")


#merged_alsaigh_filtered <- readRDS("/Users/judepops/Documents/Dissertation_Local/Output/alsaigh/alsaigh_Filtered.rds")


## Preparing the h5

# Loading the libraries
library(Seurat)
library(SeuratData)
library(SeuratDisk)

# Converting to h5

setwd("/Users/judepops/Documents/Dissertation_Local/Output/alsaigh")


SaveH5Seurat(merged_alsaigh_filtered, filename = "/Users/judepops/Documents/Dissertation_Local/Output/alsaigh/alsaigh.h5Seurat")
Convert("alsaigh.h5Seurat", dest = "h5ad")






################################################################################################################################################




# Runnign it again for python h5 object 

# script to integrate scRNA-Seq datasets to correct for batch effects
setwd("/Users/judepops/Documents/Dissertation_Local/Data")


# load libraries
library(Seurat)
library(ggplot2)
library(tidyverse)
library(gridExtra)

# get data location
dirs <- list.dirs(path = 'Alsaigh/', recursive = F, full.names = F)

for(x in dirs){
  #Assigning the name to each file as the SampleX
  name <- gsub('_filtered_feature_bc_matrix','', x)
  
  #Creating our count matrix as cts
  cts <- ReadMtx(mtx = paste0('Alsaigh/',x,'/matrix.mtx.gz'),
                 features = paste0('Alsaigh/',x,'/features.tsv.gz'),
                 cells = paste0('Alsaigh/',x,'/barcodes.tsv.gz'))
  
  # create seurat objects
  assign(name, CreateSeuratObject(counts = cts))
}




# merge datasets into one object for faster processing
# first parameter is the first object (largest) and the rest are provided as vector of objects
merged_alsaigh <- alsaigh


# QC & filtering ----------------------- (on the merged object)

View(merged_alsaigh@meta.data)


# calculate mitochondrial percentage
merged_alsaigh$mitoPercent <- PercentageFeatureSet(merged_alsaigh, pattern='^MT-')

# explore QC


# filtering
merged_alsaigh_filtered <- subset(merged_alsaigh, subset = nCount_RNA > 800 &
                                    nFeature_RNA > 500 &
                                    mitoPercent < 10)

merged_alsaigh_filtered

merged_alsaigh


# perform standard workflow steps to figure out if we see any batch effects --------
merged_alsaigh_filtered <- NormalizeData(object = merged_alsaigh_filtered)
merged_alsaigh_filtered <- FindVariableFeatures(object = merged_alsaigh_filtered)
merged_alsaigh_filtered <- ScaleData(object = merged_alsaigh_filtered)
merged_alsaigh_filtered <- RunPCA(object = merged_alsaigh_filtered)
ElbowPlot(merged_alsaigh_filtered)
merged_alsaigh_filtered <- RunUMAP(object = merged_alsaigh_filtered, dims = 1:20)


setwd('/Users/judepops/Documents/DISSERTATION/Processing/')

saveRDS(merged_alsaigh_filtered, 'AlsaighRDS.rds')




## Preparing the h5
merged_alsaigh_filtered_h5 <- merged_alsaigh_filtered

# Loading the libraries
library(Seurat)
library(SeuratData)
library(SeuratDisk)

# Converting to h5

setwd("/Users/judepops/Documents/Dissertation_Local/Output/alsaigh")


SaveH5Seurat(merged_alsaigh_filtered_h5, filename = "/Users/judepops/Documents/Dissertation_Local/Output/alsaigh/alsaighh5.h5Seurat")
Convert("alsaighh5.h5Seurat", dest = "h5ad")



