# script to integrate scRNA-Seq datasets to correct for batch effects
setwd("/Users/judepops/Documents/Dissertation_Local/Data")


# load libraries
library(Seurat)
library(ggplot2)
library(tidyverse)
library(gridExtra)




#Assigning the name to each file as the SampleX
name <- "RPE004_matrix"

#Creating our count matrix as cts
cts <- read.csv(file = '/Users/judepops/Documents/Dissertation_Local/Data/Pan/RPE004_filtered_feature_bc_matrix /matrix.txt', sep="\t", header=TRUE, row.names=1)

assign(name, CreateSeuratObject(counts = cts))




name <- "RPE005_matrix"

#Creating our count matrix as cts
cts <- read.csv(file = '/Users/judepops/Documents/Dissertation_Local/Data/Pan/RPE005_filtered_feature_bc_matrix /matrix.txt', sep="\t", header=TRUE, row.names=1)

assign(name, CreateSeuratObject(counts = cts))




name <- "RPE006_matrix"

#Creating our count matrix as cts
cts <- read.csv(file = '/Users/judepops/Documents/Dissertation_Local/Data/Pan/RPE006_filtered_feature_bc_matrix /matrix.txt', sep="\t", header=TRUE, row.names=1)

assign(name, CreateSeuratObject(counts = cts))




# merge datasets into one object for faster processing
# first parameter is the first object (largest) and the rest are provided as vector of objects
merged_pan <- merge(`RPE004_matrix`, y = c(`RPE005_matrix`,`RPE006_matrix`),
                    #when all objects are merged want to identify which cell barcodes are from which dataset
                    project = 'Als')


merged_pan

# QC & filtering ----------------------- (on the merged object)

View(merged_pan@meta.data)


# calculate mitochondrial percentage
merged_pan$mitoPercent <- PercentageFeatureSet(merged_pan, pattern='^MT-')

# explore QC


# filtering
merged_pan_filtered <- subset(merged_pan, subset = nCount_RNA > 800 &
                                nFeature_RNA > 500 &
                                mitoPercent < 10)

merged_pan_filtered

merged_pan


# perform standard workflow steps to figure out if we see any batch effects --------
merged_pan_filtered <- NormalizeData(object = merged_pan_filtered)
merged_pan_filtered <- FindVariableFeatures(object = merged_pan_filtered)
merged_pan_filtered <- ScaleData(object = merged_pan_filtered)
merged_pan_filtered <- RunPCA(object = merged_pan_filtered)
ElbowPlot(merged_pan_filtered)
merged_pan_filtered <- FindNeighbors(object = merged_pan_filtered, dims = 1:20)
merged_pan_filtered <- FindClusters(object = merged_pan_filtered)
merged_pan_filtered <- RunUMAP(object = merged_pan_filtered, dims = 1:20)

view(merged_pan_filtered)

merged_pan_filtered$Source <- 'Pan'

#merged_pan_filtered <- readRDS("/Users/judepops/Documents/Dissertation_Local/Output/Pan/Pan_Filtered.rds")
saveRDS(merged_pan_filtered, file = "/Users/judepops/Documents/Dissertation_Local/Output/Pan/Pan_Filtered_Final.rds")
#merged_pan_filtered <- readRDS(file = "/Users/judepops/Documents/Dissertation_Local/Output/Pan/Pan_Filtered.rds")

# #Converting the file to a anndata object
# 
# library(reticulate)
# 
# Pan_Filtered_ad <- Convert(from=merged_pan_filtered, to="anndata", filename="Pan_Filtered.h5ad")



#########################################################################################################################








# Creating a Pan file for h5ad python


# script to integrate scRNA-Seq datasets to correct for batch effects
setwd("/Users/judepops/Documents/Dissertation_Local/Data")


# load libraries
library(Seurat)
library(ggplot2)
library(tidyverse)
library(gridExtra)




#Assigning the name to each file as the SampleX
name <- "RPE004_matrix"

#Creating our count matrix as cts
cts <- read.csv(file = '/Users/judepops/Documents/Dissertation_Local/Data/Pan/RPE004_filtered_feature_bc_matrix /matrix.txt', sep="\t", header=TRUE, row.names=1)

assign(name, CreateSeuratObject(counts = cts))




name <- "RPE005_matrix"

#Creating our count matrix as cts
cts <- read.csv(file = '/Users/judepops/Documents/Dissertation_Local/Data/Pan/RPE005_filtered_feature_bc_matrix /matrix.txt', sep="\t", header=TRUE, row.names=1)

assign(name, CreateSeuratObject(counts = cts))




name <- "RPE006_matrix"

#Creating our count matrix as cts
cts <- read.csv(file = '/Users/judepops/Documents/Dissertation_Local/Data/Pan/RPE006_filtered_feature_bc_matrix /matrix.txt', sep="\t", header=TRUE, row.names=1)

assign(name, CreateSeuratObject(counts = cts))




# merge datasets into one object for faster processing
# first parameter is the first object (largest) and the rest are provided as vector of objects
merged_pan <- merge(`RPE004_matrix`, y = c(`RPE005_matrix`,`RPE006_matrix`),
                    #when all objects are merged want to identify which cell barcodes are from which dataset
                    project = 'Als')


merged_pan

# QC & filtering ----------------------- (on the merged object)

View(merged_pan@meta.data)


# calculate mitochondrial percentage
merged_pan$mitoPercent <- PercentageFeatureSet(merged_pan, pattern='^MT-')

# explore QC


# filtering
merged_pan_filtered <- subset(merged_pan, subset = nCount_RNA > 800 &
                                nFeature_RNA > 500 &
                                mitoPercent < 10)

merged_pan_filtered

merged_pan


# perform standard workflow steps to figure out if we see any batch effects --------
merged_pan_filtered <- NormalizeData(object = merged_pan_filtered)
merged_pan_filtered <- FindVariableFeatures(object = merged_pan_filtered)
merged_pan_filtered <- ScaleData(object = merged_pan_filtered)
merged_pan_filtered <- RunPCA(object = merged_pan_filtered)
ElbowPlot(merged_pan_filtered)
merged_pan_filtered <- RunUMAP(merged_pan_filtered, dims = 1:20, reduction = 'pca')


setwd('/Users/judepops/Documents/DISSERTATION/Processing/')

saveRDS(merged_pan_filtered, 'PanRDS.rds')
# Saving merged_pan_filtered_h5

merged_pan_filtered_h5 <- merged_pan_filtered


#saving file


#Converting the file to adnndata

library(Seurat)
library(SeuratData)
library(SeuratDisk)

setwd("/Users/judepops/Documents/Dissertation_Local/Output/Pan")


SaveH5Seurat(merged_pan_filtered_h5, filename = "/Users/judepops/Documents/Dissertation_Local/Output/Pan/Panh5.h5Seurat")
Convert("Panh5.h5Seurat", dest = "h5ad")


