# script to integrate scRNA-Seq datasets to correct for batch effects
setwd("/Users/judepops/Documents/Dissertation_Local/Data")


# load libraries
library(Seurat)
library(ggplot2)
library(tidyverse)
library(gridExtra)

# get data location
dirs <- list.dirs(path = 'Fernandes/', recursive = F, full.names = F)

for(x in dirs){
  #Assigning the name to each file as the SampleX
  name <- gsub('_filtered_feature_bc_matrix','', x)
  
  #Creating our count matrix as cts
  cts <- ReadMtx(mtx = paste0('Fernandes/',x,'/matrix.mtx.gz'),
                 features = paste0('Fernandes/',x,'/genes.tsv.gz'),
                 cells = paste0('Fernandes/',x,'/barcodes.tsv.gz'))
  
  # create seurat objects
  assign(name, CreateSeuratObject(counts = cts))
}




# merge datasets into one object for faster processing
# first parameter is the first object (largest) and the rest are provided as vector of objects
merged_fernandes <- merge(Sample1, y = c(Sample2, `Sample3 ` , `Sample3A `, `Sample4  `,
                                      `Sample5 `),
                       add.cell.ids = ls()[12:17], #when all objects are merged want to identify which cell barcodes are from which dataset
                       project = 'Fern')


merged_fernandes

# QC & filtering ----------------------- (on the merged object)

View(merged_fernandes@meta.data)


# calculate mitochondrial percentage
merged_fernandes$mitoPercent <- PercentageFeatureSet(merged_fernandes, pattern='^MT-')

# explore QC


# filtering
merged_fernandes_filtered <- subset(merged_fernandes, subset = nCount_RNA > 800 &
                                   nFeature_RNA > 500 &
                                   mitoPercent < 10)

merged_fernandes_filtered

merged_fernandes


# perform standard workflow steps to figure out if we see any batch effects --------
merged_fernandes_filtered <- NormalizeData(object = merged_fernandes_filtered)
merged_fernandes_filtered <- FindVariableFeatures(object = merged_fernandes_filtered)
merged_fernandes_filtered <- ScaleData(object = merged_fernandes_filtered)
merged_fernandes_filtered <- RunPCA(object = merged_fernandes_filtered)
ElbowPlot(merged_fernandes_filtered)
merged_fernandes_filtered <- FindNeighbors(object = merged_fernandes_filtered, dims = 1:20)
merged_fernandes_filtered <- FindClusters(object = merged_fernandes_filtered)
merged_fernandes_filtered <- RunUMAP(object = merged_fernandes_filtered, dims = 1:20)

view(merged_fernandes_filtered)

#merged_fernandes_filtered <- readRDS("/Users/judepops/Documents/Dissertation_Local/Output/fernandes/fernandes_Filtered.rds")
merged_fernandes_filtered$Source <- 'fernandes'
saveRDS(merged_fernandes_filtered, file = "/Users/judepops/Documents/Dissertation_Local/Output/fernandes/Fernandes_Filtered_Final.rds")
#merged_fernandes_filtered <- readRDS("/Users/judepops/Documents/Dissertation_Local/Output/fernandes/fernandes_Filtered.rds")




################################################################################################################################################




# Runnign it again for python h5 object 
setwd("/Users/judepops/Documents/Dissertation_Local/Data")



merged_fernandes

# QC & filtering ----------------------- (on the merged object)

View(merged_fernandes@meta.data)


# calculate mitochondrial percentage
merged_fernandes$mitoPercent <- PercentageFeatureSet(merged_fernandes, pattern='^MT-')

# explore QC


# filtering
merged_fernandes_filtered <- subset(merged_fernandes, subset = nCount_RNA > 800 &
                                      nFeature_RNA > 500 &
                                      mitoPercent < 10)

merged_fernandes_filtered

merged_fernandes


# perform standard workflow steps to figure out if we see any batch effects --------
merged_fernandes_filtered <- NormalizeData(object = merged_fernandes_filtered)
merged_fernandes_filtered <- FindVariableFeatures(object = merged_fernandes_filtered)
merged_fernandes_filtered <- ScaleData(object = merged_fernandes_filtered)
merged_fernandes_filtered <- RunPCA(object = merged_fernandes_filtered)
ElbowPlot(merged_fernandes_filtered)
merged_fernandes_filtered <- RunUMAP(object = merged_fernandes_filtered, dims = 1:20)



merged_fernandes_filtered_h5 <- merged_fernandes_filtered

## Preparing the h5



# Loading the libraries
library(Seurat)
library(SeuratData)
library(SeuratDisk)

# Converting to h5

setwd("/Users/judepops/Documents/Dissertation_Local/Output/fernandes")


SaveH5Seurat(merged_fernandes_filtered_h5, filename = "/Users/judepops/Documents/Dissertation_Local/Output/fernandes/fernandesh5.h5Seurat")
Convert("fernandesh5.h5Seurat", dest = "h5ad")


