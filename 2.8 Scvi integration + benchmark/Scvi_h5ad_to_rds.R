#Working with loom

#Packages
# install scater https://bioconductor.org/packages/release/bioc/html/scater.html
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("scater")

library(scater)
library(Seurat)
# install SeuratDisk from GitHub using the remotes package remotes::install_github(repo =
# 'mojaveazure/seurat-disk', ref = 'develop')
library(SeuratDisk)
library(SeuratData)
library(patchwork)

# Install devtools from CRAN
install.packages("devtools")
# Use devtools to install hdf5r and loomR from GitHub



library(loomR)
#File directory
lfile <- connect(filename = "/Users/judepops/Downloads/scvi_integration.loom", mode = "r+", skip.validate = TRUE) 
seurat_loom <- as.Seurat(lfile)




# Converting a file from h5ad to suerat

# Loading the libraries
library(Seurat)
library(SeuratData)
library(SeuratDisk)

# Converting to h5

setwd("/Users/judepops/Documents/DISSERTATION/Scvi/")



Convert("test.h5ad", dest = "h5seurat", overwrite = TRUE)
scvi <- LoadH5Seurat("test.h5seurat", assays = 'RNA', meta.data = FALSE, misc = FALSE)
scvi

#Adding metadata to the object
scvi
metadata <- read.csv('2_scvi_metadata.csv')
scvi@meta.data <-cbind(scvi@meta.data,metadata)
scvi@meta.data
DimPlot(scvi, group.by = 'leiden')


#Saving the object
saveRDS(scvi, 'Scvi_Integrated_Leiden_Louvain.rds')

DimPlot(scvi, group.by = 'louvain_0.8')










#WORKS ######################################################################################################


# trying it again 

setwd("/Users/judepops/Documents/DISSERTATION/Scvi/")

Convert("4_scvi_integrated.h5ad", dest = "h5seurat", overwrite = TRUE)
scvi <- LoadH5Seurat("4_scvi_integrated.h5seurat", assays = 'RNA')
scvi

metadata <- read.csv('final_scvi_metadata.csv')
scvi@meta.data <-cbind(scvi@meta.data,metadata)
scvi@meta.data
DimPlot(scvi, group.by = 'leiden')


saveRDS(scvi, 'FINAL_Scvi_Integrated_Leiden_Louvain.rds')



rownames(scvi)








ERROR ######################################################################################################

setwd("/Users/judepops/Documents/DISSERTATION/Scvi/")

Convert("5_scvi_integrated.h5ad", dest = "h5seurat", overwrite = TRUE)
scvi <- LoadH5Seurat("5_scvi_integrated.h5seurat", assays = 'RNA')
scvi
scvi@meta.data





metadata <- read.csv('final_scvi_metadata.csv')
scvi@meta.data <-cbind(scvi@meta.data,metadata)
scvi@meta.data
DimPlot(scvi, group.by = 'leiden')


saveRDS(scvi, 'FINAL_Scvi_Integrated_Leiden_Louvain.rds')
