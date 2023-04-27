install.packages("devtools")
devtools::install_github("immunogenomics/lisi")

library(lisi)

?compute_lisi


## Example with 400 cells. 
library(lisi)
library(dplyr)
library(tidyr)
library(ggplot2)
library(magrittr)
library(stringr)


#Calculate LISI metric for each integration using a) Leiden, b) default Seurat metric. 
#to do this you need to use the codes and instructions of the following link: https://github.com/immunogenomics/LISI. 
#This takes as input a) one matrix with the UMAP dimensions of the integrated dataset 
#b) a data frame with categorical variables (one row for each cell) 
#(in your case it will have one column with the clustering information 
#and one column with the dataset that each cell comes from (eg. 1 for Pan, 2 for Alsaigh and 3 for Fernandez) 
#c) a vector with the columns names of the dataframe of the second argument. 
#The compute_lisi method will return you a dataframe with 2 columns and rows 
#equal with the number of cells in the integrated object. You should be calculating 
#first the average per row and the final LISI score of the integration will be the average of the row averages

#Computing LISI for the default seurat clustering (0.8)


setwd('C:/Users/Home/Documents/Dissertation/Carotid_Datasets/Harmony/')

# Integrated <- readRDS('Carotid_Harmonised.rds')
# 
# filenames <- list.files(pattern = ".rds", full.names = TRUE)
# 
# 
# # Getting the UMAP dimesions from the integrated data
# 
# Umap <- as.matrix(Integrated[["umap"]]@cell.embeddings)
# 
# # Clustering information of each cell
# 
# metadata <- as.data.frame(Integrated@meta.data)
# metadata <- metadata %>% select(RNA_snn_res.0.8, Source)
# # metadata$Source <- str_replace(metadata$Source, "fernandes", "1")
# # metadata$Source <- str_replace(metadata$Source, "alsaigh", "2")
# # metadata$Source <- str_replace(metadata$Source, "Pan", "3")
# # metadata$Source <- as.character(metadata$Source)
# 
# 
# # Column names for the metadata argument
# 
# vector <- as.vector(colnames(metadata))
# 
# #Computing LISI
# res <- compute_lisi(Umap, metadata, vector)

################## ######### ######### ######### ######### ######### ######### 


#For loop

setwd('C:/Users/Home/Documents/Dissertation/FInal_Integrated/')

harmony <- readRDS('Harmony_Integrated_Leiden_Louvain.rds')
scanorama <- readRDS('Panorama_Integrated_Leiden_Louvain.rds')
rPCA <- readRDS('rPCA_Integrated_Leiden_Louvain.rds')


names = c('harmony', 'scanorama', 'rPCA')

for(i in names) {
  Umap <- as.matrix(get(i)[["umap"]]@cell.embeddings)
  metadata <- as.data.frame(get(i)@meta.data)
  metadata <- metadata %>% select(RNA_snn_res.0.8, Source)
  res <- compute_lisi(Umap, metadata, c('RNA_snn_res.0.8', 'Source'))

}



# Writing a funciton for LISI

#Harmony

do.LISI.harm <- function(obj) {
  Umap <- as.matrix(obj[["umap"]]@cell.embeddings)
  metadata <- as.data.frame(obj@meta.data)
  metadata <- metadata %>% select(RNA_snn_res.0.8, Source)
  res <- compute_lisi(Umap, metadata, c('RNA_snn_res.0.8', 'Source'))
  res
}

LISI_harmony_0.8 <- do.LISI.harm(harmony)






#Scanorama


#rPCA


S1 <- do.LISI(S1)
S2 <- do.LISI(S2)



do.seurat <- function(obj) {
  obj[['percent.mt']] <- PercentageFeatureSet(...)
  ...  # rest of "traditional Seurat analysis"
  obj
}

S1 <- do.seurat(S1)
S2 <- do.seurat(S2)



# Getting the UMAP dimesions from the integrated data

X <- as.matrix(harmony[["umap"]]@cell.embeddings)
# Clustering information of each cell

metadata <- as.data.frame(harmony@meta.data)
metadata <- metadata %>% select(Source, RNA_snn_res.0.8)
metadata$Source <- str_replace(metadata$Source, "fernandes", "1")
metadata$Source <- str_replace(metadata$Source, "alsaigh", "2")
metadata$Source <- str_replace(metadata$Source, "Pan", "3")
metadata$Source <- as.integer(metadata$Source)



table(metadata$RNA_snn_res.0.8)
table(metadata$Source)

head(metadata)
# Column names for the metadata argument

vector <- as.vector(colnames(metadata))

#Computing LISI
res <- compute_lisi(X, metadata, c('Source', 'RNA_snn_res.0.8'))


#Calculating row averages

res$Avg_LISI <- rowMeans(res, na.rm=TRUE)


# Mean

LISI <- mean(res$Avg_LISI)

















#########################################################  louvain Function

LISI.harm <- function(obj) {
  Umap <- as.matrix(obj[["umap"]]@cell.embeddings)
  metadata <- as.data.frame(obj@meta.data)
  metadata <- metadata %>% select(RNA_snn_res.0.8, Source)
  res <- compute_lisi(Umap, metadata, c('RNA_snn_res.0.8', 'Source'))
  res$Avg_LISI <- rowMeans(res, na.rm=TRUE)
  output <- mean(res$Avg_LISI)
  output
}

LISI.scan <- function(obj) {
  Umap <- as.matrix(obj[["umap"]]@cell.embeddings)
  metadata <- as.data.frame(obj@meta.data)
  metadata <- metadata %>% select(pano_snn_res.0.8, Source)
  res <- compute_lisi(Umap, metadata, c('pano_snn_res.0.8', 'Source'))
  res$Avg_LISI <- rowMeans(res, na.rm=TRUE)
  output <- mean(res$Avg_LISI)
  output
}

LISI.rPCA <- function(obj) {
  Umap <- as.matrix(obj[["umap"]]@cell.embeddings)
  metadata <- as.data.frame(obj@meta.data)
  metadata <- metadata %>% select(integrated_snn_res.0.8, Source)
  res <- compute_lisi(Umap, metadata, c('integrated_snn_res.0.8', 'Source'))
  res$Avg_LISI <- rowMeans(res, na.rm=TRUE)
  output <- mean(res$Avg_LISI)
  output
}


LISI_harmony_0.8 <- LISI.harm(harmony)
LISI_scanorama_0.8 <- LISI.scan(scanorama)
LISI_rPCA_0.8 <- LISI.rPCA(rPCA)


#########################################################  leiden Function

LISI.harm.leiden <- function(obj) {
  Umap <- as.matrix(obj[["umap"]]@cell.embeddings)
  metadata <- as.data.frame(obj@meta.data)
  metadata <- metadata %>% select(leiden, Source)
  res <- compute_lisi(Umap, metadata, c('leiden', 'Source'))
  res$Avg_LISI <- rowMeans(res, na.rm=TRUE)
  output <- mean(res$Avg_LISI)
  output
}

LISI.scan.leiden <- function(obj) {
  Umap <- as.matrix(obj[["umap"]]@cell.embeddings)
  metadata <- as.data.frame(obj@meta.data)
  metadata <- metadata %>% select(leiden, Source)
  res <- compute_lisi(Umap, metadata, c('leiden', 'Source'))
  res$Avg_LISI <- rowMeans(res, na.rm=TRUE)
  output <- mean(res$Avg_LISI)
  output
}

LISI.rPCA.leiden <- function(obj) {
  Umap <- as.matrix(obj[["umap"]]@cell.embeddings)
  metadata <- as.data.frame(obj@meta.data)
  metadata <- metadata %>% select(leiden, Source)
  res <- compute_lisi(Umap, metadata, c('leiden', 'Source'))
  res$Avg_LISI <- rowMeans(res, na.rm=TRUE)
  output <- mean(res$Avg_LISI)
  output
}


LISI_harmony_leiden <- LISI.harm.leiden(harmony)
LISI_scanorama_leiden <- LISI.scan.leiden(scanorama)
LISI_rPCA_leiden <- LISI.rPCA.leiden(rPCA)


# Presenting results
harmony_avg <- mean(c(LISI_harmony_0.8, LISI_harmony_leiden))
scanorama_avg <- mean(c(LISI_scanorama_0.8, LISI_scanorama_leiden))
rPCA_avg <- mean(c(LISI_rPCA_0.8, LISI_rPCA_leiden))

Louvain <- c(LISI_harmony_0.8,LISI_scanorama_0.8,LISI_rPCA_0.8)
Leiden <- c(LISI_harmony_leiden,LISI_scanorama_leiden,LISI_rPCA_leiden)
Average <- c(harmony_avg,scanorama_avg, rPCA_avg)
  
Integration <- c('harmony', 'scanorama', 'rPCA')


LISI_RESULTS <- data.frame(Integration, Leiden,Louvain, Average)


# Saving results

setwd('C:/Users/Home/Documents/Dissertation/Benchmarking/')

write.csv(LISI_RESULTS, 'LISI_RESULTS')








