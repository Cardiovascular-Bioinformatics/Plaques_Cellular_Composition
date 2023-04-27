# Computing silhouette metric R

#libraries
library(dplyr)
library(tidyr)
library(ggplot2)
library(magrittr)
library(stringr)
library(Seurat)
library(tidyverse)
library(cluster, quietly = TRUE)


#Datasets 
setwd('C:/Users/Home/Documents/Dissertation/FInal_Integrated/')

harmony <- readRDS('Harmony_Integrated_Leiden_Louvain.rds')
scanorama <- readRDS('Panorama_Integrated_Leiden_Louvain.rds')
rPCA <- readRDS('rPCA_Integrated_Leiden_Louvain.rds')



# silhouette metric
dims <- 1:10
reduction <- 'pca'
dist.matrix <- dist(x = Embeddings(object = harmony[[reduction]])[, dims])
clusters <- harmony$RNA_snn_res.0.8
sil.harmony <- silhouette(x = as.numeric(x = as.factor(x = clusters)), dist = dist.matrix)
harmony$sil <- sil.harmony[, 3]  ##THIS IS THE CALCULATED SILHOUETTE SCORE OF THE CLUSTERING

dims <- 1:10
reduction <- 'pca'
dist.matrix <- dist(x = Embeddings(object = scanorama[[reduction]])[, dims])
clusters <- scanorama$pano_snn_res.0.8
sil.scanorama <- silhouette(x = as.numeric(x = as.factor(x = clusters)), dist = dist.matrix)
scanorama$sil <- sil.scanorama[, 3]  ##THIS IS THE CALCULATED SILHOUETTE SCORE OF THE CLUSTERING


dims <- 1:10
reduction <- 'pca'
dist.matrix <- dist(x = Embeddings(object = rPCA[[reduction]])[, dims])
clusters <- rPCA$integrated_snn_res.0.8
sil.rPCA <- silhouette(x = as.numeric(x = as.factor(x = clusters)), dist = dist.matrix)
rPCA$sil <- sil.rPCA[, 3]  ##THIS IS THE CALCULATED SILHOUETTE SCORE OF THE CLUSTERING

setwd('C:/Users/Home/Documents/Dissertation/Benchmarking/')

write.csv(harmony@meta.data, 'harmony_sil.csv')
write.csv(scanorama@meta.data, 'scanorama_sil.csv')
write.csv(rPCA@meta.data, 'rPCA_sil.csv')

# Calculating ASW (average)

ASW_harmony <- mean(harmony@meta.data$sil)
ASW_scanorama <- mean(scanorama@meta.data$sil)
ASW_rPCA <- mean(rPCA@meta.data$sil)


# Creating a dataframe with results


ASW <- c(ASW_harmony,ASW_scanorama, ASW_rPCA)

Integration <- c('harmony', 'scanorama', 'rPCA')


ASW_RESULTS <- data.frame(Integration, ASW)

# Saving file

setwd('C:/Users/Home/Documents/Dissertation/Benchmarking/')

saveRDS(ASW_RESULTS, 'ASW_RESULTS.rds')


