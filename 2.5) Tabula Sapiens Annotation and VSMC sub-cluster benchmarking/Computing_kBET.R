#Running Silhouette score metric
library(Seurat)
library(lisi)
library(dplyr)
library(tidyr)
library(ggplot2)
library(magrittr)
library(stringr)
library(tidyverse)
#Kbet

library(kBET)


#data

setwd('C:/Users/Home/Documents/Dissertation/FInal_Integrated/')

harmony <- readRDS('Harmony_Integrated_Leiden_Louvain.rds')
scanorama <- readRDS('Panorama_Integrated_Leiden_Louvain.rds')
rPCA <- readRDS('rPCA_Integrated_Leiden_Louvain.rds')


# Preparing data
matrix <- harmony@assays$RNA@scale.data
data <- t(matrix)
# Preparing batch
batch <- harmony@meta.data$Source
#kBET
subset_size <- 0.1 #subsample to 50% of the data
subset_id <- sample.int(n = length(batch), size = floor(subset_size * length(batch)), replace=FALSE)
batch.estimate <- kBET(data[subset_id,], batch[subset_id])
#batch.estimate <- kBET(data, batch)
batch.estimate$stats$kBET.observed 


# Making a function

kBET.harmony <- function(obj) {
  matrix <- obj@assays$RNA@scale.data
  data <- t(matrix)
  # Preparing batch
  batch <- obj@meta.data$Source
  #kBET
  subset_size <- 0.5
  subset_id <- sample.int(n = length(batch), size = floor(subset_size * length(batch)), replace=FALSE)
  batch.estimate <- kBET(data[subset_id,], batch[subset_id])
  batch.estimate
}

kBET.scanorama <- function(obj) {
  matrix <- obj@assays$pano@scale.data
  data <- t(matrix) #transposing the matrix so that 
  # Preparing batch
  batch <- obj@meta.data$Source
  #kBET
  subset_size <- 0.5
  subset_id <- sample.int(n = length(batch), size = floor(subset_size * length(batch)), replace=FALSE)
  batch.estimate <- kBET(data[subset_id,], batch[subset_id])
  batch.estimate
}

kBET.rPCA <- function(obj) {
  matrix <- rPCA@assays$integrated@scale.data
  data <- t(matrix)
  # Preparing batch
  batch <- obj@meta.data$Source
  #kBET
  subset_size <- 0.5
  subset_id <- sample.int(n = length(batch), size = floor(subset_size * length(batch)), replace=FALSE)
  batch.estimate <- kBET(data[subset_id,], batch[subset_id])
  batch.estimate
}

# Calculating kbet

harmony_kBET <- kBET.harmony(harmony)
scanorama_kBET <- kBET.scanorama(scanorama)
rPCA_kBET <- kBET.rPCA(rPCA)


harmony_kBET_Avg <- mean(harmony_kBET$stats$kBET.observed)
scanorama_kBET_Avg <- mean(harmony_kBET$stats$kBET.observed)
rPCA_kBET_Avg <- mean(rPCA_kBET$stats$kBET.observed)

  

# Presenting results

kBET.observed <- c(harmony_kBET_Avg,scanorama_kBET_Avg, rPCA_kBET_Avg)

Integration <- c('harmony', 'scanorama', 'rPCA')


KBET_RESULTS <- data.frame(Integration, kBET.observed)


# Saving results

setwd('C:/Users/Home/Documents/Dissertation/Benchmarking/')

write.csv(KBET_RESULTS, 'KBET_RESULTS')














matrix <- scanorama@assays$pano@data
data <- t(matrix)




harmony_kBET$stats$kBET.observed 
rPCA_kBET$stats$kBET.observed 







matrix <- rPCA@assays$integrated@scale.data
data <- t(matrix)








