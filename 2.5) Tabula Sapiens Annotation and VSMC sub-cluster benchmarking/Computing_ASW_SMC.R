# Computing silhouette metric R

# libraries
library(dplyr)
library(tidyr)
library(ggplot2)
library(magrittr)
library(stringr)
library(Seurat)
library(tidyverse)
library(cluster, quietly = TRUE)
library(ggplot2)

# Dataset

setwd('C:/Users/Home/Documents/Dissertation/Subset')
obj.smc <- readRDS('rPCA_SMC_clustered.rds')



# Performing sillhouette on each one

# 0.3
dims <- 1:20
reduction <- 'pca'
dist.matrix <- dist(x = Embeddings(object = obj.smc[[reduction]])[, dims])
clusters <- obj.smc$integrated_snn_res.0.3 # Changed for each one
sil <- silhouette(x = as.numeric(x = as.factor(x = clusters)), dist = dist.matrix)
obj.smc$sil.0.3 <- sil[, 3]  


# 0.4
dims <- 1:20
reduction <- 'pca'
dist.matrix <- dist(x = Embeddings(object = obj.smc[[reduction]])[, dims])
clusters <- obj.smc$integrated_snn_res.0.4 # Changed for each one
sil <- silhouette(x = as.numeric(x = as.factor(x = clusters)), dist = dist.matrix)
obj.smc$sil.0.4 <- sil[, 3]


# 0.5
dims <- 1:20
reduction <- 'pca'
dist.matrix <- dist(x = Embeddings(object = obj.smc[[reduction]])[, dims])
clusters <- obj.smc$integrated_snn_res.0.5 # Changed for each one
sil <- silhouette(x = as.numeric(x = as.factor(x = clusters)), dist = dist.matrix)
obj.smc$sil.0.5 <- sil[, 3]

# 0.6
dims <- 1:20
reduction <- 'pca'
dist.matrix <- dist(x = Embeddings(object = obj.smc[[reduction]])[, dims])
clusters <- obj.smc$integrated_snn_res.0.6 # Changed for each one
sil <- silhouette(x = as.numeric(x = as.factor(x = clusters)), dist = dist.matrix)
obj.smc$sil.0.6 <- sil[, 3]

# 0.7
dims <- 1:20
reduction <- 'pca'
dist.matrix <- dist(x = Embeddings(object = obj.smc[[reduction]])[, dims])
clusters <- obj.smc$integrated_snn_res.0.7 # Changed for each one
sil <- silhouette(x = as.numeric(x = as.factor(x = clusters)), dist = dist.matrix)
obj.smc$sil.0.7 <- sil[, 3]

# 0.8
dims <- 1:20
reduction <- 'pca'
dist.matrix <- dist(x = Embeddings(object = obj.smc[[reduction]])[, dims])
clusters <- obj.smc$integrated_snn_res.0.8 # Changed for each one
sil <- silhouette(x = as.numeric(x = as.factor(x = clusters)), dist = dist.matrix)
obj.smc$sil.0.8 <- sil[, 3]

# 0.9
dims <- 1:20
reduction <- 'pca'
dist.matrix <- dist(x = Embeddings(object = obj.smc[[reduction]])[, dims])
clusters <- obj.smc$integrated_snn_res.0.9 # Changed for each one
sil <- silhouette(x = as.numeric(x = as.factor(x = clusters)), dist = dist.matrix)
obj.smc$sil.0.9 <- sil[, 3]

# 1.0
dims <- 1:20
reduction <- 'pca'
dist.matrix <- dist(x = Embeddings(object = obj.smc[[reduction]])[, dims])
clusters <- obj.smc$integrated_snn_res.1 # Changed for each one
sil <- silhouette(x = as.numeric(x = as.factor(x = clusters)), dist = dist.matrix)
obj.smc$sil.1.0 <- sil[, 3]

# 1.1
dims <- 1:20
reduction <- 'pca'
dist.matrix <- dist(x = Embeddings(object = obj.smc[[reduction]])[, dims])
clusters <- obj.smc$integrated_snn_res.1.1 # Changed for each one
sil <- silhouette(x = as.numeric(x = as.factor(x = clusters)), dist = dist.matrix)
obj.smc$sil.1.1 <- sil[, 3]

# 1.2
dims <- 1:20
reduction <- 'pca'
dist.matrix <- dist(x = Embeddings(object = obj.smc[[reduction]])[, dims])
clusters <- obj.smc$integrated_snn_res.1.2 # Changed for each one
sil <- silhouette(x = as.numeric(x = as.factor(x = clusters)), dist = dist.matrix)
obj.smc$sil.1.2 <- sil[, 3]

# 1.3
dims <- 1:20
reduction <- 'pca'
dist.matrix <- dist(x = Embeddings(object = obj.smc[[reduction]])[, dims])
clusters <- obj.smc$integrated_snn_res.1.3 # Changed for each one
sil <- silhouette(x = as.numeric(x = as.factor(x = clusters)), dist = dist.matrix)
obj.smc$sil.1.3 <- sil[, 3]

# 1.4
dims <- 1:20
reduction <- 'pca'
dist.matrix <- dist(x = Embeddings(object = obj.smc[[reduction]])[, dims])
clusters <- obj.smc$integrated_snn_res.1.4 # Changed for each one
sil <- silhouette(x = as.numeric(x = as.factor(x = clusters)), dist = dist.matrix)
obj.smc$sil.1.4 <- sil[, 3]




# Viewing scores
obj.smc@meta.data

obj.smc.0.3 <- mean(obj.smc$sil.0.3)
obj.smc.0.4 <- mean(obj.smc$sil.0.4)
obj.smc.0.5 <- mean(obj.smc$sil.0.5)
obj.smc.0.6 <- mean(obj.smc$sil.0.6)
obj.smc.0.7 <- mean(obj.smc$sil.0.7)
obj.smc.0.8 <- mean(obj.smc$sil.0.8)
obj.smc.0.9 <- mean(obj.smc$sil.0.9)
obj.smc.1.0 <- mean(obj.smc$sil.1.0)
obj.smc.1.1 <- mean(obj.smc$sil.1.1)
obj.smc.1.2 <- mean(obj.smc$sil.1.2)
obj.smc.1.3 <- mean(obj.smc$sil.1.3)
obj.smc.1.4 <- mean(obj.smc$sil.1.4)


# Creating a dataframe with results


ASW <- c(obj.smc.0.3,obj.smc.0.4, obj.smc.0.5,obj.smc.0.6,obj.smc.0.7,
         obj.smc.0.8,obj.smc.0.9,obj.smc.1.0,obj.smc.1.1,obj.smc.1.2,obj.smc.1.3,obj.smc.1.4)

Resolution <- c('0.3','0.4', '0.5','0.6','0.7',
                           '0.8','0.9','1.0','1.1','1.2','1.3','1.4')


ASW_RESULTS_SMC <- data.frame(Clustering.resolution, ASW)
setwd('C:/Users/Home/Documents/Dissertation/Benchmarking/')

saveRDS(ASW_RESULTS_SMC, 'ASW_RESULTS_SMC.rds')


ggplot(data=ASW_RESULTS_SMC, aes(x=Resolution, y=ASW, group=1)) +
  geom_line()+
  geom_point() + theme_classic()








