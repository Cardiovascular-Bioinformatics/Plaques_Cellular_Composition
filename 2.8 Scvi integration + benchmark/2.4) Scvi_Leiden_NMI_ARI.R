
# Importing csv files
library(tidyverse)
library(Seurat)

# Data

fernandes <- readRDS('/Users/judepops/Documents/DISSERTATION/Source_Leiden_Louvain0.8/Fernandes_Leiden_Louvain.rds')
alsaigh <- readRDS('/Users/judepops/Documents/DISSERTATION/Source_Leiden_Louvain0.8/Alsaigh_Leiden_Louvain.rds')
pan <- readRDS('/Users/judepops/Documents/DISSERTATION/Source_Leiden_Louvain0.8/Pan_Leiden_Louvain.rds')


# Extracting 0.8 clustering (default)

harmony <- readRDS('/Users/judepops/Documents/DISSERTATION/Integrated_Leiden_Louvain0.8/Harmony_Integrated_Leiden_Louvain.rds')
rPCA <- readRDS('/Users/judepops/Documents/DISSERTATION/Integrated_Leiden_Louvain0.8/rPCA_Integrated_Leiden_Louvain.rds')
scanorama <- readRDS('/Users/judepops/Documents/DISSERTATION/Integrated_Leiden_Louvain0.8/Panorama_Integrated_Leiden_Louvain.rds')


# Extracting metadata

fernandes <- fernandes@meta.data
alsaigh <- alsaigh@meta.data
pan <- pan@meta.data
harmony <- harmony@meta.data
rPCA <- rPCA@meta.data
scanorama <- scanorama@meta.data


# Function to create louvain cluster vectors for source files

source.clusters <- function(obj) {
  obj <- obj %>% rownames_to_column(var = 'X')
  obj <- obj[order(obj$X),]
  obj <- obj %>% select(leiden)
  obj
}

fernandes <- source.clusters(fernandes)
pan <- source.clusters(pan)
alsaigh <- source.clusters(alsaigh)

#Column to rownames for each one
harmony <- harmony %>% rownames_to_column(var = 'X')
scanorama <- scanorama %>% rownames_to_column(var = 'X')
rPCA <- rPCA %>% rownames_to_column(var = 'X')
scvi <- scvi %>% rownames_to_column(var = 'X')

## Subsetting the dataframe for each experiment source

#Harmony

#Subsetting
harmony_fernandes <- harmony[harmony$Source == 'fernandes',]
harmony_pan <- harmony[harmony$Source == 'Pan',]
harmony_alsaigh <- harmony[harmony$Source == 'alsaigh',]
#ordering
harmony_fernandes <- harmony_fernandes[order(harmony_fernandes$X),]
harmony_pan <- harmony_pan[order(harmony_pan$X),]
harmony_alsaigh <- harmony_alsaigh[order(harmony_alsaigh$X),]
#Selecting
harmony_fernandes <- harmony_fernandes %>% select(leiden)
harmony_pan <- harmony_pan %>% select(leiden)
harmony_alsaigh<- harmony_alsaigh %>% select(leiden)


# Scanorama

#Subsetting
scanorama_fernandes <- scanorama[scanorama$Source == 'fernandes',]
scanorama_pan <- scanorama[scanorama$Source == 'Pan',]
scanorama_alsaigh <- scanorama[scanorama$Source == 'alsaigh',]
#Ordering
scanorama_fernandes <- scanorama_fernandes[order(scanorama_fernandes$X),]
scanorama_pan <- scanorama_pan[order(scanorama_pan$X),]
scanorama_alsaigh <- scanorama_alsaigh[order(scanorama_alsaigh$X),]
#Selecting
scanorama_fernandes <- scanorama_fernandes %>% select(leiden)
scanorama_pan <- scanorama_pan %>% select(leiden)
scanorama_alsaigh<- scanorama_alsaigh %>% select(leiden)

#rPCA

#Subsetting
rPCA_fernandes <- rPCA[rPCA$Source == 'fernandes',]
rPCA_pan <- rPCA[rPCA$Source == 'Pan',]
rPCA_alsaigh <- rPCA[rPCA$Source == 'alsaigh',]
#Ordering
rPCA_fernandes <- rPCA_fernandes[order(rPCA_fernandes$X),]
rPCA_pan <- rPCA_pan[order(rPCA_pan$X),]
rPCA_alsaigh <- rPCA_alsaigh[order(rPCA_alsaigh$X),]
#Selecting
rPCA_fernandes <- rPCA_fernandes %>% select(leiden)
rPCA_pan <- rPCA_pan %>% select(leiden)
rPCA_alsaigh<- rPCA_alsaigh %>% select(leiden)




## Integrated


##############################################################################################


# Computing NMI
library(aricode)


### 1) HARMONY

## fernandes
harmony_fernandes_vector <- as.vector(harmony_fernandes['leiden'])
fernandes_vector <- as.vector(fernandes['leiden'])

harmony_fernandes_vector <- unlist(harmony_fernandes_vector) 
fernandes_vector <- unlist(fernandes_vector)

harmony_NMI_fernandes <- print(aricode::NMI(harmony_fernandes_vector , fernandes_vector))   
harmony_ARI_fernandes <- print(aricode::ARI(harmony_fernandes_vector , fernandes_vector))  

## alsaigh

harmony_alsaigh_vector <- as.vector(harmony_alsaigh['leiden'])
alsaigh_vector <- as.vector(alsaigh['leiden'])

harmony_alsaigh_vector <- unlist(harmony_alsaigh_vector) 
alsaigh_vector <- unlist(alsaigh_vector)


harmony_NMI_alsaigh <- print(aricode::NMI(harmony_alsaigh_vector , alsaigh_vector))   
harmony_ARI_alsaigh <- print(aricode::ARI(harmony_alsaigh_vector , alsaigh_vector))  

## pan

harmony_pan_vector <- as.vector(harmony_pan['leiden'])
pan_vector <- as.vector(pan['leiden'])

harmony_pan_vector <- unlist(harmony_pan_vector) 
pan_vector <- unlist(pan_vector)


harmony_NMI_pan <- print(aricode::NMI(harmony_pan_vector , pan_vector))  
harmony_ARI_pan <- print(aricode::ARI(harmony_pan_vector , pan_vector))  





### 2) Scanorama

## fernandes
scanorama_fernandes_vector <- as.vector(scanorama_fernandes['leiden'])
fernandes_vector <- as.vector(fernandes['leiden'])
scanorama_fernandes_vector <- unlist(scanorama_fernandes_vector) 
fernandes_vector <- unlist(fernandes_vector)

scanorama_NMI_fernandes <- print(aricode::NMI(scanorama_fernandes_vector , fernandes_vector))   
scanorama_ARI_fernandes <- print(aricode::ARI(scanorama_fernandes_vector , fernandes_vector))  

## alsaigh

scanorama_alsaigh_vector <- as.vector(scanorama_alsaigh['leiden'])
alsaigh_vector <- as.vector(alsaigh['leiden'])
scanorama_alsaigh_vector <- unlist(scanorama_alsaigh_vector) 
alsaigh_vector <- unlist(alsaigh_vector)


scanorama_NMI_alsaigh <- print(aricode::NMI(scanorama_alsaigh_vector , alsaigh_vector))   
scanorama_ARI_alsaigh <- print(aricode::ARI(scanorama_alsaigh_vector , alsaigh_vector))  

## pan

scanorama_pan_vector <- as.vector(scanorama_pan['leiden'])
pan_vector <- as.vector(pan['leiden'])
scanorama_pan_vector <- unlist(scanorama_pan_vector) 
pan_vector <- unlist(pan_vector)


scanorama_NMI_pan <- print(aricode::NMI(scanorama_pan_vector , pan_vector))  
scanorama_ARI_pan <- print(aricode::ARI(scanorama_pan_vector , pan_vector))  


### rPCA

## fernandes
rPCA_fernandes_vector <- as.vector(rPCA_fernandes['leiden'])
fernandes_vector <- as.vector(fernandes['leiden'])
rPCA_fernandes_vector <- unlist(rPCA_fernandes_vector) 
fernandes_vector <- unlist(fernandes_vector)

rPCA_NMI_fernandes <- aricode::NMI(rPCA_fernandes_vector , fernandes_vector)  
rPCA_ARI_fernandes <- aricode::ARI(rPCA_fernandes_vector , fernandes_vector) 

## alsaigh

rPCA_alsaigh_vector <- as.vector(rPCA_alsaigh['leiden'])
alsaigh_vector <- as.vector(alsaigh['leiden'])
rPCA_alsaigh <- unlist(rPCA_alsaigh) 
alsaigh_vector <- unlist(alsaigh_vector)


rPCA_NMI_alsaigh <- aricode::NMI(rPCA_alsaigh , alsaigh_vector) 
rPCA_ARI_alsaigh <- aricode::ARI(rPCA_alsaigh , alsaigh_vector)

## pan

rPCA_pan_vector <- as.vector(rPCA_pan['leiden'])
pan_vector <- as.vector(pan['leiden'])
rPCA_pan_vector <- unlist(rPCA_pan_vector) 
pan_vector <- unlist(pan_vector)

rPCA_NMI_pan <- aricode::NMI(rPCA_pan_vector , pan_vector)
rPCA_ARI_pan <- aricode::ARI(rPCA_pan_vector , pan_vector)


##############################################################################################


## Averages ARI

average_harmony_NMI <- mean(c(harmony_NMI_fernandes, harmony_NMI_alsaigh, harmony_NMI_pan))
average_scanorama_NMI <- mean(c(scanorama_NMI_fernandes, scanorama_NMI_alsaigh, scanorama_NMI_pan))
average_rPCA_NMI <- mean(c(rPCA_NMI_fernandes,rPCA_NMI_alsaigh,rPCA_NMI_pan))



# Making a dataframe from results

# Create a, b, c, d variables
Source <- c('fernandez', 'alsaigh', 'pan', 'average')
Harmony_NMI <- c(harmony_NMI_fernandes, harmony_NMI_alsaigh, harmony_NMI_pan, average_harmony_NMI)
rPCA_NMI <- c(rPCA_NMI_fernandes,rPCA_NMI_alsaigh,rPCA_NMI_pan, average_rPCA_NMI)
Scanorama_NMI <- c(scanorama_NMI_fernandes, scanorama_NMI_alsaigh, scanorama_NMI_pan, average_scanorama_NMI)


# Join the variables to create a data frame
NMI_RESULTS_LEIDEN <- data.frame(Source,Harmony_NMI,rPCA_NMI,Scanorama_NMI)




##############################################################################################


## Averages ARI
average_harmony_ARI <- mean(c(harmony_ARI_fernandes, harmony_ARI_alsaigh, harmony_ARI_pan))
average_scanorama_ARI <- mean(c(scanorama_ARI_fernandes, scanorama_ARI_alsaigh, scanorama_ARI_pan))
average_rPCA_ARI <- mean(c(rPCA_ARI_fernandes,rPCA_ARI_alsaigh,rPCA_ARI_pan))


# Making a dataframe from results

# Create variables
Source <- c('fernandez', 'alsaigh', 'pan', 'average')
Harmony_ARI <- c(harmony_ARI_fernandes, harmony_ARI_alsaigh, harmony_ARI_pan, average_harmony_ARI)
rPCA_ARI <- c(rPCA_ARI_fernandes,rPCA_ARI_alsaigh,rPCA_ARI_pan, average_rPCA_ARI)
Scanorama_ARI <- c(scanorama_ARI_fernandes, scanorama_ARI_alsaigh, scanorama_ARI_pan, average_scanorama_ARI)


# Join the variables to create a data frame
ARI_RESULTS_LEIDEN <- data.frame(Source,Harmony_ARI,rPCA_ARI,Scanorama_ARI)



# setwd('/Users/judepops/Documents/DISSERTATION/Benchmarking/')
# 
# write_csv(df1, 'NMI_results.csv')
# write_csv(df2, 'ARI_results.csv')


