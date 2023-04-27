# libraries
library(Seurat)
library(tidyverse)
# Integrating with tabula sapiens
setwd('C:/Users/Home/Documents/Dissertation/Carotid_Annotation/')
# Loading data

obj <- readRDS('TabulaSapiensVasculature.rds')

# Setting wd

setwd('C:/Users/Home/Documents/Dissertation/FInal_Integrated/')

Integrated <- readRDS('rPCA_Integrated_Leiden_Louvain.rds')

# Finding anchors to annotate

tabula.anchors <- FindTransferAnchors(reference = obj, query = Integrated,
                                      dims = 1:30, reference.reduction = "pca")


# Calculating Predictions

predictions <- TransferData(anchorset = tabula.anchors, refdata = obj$cell_ontology_class,
                            dims = 1:30)

# Adding the metadata

Integrated <- AddMetaData(Integrated, metadata = predictions)

# Running UMAP

DimPlot(Integrated, reduction = 'umap', label = TRUE)
Integrated@meta.data

Integrated@meta.data %>% group_by(RNA_snn_res.0.5)


DimPlot(Integrated, reduction = "umap", group.by = "predicted.id", label = TRUE, label.size = 3, repel = TRUE) + NoLegend()

# Saving the file 
saveRDS(Integrated, file = "rPCA_Integration_Labelled.rds")







SMC_clusters <- as.character(c("2", "13", "6", "14", "7", '15'))

# Subset the dataframe for clusters in the smooth muscle cell list

Idents(Integrated)
obj.smc <- subset(x = Integrated, idents = SMC_clusters) 


# Normalising

obj.smc <- NormalizeData(object = obj.smc)
