# Packages
library(Seurat)
library(tidyverse)
library(ggplot2)
library(RColorBrewer)
library(cowplot)
install.packages("RColorBrewer")

# Data
setwd('/Users/judepops/Documents/DISSERTATION/Subsetting/')
obj.smc.clustered <- readRDS('rPCA_SMC_clustered.rds')

# Setting idents to the best clustering parameter

Idents(obj.smc.clustered)
Idents(object = obj.smc.clustered) <- "integrated_snn_res.0.6"

DimPlot(obj.smc.clustered, reduction = 'umap', group.by = 'integrated_snn_res.0.6')

# Setting the key features for the FeaturePlot

smc.markers <- as.vector(c('TAGLN2', 'ACTA2', 'CNN1', 'MGP', 'VCAN', 'KLF4'))
x <- FeaturePlot(object = obj.smc.clustered, features = smc.markers, keep.scale = 'all')  

FeaturePlot(object = obj.smc.clustered, features = 'KLF4', keep.scale = 'all')

# For loop to plot the feature plots


myfeatures <-  as.vector(c('TAGLN2', 'ACTA2', 'CNN1', 'MGP', 'VCAN', 'ACAN'))
plot_list <- list()
for (i in myfeatures) {
  plot_list[[i]] <- FeaturePlot(obj.smc.clustered, reduction = "umap", features = i,
                                ncol = 3, order = T, keep.scale = 'all') + NoAxes() + NoGrid() & scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu"))) 
}
plot_grid(ncol = 3, plotlist = plot_list)






setwd('/Users/judepops/Documents/DISSERTATION/Figures/')

jpeg(filename="Loop_SMC_featureplots.jpg",units="cm",width=60,height=20, res=600)

plot_grid(ncol = 3, plotlist = plot_list)


dev.off()


setwd('/Users/judepops/Documents/DISSERTATION/Figures/')

jpeg(filename="KLF4_SMC_featureplots_2.jpg",units="cm",width=20,height=30, res=600)

smc.fplot

dev.off()



