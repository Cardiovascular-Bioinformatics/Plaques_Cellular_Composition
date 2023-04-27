# Packages
library(Seurat)
library(tidyverse)
library(ggrepel)
library(ggplot2)


# Loading the carotid dataset 
obj <- readRDS(file = '/Users/judepops/Documents/DISSERTATION/Integrated_Leiden_Louvain0.8/rPCA_Integration_Labelled.rds')

# Visualising the clusters 

Idents(obj)

Carotid_umap <- DimPlot(obj, reduction = "umap", label = TRUE)
Carotid_umap_id <- DimPlot(obj, reduction = "umap", group.by = 'predicted.id', label = TRUE)
Carotid_umap_source <- DimPlot(obj, reduction = "umap", group.by = 'Source', label = TRUE)

ggplot(obj@meta.data, aes(x = Source, fill = predicted.id)) + geom_bar() +
  theme_classic()




# Use the custom color palette in the DimPlot

plot1 <- DimPlot(obj, reduction = "umap", group.by = 'predicted.id', label = FALSE,
                 cols = c('smooth muscle cell' = '#000099', 'pericyte cell' = '#0000CC', 
                          'fibroblast' = '#0000FF', 
                          'macrophage' = '#990000', 't cell' = '#CC0000', 'b cell' = '#FF0000', 
                          'nk cell' = '#FF3333', 'plasma cell' = '#FF6666', 'mast cell' = '#660000', 'artery endothelial cell' = '#006600', 'endothelial cell' = '#009900', 'erythrocyte' = '#00CC00')) + NoAxes() + NoLegend()

LabelClusters(plot1, id = "predicted.id", size = 8, repel = T,  box.padding = 1)



# Making a ggplot

cols = c('smooth muscle cell' = '#000099', 'pericyte cell' = '#0000CC', 
         'fibroblast' = '#0000FF', 
         'macrophage' = '#990000', 't cell' = '#CC0000', 'b cell' = '#FF0000', 
         'nk cell' = '#FF3333', 'plasma cell' = '#FF6666', 'mast cell' = '#660000', 'artery endothelial cell' = '#006600', 'endothelial cell' = '#009900', 'erythrocyte' = '#00CC00')


plot_rgb <- ggplot(obj@meta.data, aes(x = integrated_snn_res.0.8, fill = predicted.id)) +
  geom_bar() +
  scale_fill_manual(values = cols) +  # Add custom colors
  theme_classic() +
  theme(text = element_text(size = 20)) 


cols = c('smooth muscle cell' = '#000099', 'pericyte cell' = '#0000CC', 
         'fibroblast' = '#0000FF', 'mast cell' = '#660000',
         'macrophage' = '#990000', 't cell' = '#CC0000', 'b cell' = '#FF0000', 
         'nk cell' = '#FF3333', 'plasma cell' = '#FF6666', 'artery endothelial cell' = '#006600', 'endothelial cell' = '#009900', 'erythrocyte' = '#00CC00')

# Set up a plot
plot(1, type='n', xlim=c(0, 1), ylim=c(0, 2), axes=FALSE, xlab='', ylab='')

# Add rectangles with the colors and labels from the 'cols' vector
rect(xleft=0.1, ybottom=0.95-seq_along(cols)/20, xright=0.2, ytop=0.96-seq_along(cols)/20, col=cols, border=NA)
text(x=0.3, y=0.95-seq_along(cols)/20, labels=names(cols), adj=0, cex=0.8)

# Add a legend with the same colors and labels
legend('center', fill=cols, legend=names(cols), bty='n')




setwd('/Users/judepops/Documents/DISSERTATION/Figures/')

jpeg(filename="rPCA_TABULA_SAPIENS_labelled.jpg",units="cm",width=35,height=35, res=600)

plot1 <- DimPlot(obj, reduction = "umap", group.by = 'predicted.id', label = FALSE) + NoAxes() + NoLegend()
LabelClusters(plot1, id = "predicted.id", size = 8, repel = T,  box.padding = 1)

dev.off()


setwd('/Users/judepops/Documents/DISSERTATION/Figures/')

jpeg(filename="RGB_rPCA_TABULA_SAPIENS_labelled.jpg",units="cm",width=35,height=35, res=600)

plot1 <- DimPlot(obj, reduction = "umap", group.by = 'predicted.id', label = FALSE,
                  cols = c('smooth muscle cell' = '#000099', 'pericyte cell' = '#0000CC', 
                           'fibroblast' = '#0000FF', 
                           'macrophage' = '#990000', 't cell' = '#CC0000', 'b cell' = '#FF0000', 
                           'nk cell' = '#FF3333', 'plasma cell' = '#FF6666', 'mast cell' = '#660000', 'artery endothelial cell' = '#006600', 'endothelial cell' = '#009900', 'erythrocyte' = '#00CC00')) + NoAxes() + NoLegend()

LabelClusters(plot1, id = "predicted.id", size = 8, repel = T,  box.padding = 1)

dev.off()




setwd('/Users/judepops/Documents/DISSERTATION/Figures/')

jpeg(filename="rPCA_TABULA_SAPIENS_populations_bar.jpg",units="cm",width=40,height=25, res=600)

ggplot(obj@meta.data, aes(x = integrated_snn_res.0.8, fill = predicted.id)) + geom_bar() +
  theme_classic() + theme(text = element_text(size=20)) + NoLegend()

dev.off()



setwd('/Users/judepops/Documents/DISSERTATION/Figures/')

jpeg(filename="RGB_rPCA_TABULA_SAPIENS_populations_bar.jpg",units="cm",width=40,height=25, res=600)

ggplot(obj@meta.data, aes(x = integrated_snn_res.0.8, fill = predicted.id)) +
  geom_bar() +
  scale_fill_manual(values = cols) +  # Add custom colors
  theme_classic() +
  theme(text = element_text(size = 20)) +
  NoLegend()

dev.off()



setwd('/Users/judepops/Documents/DISSERTATION/Figures/')

jpeg(filename="rPCA_25_Louvain_Clusters.jpg",units="cm",width=35,height=35, res=600)

plot2 <- DimPlot(obj, reduction = "umap", group.by = 'integrated_snn_res.0.8', label = FALSE) + NoAxes() + NoLegend() 
LabelClusters(plot2, id = "integrated_snn_res.0.8", size = 8, repel = T,  box.padding = 1)


dev.off()



###################################################################### SUBSETTING


# Subset Clusters that represent smooth muscle cells, pericytes and fibroblasts 
# Using annotations (TABULA SAPIENS)

SMC_clusters <- as.character(c("2", "13", "6", "14", "7", '15'))

# Subset the dataframe for clusters in the smooth muscle cell list

Idents(obj)
obj.smc <- subset(x = obj, idents = SMC_clusters)      


# Visualising new subsetted dataframe

SMC_umap <- DimPlot(obj.smc, reduction = "umap", label = TRUE)


###################################################################### FILTERING


# Removing cells outside UMAP dimensions

SMC_embeddings <- obj.smc[["umap"]]@cell.embeddings 
SMC_embeddings <- as.data.frame(SMC_embeddings)
UMAP_1 <- SMC_embeddings %>% filter(UMAP_1 > 0)
UMAP_2 <- SMC_embeddings %>% filter(UMAP_2 > (-2))

# Creating a dataframe of bad cells to be removed
Bad_Cells <- rbind(UMAP_1, UMAP_2) %>% 
  unique()

# Cells to remove as a vector
toRemove <- as.vector(rownames(Bad_Cells))

## filter them out:
obj.smc.filtered <- obj.smc[,!colnames(obj.smc) %in% toRemove]

#umap
SMC_umap_filtered <- DimPlot(obj.smc.filtered, reduction = "umap")
SMC_umap_filtered_source <-DimPlot(obj.smc.filtered, reduction = "umap", group.by = 'Source', label = TRUE)


setwd('/Users/judepops/Documents/DISSERTATION/Subsetting/')
saveRDS(obj.smc, 'rPCA_SMC.rds')
saveRDS(obj.smc.filtered, 'rPCA_SMC_filtered.rds')



###################################################################### CLUSTERING
obj.smc.filtered <- readRDS('rPCA_SMC_filtered.rds')
DimPlot(obj.smc.filtered, reduction = "umap")

# Re-clustering the subset (calculating PCs, SNN etc.)

# Re-running processing, normalisation and scaling but using ScTransform for higher res
obj.smc.filtered <- ScaleData(object = obj.smc.filtered)
# Running and visualising PCA
obj.smc.filtered <- RunPCA(object = obj.smc.filtered)
ElbowPlot(obj.smc.filtered, ndims = 30)

## Trying different resolutions (0.4-->1.4)

obj.smc.filtered <- FindNeighbors(object = obj.smc.filtered, 
                                   dims = 1:20)

obj.smc.filtered <- FindClusters(object = obj.smc.filtered,
                                  resolution = c(0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4))


setwd('/Users/judepops/Documents/DISSERTATION/Subsetting/')
saveRDS(obj.smc.filtered, 'rPCA_SMC_clustered.rds')

# Exploring the resolutions

obj.smc.clustered <- readRDS('rPCA_SMC_clustered.rds')

obj.smc.clustered@meta.data %>% 
  View()

# Visualising the resolutions
Idents(object = obj.smc.clustered) <- "integrated_snn_res.1.3"
DimPlot(obj.smc.clustered,
                   reduction = "umap",
                   label = TRUE,
                   label.size = 6)


Idents(object = obj.smc.clustered) <- "integrated_snn_res.0.6"
SMC_0.6 <- DimPlot(obj.smc.clustered,
        reduction = "umap",
        label = TRUE,
        label.size = 6)

# Visualising the 0.9 resolution
Idents(object = obj.smc.clustered) <- "integrated_snn_res.0.9"
SMC_0.9 <- DimPlot(obj.smc.clustered,
        reduction = "umap",
        label = TRUE,
        label.size = 6)


# Visualising the plots together

SMC_0.5 | SMC_0.9 | SMC_1.4

# Plotting
DimPlot(obj.smc.clustered, reduction = "umap", label.box = F, label = T) & NoLegend() 


Idents(object = obj.smc.clustered) <- "integrated_snn_res.0.6"
SMC_0.6 <- DimPlot(obj.smc.clustered,
                   reduction = "umap",
                   label = TRUE,
                   label.size = 6)


Idents(object = obj.smc.clustered) <- "integrated_snn_res.0.7"
SMC_0.7 <- DimPlot(obj.smc.clustered,
                   reduction = "umap",
                   label = TRUE,
                   label.size = 6)

Idents(object = obj.smc.clustered) <- "integrated_snn_res.0.8"
SMC_0.8 <- DimPlot(obj.smc.clustered,
                   reduction = "umap",
                   label = TRUE,
                   label.size = 6)


obj.smc.clustered@meta.data

DimPlot(obj.smc.clustered, reduction = "umap", label.box = F, label = T, group.by = 'predicted.id') & NoLegend() 


# Benchmark Data


setwd('/Users/judepops/Documents/DISSERTATION/Subsetting/')

smc.asw <- readRDS('ASW_RESULTS_SMC.rds')
smc.asw <- rename(smc.asw,c("Resolution" = "Clustering.resolution"))

ggplot(data=smc.asw, aes(x=Resolution, y=ASW, group=1)) +
  geom_line() +
  geom_point() + theme_classic()

My_Theme = theme(
  axis.title.x = element_text(size = 16),
  axis.text.x = element_text(size = 14),
  axis.title.y = element_text(size = 16))

DimPlot(obj.smc.clustered,
        reduction = "umap",
        label = TRUE,
        label.size = 6, label.color = 'black') + NoLegend() + NoAxes() 
# Figures


setwd('/Users/judepops/Documents/DISSERTATION/Figures/')
Idents(object = obj.smc.clustered) <- "integrated_snn_res.0.6"

jpeg(filename="SMC_subset_0.6_labelled.jpg",units="cm",width=20,height=20, res=600)

plot_smc_clusters <- DimPlot(obj.smc.clustered,
                   reduction = "umap",
                   label = FALSE,
                   group.by = 'integrated_snn_res.0.6') + NoLegend() + NoAxes() + scale_color_manual(values = my_palette)

LabelClusters(plot_smc_clusters, id = "integrated_snn_res.0.6", size = 8, repel = T,  box.padding = 1)


dev.off()



#Creating blue colour pallette
my_palette <- colorRampPalette(c("#1f78b4", "#164f86", "#08306b", "#8c6bb1", "#88419d", "#6e016b"))(14)

# display the palette
pie(rep(1, 14), col = my_palette, main = "Blue and Purple Palette")




setwd('/Users/judepops/Documents/DISSERTATION/Figures/')

jpeg(filename="SMC_ASW.jpg",units="cm",width=15,height=15, res=600)

ggplot(data=smc.asw, aes(x=Resolution, y=ASW, group=1)) +
  geom_line()+
  geom_point() + theme_classic() + theme(
    axis.title.x = element_text(size = 16),
    axis.text.x = element_text(size = 14),
    axis.title.y = element_text(size = 16),
    axis.text.y = element_text(size = 14))


dev.off()

setwd('/Users/judepops/Documents/DISSERTATION/Figures/')

jpeg(filename="SMC_ASW_BLACK.jpg",units="cm",width=25,height=25, res=600)

ggplot(data=smc.asw, aes(x=Resolution, y=ASW, group=1)) +
  geom_line(color = "black", size = 1.5) +  # set the line color and size
  geom_point() + theme_classic() + theme(
    axis.title.x = element_text(size = 16),
    axis.text.x = element_text(size = 14),
    axis.title.y = element_text(size = 16),
    axis.text.y = element_text(size = 14))


dev.off()


my_palette <- rainbow(14)

setwd('/Users/judepops/Documents/DISSERTATION/Figures/')

jpeg(filename="Rainbow_SMC_subset_0.6_labelled.jpg",units="cm",width=20,height=20, res=600)

plot_smc_clusters <- DimPlot(obj.smc.clustered,
                             reduction = "umap",
                             label = FALSE,
                             group.by = 'integrated_snn_res.0.6') + NoAxes() + scale_color_manual(values = my_palette)

LabelClusters(plot_smc_clusters, id = "integrated_snn_res.0.6", size = 8, repel = T,  box.padding = 1)


dev.off()



my_palette <- rainbow(14)

setwd('/Users/judepops/Documents/DISSERTATION/Figures/')

jpeg(filename="Rainbow_SMC_subset_0.6_labelled.jpg",units="cm",width=20,height=20, res=600)

plot_smc_clusters <- DimPlot(obj.smc.clustered,
                             reduction = "umap",
                             label = FALSE,
                             group.by = 'integrated_snn_res.0.6') + NoAxes() + scale_color_manual(values = my_palette)

LabelClusters(plot_smc_clusters, id = "integrated_snn_res.0.6", size = 8, repel = T,  box.padding = 1)


dev.off()








setwd('/Users/judepops/Documents/DISSERTATION/Figures/')

jpeg(filename="BLUE_SMC_subset_0.6_labelled.jpg",units="cm",width=20,height=20, res=600)

plot_subset <- DimPlot(obj.smc.clustered, reduction = "umap", group.by = 'predicted.id', label = FALSE,
                       cols = c('smooth muscle cell' = '#000099', 'pericyte cell' = '#0000CC', 
                                'fibroblast' = '#0000FF')) + NoLegend() + NoAxes()
plot_subset

dev.off()





setwd('/Users/judepops/Documents/DISSERTATION/Figures/')

jpeg(filename="RGB_rPCA_TABULA_SAPIENS.jpg",units="cm",width=35,height=35, res=600)

plot1 <- DimPlot(obj, reduction = "umap", group.by = 'predicted.id', label = FALSE,
                 cols = c('smooth muscle cell' = '#000099', 'pericyte cell' = '#0000CC', 
                          'fibroblast' = '#0000FF', 
                          'macrophage' = '#990000', 't cell' = '#CC0000', 'b cell' = '#FF0000', 
                          'nk cell' = '#FF3333', 'plasma cell' = '#FF6666', 'mast cell' = '#660000', 'artery endothelial cell' = '#006600', 'endothelial cell' = '#009900', 'erythrocyte' = '#00CC00')) + NoAxes() + NoLegend()

plot1
dev.off()

setwd('/Users/judepops/Documents/DISSERTATION/Figures/')

jpeg(filename="Cell_Populations_Per_Source_rPCA.jpg",units="cm",width=20,height=25, res=600)


cols = c('smooth muscle cell' = '#000099', 'pericyte cell' = '#0000CC', 
         'fibroblast' = '#0000FF', 
         'macrophage' = '#990000', 't cell' = '#CC0000', 'b cell' = '#FF0000', 
         'nk cell' = '#FF3333', 'plasma cell' = '#FF6666', 'mast cell' = '#660000', 'artery endothelial cell' = '#006600', 'endothelial cell' = '#009900', 'erythrocyte' = '#00CC00')


ggplot(obj@meta.data, aes(x = Source, fill = predicted.id)) +
  geom_bar() +
  scale_fill_manual(values = cols) +  # Add custom colors
  theme_classic() +
  theme(text = element_text(size = 20)) +
  NoLegend()


dev.off()



