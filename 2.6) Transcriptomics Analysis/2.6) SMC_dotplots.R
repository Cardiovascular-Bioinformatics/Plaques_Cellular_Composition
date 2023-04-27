# Packages
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(cowplot)
  library(ggplot2)
  library(pheatmap)
  library(enrichR)
  library(rafalib)
  library(Matrix)
  library(edgeR)
  library(MAST)
  library(tidyverse)
})

# Data
setwd('/Users/judepops/Documents/DISSERTATION/Subsetting/')
obj.smc.clustered <- readRDS('rPCA_SMC_clustered.rds')

# Set the identity as louvain with resolution 0.5
sel.clust = "integrated_snn_res.0.6"

obj.smc.clustered <- SetIdent(obj.smc.clustered, value = sel.clust)
table(obj.smc.clustered@active.ident)


# plot this clustering
plot_grid(ncol = 3, DimPlot(obj.smc.clustered, label = T) + NoAxes(), 
          DimPlot(obj.smc.clustered, group.by = "Source") +
            NoAxes(), DimPlot(obj.smc.clustered, group.by = "predicted.id") + NoAxes())



# Finding marker genes

markers_genes <- FindAllMarkers(obj.smc.clustered, log2FC.threshold = 0.2, 
                                test.use = "wilcox",
                                min.pct = 0.1, min.diff.pct = 0.2, only.pos = TRUE,
                                assay = "RNA")

setwd('/Users/judepops/Documents/DISSERTATION/SMCFeatures/')
saveRDS(markers_genes, 'SMC_markers.rds')
markers_genes <- readRDS('SMC_markers.rds')
# Selecting the top 5 marker genes

top5 <- markers_genes %>% group_by(cluster) %>% slice_min(order_by = p_val_adj, n = 5)
top3 <- markers_genes %>% group_by(cluster) %>% slice_min(order_by = p_val_adj, n = 3)


top5_random <- top5 %>% 
  group_by(cluster) %>% slice_sample(n=5)

top3_random <- top3 %>% 
  group_by(cluster) %>% slice_sample(n=3)



setwd('/Users/judepops/Documents/DISSERTATION/SMCFeatures/')
write.csv(top5, 'top5_per_SMC_cluster.csv')
write.csv(top5_random, 'random_top5_per_SMC_cluster.csv')

# create a scale.data slot for the selected genes 
obj.smc.clustered <- ScaleData(obj.smc.clustered, features = as.character(unique(top5_random$gene)), assay = "RNA")


#Visualising with a dotplot
smc.top5.dotplot <- DotPlot(obj.smc.clustered, features = rev(as.character(unique(top5_random$gene))), group.by = sel.clust,
        assay = "RNA") + coord_flip() + theme(axis.text=element_text(size=6)) + theme_classic()


#Top 3 dotplot for each cluster
obj.smc.clustered <- ScaleData(obj.smc.clustered, features = as.character(unique(top3_random$gene)), assay = "RNA")
# plotting heatmap for the selected genes which have been scaled

smc.top3.heatmap <- DoHeatmap(obj.smc.clustered, features = as.character(unique(top3_random$gene)), 
                              group.by = sel.clust,assay = "RNA") + theme(axis.text=element_text(size=10))




top5_processed <- read.csv('/Users/judepops/Documents/DISSERTATION/SMCFeatures/Copy of top5_per_SMC_cluster_processed.csv')
top5_processed %>% column_to_rownames(X)
top5_processed <- top5_processed[,-1]

smc.top5.dotplot.processed <- DotPlot(obj.smc.clustered, features = rev(as.character(unique(top5_processed$gene))), group.by = sel.clust,
                            assay = "RNA") + coord_flip() + theme(axis.text=element_text(size=6)) + theme_classic()



# Printing image 
setwd('/Users/judepops/Documents/DISSERTATION/Figures/')

jpeg(filename="SMC_top5.jpg",units="cm",width=30,height=40, res=600) 

smc.top5.dotplot

dev.off()


setwd('/Users/judepops/Documents/DISSERTATION/Figures/')

jpeg(filename="SMC_top5_PROCESSED.jpg",units="cm",width=20,height=25, res=600) 

smc.top5.dotplot.processed

dev.off()



setwd('/Users/judepops/Documents/DISSERTATION/Figures/')

jpeg(filename="SMC_top3_Heatmap.jpg",units="cm",width=50,height=30, res=600)

smc.top3.heatmap

dev.off()





top5_processed <- read.csv('/Users/judepops/Documents/DISSERTATION/SMCFeatures/Copy of top5_per_SMC_cluster_processed.csv')
top5_processed %>% column_to_rownames(X)
top5_processed <- top5_processed[,-1]

smc.top5.dotplot.processed <- DotPlot(obj.smc.clustered, features = rev(as.character(unique(top5_processed$gene))), group.by = sel.clust,
                                      assay = "RNA") + coord_flip() + theme(axis.text=element_text(size=6)) + theme_classic()



# All genes with p value less than 0.05

markers_genes_sig <- markers_genes %>% filter(p_val_adj < 0.05)


setwd('/Users/judepops/Documents/DISSERTATION/SMCFeatures/')
write.csv(markers_genes_sig, 'Significant_Marker_Genes_SMC.csv')
write.csv(, 'ALL_Genes_SMC.csv')

markers_genes_sig <- read.csv('Significant_Marker_Genes_SMC.csv')




# Rearranging the levels

DimPlot(obj.smc.clustered, reduction = 'umap')



myLevels <- c('1','3','5','12','2','9','6','0','7','4','10','13','8','11')

factor(Idents(obj.smc.clustered), levels= myLevels)
Idents(obj.smc.clustered) <- factor(Idents(obj.smc.clustered), levels= myLevels)
Idents(obj.smc.clustered)


smc.top5.dotplot.processed <- DotPlot(obj.smc.clustered, features = rev(as.character(unique(top5_processed$gene))), 
                                      assay = "RNA") + coord_flip() + theme(axis.text=element_text(size=6)) + theme_classic()


smc.top5.dotplot.processed.rgb <- DotPlot(obj.smc.clustered, 
        features = rev(as.character(unique(top5_processed$gene))), 
        assay = "RNA") + 
  coord_flip() + 
  theme(axis.text=element_text(size=10)) + 
  theme_classic() +  scale_colour_gradient2(low = 'green', mid = "blue", high = "red")



#### Plotting

setwd('/Users/judepops/Documents/DISSERTATION/Figures/')
jpeg(filename="Rearranged_SMC_top5_PROCESSED.jpg",units="cm",width=20,height=30, res=600)

smc.top5.dotplot.processed 

dev.off()

setwd('/Users/judepops/Documents/DISSERTATION/Figures/')
jpeg(filename="RGB_Rearranged_SMC_top5_PROCESSED.jpg",units="cm",width=20,height=30, res=600)

smc.top5.dotplot.processed.rgb 

dev.off()


setwd('/Users/judepops/Documents/DISSERTATION/Figures/')
jpeg(filename="RGB_MARKER_Rearranged_SMC_top5_PROCESSED.jpg",units="cm",width=20,height=30, res=600)

smc.top5.dotplot.processed.rgb 

dev.off()



# Finding gene ontologies
# Load additional packages



DimPlot(obj.smc.clustered, reduction = 'umap')



myLevels <- c('1','3','5','12','2','9','6','0','7','4','10','13','8','11')

factor(Idents(obj.smc.clustered), levels= myLevels)
Idents(obj.smc.clustered) <- factor(Idents(obj.smc.clustered), levels= myLevels)
Idents(obj.smc.clustered)

top5_final<- read.csv('/Users/judepops/Documents/DISSERTATION/SMCFeatures/Literature_Rearranged_Top5_Markers.csv')
top5_final %>% column_to_rownames(X)
top5_final <- top5_final[,-1]

smc.top5.final <- DotPlot(obj.smc.clustered, 
                                          features = rev(as.character(unique(top5_final$gene))), 
                                          assay = "RNA") + 
  coord_flip() + 
  theme(axis.text=element_text(size=10)) + 
  theme_classic() +  scale_colour_gradient2(low = 'green', mid = "blue", high = "red")


setwd('/Users/judepops/Documents/DISSERTATION/Figures/')
jpeg(filename="RGB_MARKER_Rearranged_SMC_top5_PROCESSED.jpg",units="cm",width=20,height=30, res=600)

smc.top5.dotplot.processed.rgb 

dev.off()





# Dotpots for the featureplot genes

genes <- c('TAGLN2', 'ACTA2', 'CNN1', 'MGP', 'VCAN', 'KLF4')

smc.featureplot.dotplot <- DotPlot(obj.smc.clustered, 
                          features = rev(genes), 
                          assay = "RNA") + 
  coord_flip() + 
  theme(axis.text=element_text(size=20)) + 
  theme_classic() +  scale_colour_gradient2(low = 'green', mid = "blue", high = "red")

setwd('/Users/judepops/Documents/DISSERTATION/Figures/')
jpeg(filename="FeaturePlot_Dotplot.jpg",units="cm",width=15,height=25, res=600)

smc.featureplot.dotplot 

dev.off()






