# Libraries
library(ggplot2)
library(ggrepel)

# Data

proteomics.calcified <- read.csv('/Users/judepops/Documents/DISSERTATION/Proteomics/Proteomics_Calcification.csv')

# Volcano Plots

# Set color palette
colors <- c("blue", "red")

# Define threshold for significant genes
sig_threshold <- -log10(0.05)

# Create volcano plot
ggplot(data = proteomics.calcified, aes(x = Log2.Fold.Change, y = -log10(P.Value))) +
  geom_point(aes(color = ifelse(Log2.Fold.Change > 0, colors[2], colors[1])),
             size = 2) +
  scale_color_identity(guide = "legend", labels = c("Down-regulated", "Up-regulated"),
                       breaks = c(colors[1], colors[2])) +
  geom_hline(yintercept = sig_threshold, linetype = "dashed", color = "gray") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray") +
  xlab("Log2 Fold Change") + ylab("-log10 P-value") +
  ggtitle("Volcano Plot") +
  theme_classic()

# Add gene labels to the plot
proteomics.calcified.dotplot <- ggplot(data = proteomics.calcified, aes(x = Log2.Fold.Change, y = -log10(P.Value))) +
  geom_point(aes(color = ifelse(Log2.Fold.Change > 0, colors[2], colors[1])),
             size = 2) +
  scale_color_identity(guide = "legend", labels = c("Down-regulated", "Up-regulated"),
                       breaks = c(colors[1], colors[2])) +
  geom_hline(yintercept = sig_threshold, linetype = "dashed", color = "gray") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray") +
  xlab("Log2 Fold Change") + ylab("-log10 P-value") +
  ggtitle("Volcano Plot") +
  theme_classic() +
  geom_text_repel(aes(label = Gene), size = 3, point.padding = 0.5,
                  segment.alpha = 0, box.padding = 0.5) + NoLegend() 


# Saving the image

setwd('/Users/judepops/Documents/DISSERTATION/Proteomics')

jpeg(filename="SMC_CALCIFICAITON_VOLCANO.jpg",units="cm",width=20,height=10, res=600)

proteomics.calcified.dotplot

dev.off()


# Symptoms Volcano Plot

proteomics.symptoms <- read.csv('/Users/judepops/Documents/DISSERTATION/Proteomics/Proteomics_Symptoms.csv')

# Create volcano plot
ggplot(data = proteomics.symptoms, aes(x = Log2.Fold.Change, y = -log10(P.Value))) +
  geom_point(aes(color = ifelse(Log2.Fold.Change > 0, colors[2], colors[1])),
             size = 2) +
  scale_color_identity(guide = "legend", labels = c("Down-regulated", "Up-regulated"),
                       breaks = c(colors[1], colors[2])) +
  geom_hline(yintercept = sig_threshold, linetype = "dashed", color = "gray") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray") +
  xlab("Log2 Fold Change") + ylab("-log10 P-value") +
  ggtitle("Volcano Plot") +
  theme_classic()

# Add gene labels to the plot
proteomics.symptoms.dotplot <- ggplot(data = proteomics.symptoms, aes(x = Log2.Fold.Change, y = -log10(P.Value))) +
  geom_point(aes(color = ifelse(Log2.Fold.Change > 0, colors[2], colors[1])),
             size = 2) +
  scale_color_identity(guide = "legend", labels = c("Down-regulated", "Up-regulated"),
                       breaks = c(colors[1], colors[2])) +
  geom_hline(yintercept = sig_threshold, linetype = "dashed", color = "gray") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray") +
  xlab("Log2 Fold Change") + ylab("-log10 P-value") +
  ggtitle("Volcano Plot") +
  theme_classic() +
  geom_text_repel(aes(label = Gene), size = 3, point.padding = 0.5,
                  segment.alpha = 0, box.padding = 0.5) + NoLegend() 



setwd('/Users/judepops/Documents/DISSERTATION/Proteomics')

jpeg(filename="SMC_SYMPTOMS_VOLCANO.jpg",units="cm",width=20,height=10, res=600)

proteomics.symptoms.dotplot

dev.off()



# Sex Diferences
proteomics.sex <- read.csv('/Users/judepops/Documents/DISSERTATION/Proteomics/Proteomics_Sex.csv')

ggplot(data = proteomics.sex, aes(x = Log2.Fold.Change, y = -log10(P.Value))) +
  geom_point(aes(color = ifelse(Log2.Fold.Change > 0, colors[2], colors[1])),
             size = 2) +
  scale_color_identity(guide = "legend", labels = c("Down-regulated", "Up-regulated"),
                       breaks = c(colors[1], colors[2])) +
  geom_hline(yintercept = sig_threshold, linetype = "dashed", color = "gray") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray") +
  xlab("Log2 Fold Change") + ylab("-log10 P-value") +
  ggtitle("Volcano Plot") +
  theme_classic()

# Add gene labels to the plot
proteomics.sex.dotplot <- ggplot(data = proteomics.sex, aes(x = Log2.Fold.Change, y = -log10(P.Value))) +
  geom_point(aes(color = ifelse(Log2.Fold.Change > 0, colors[2], colors[1])),
             size = 2) +
  scale_color_identity(guide = "legend", labels = c("Down-regulated", "Up-regulated"),
                       breaks = c(colors[1], colors[2])) +
  geom_hline(yintercept = sig_threshold, linetype = "dashed", color = "gray") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray") +
  xlab("Log2 Fold Change") + ylab("-log10 P-value") +
  ggtitle("Volcano Plot") +
  theme_classic() +
  geom_text_repel(aes(label = Gene), size = 3, point.padding = 0.5,
                  segment.alpha = 0, box.padding = 0.5) + NoLegend() 


setwd('/Users/judepops/Documents/DISSERTATION/Proteomics')

jpeg(filename="SMC_SEX_VOLCANO.jpg",units="cm",width=20,height=10, res=600)

proteomics.sex.dotplot

dev.off()


# Core periphery volcano plot

proteomics.coreperiph <- read.csv('/Users/judepops/Documents/DISSERTATION/Proteomics/Proteomics_CorePeriph.csv')

ggplot(data = proteomics.coreperiph, aes(x = Log2.Fold.Change, y = -log10(P.Value))) +
  geom_point(aes(color = ifelse(Log2.Fold.Change > 0, colors[2], colors[1])),
             size = 2) +
  scale_color_identity(guide = "legend", labels = c("Down-regulated", "Up-regulated"),
                       breaks = c(colors[1], colors[2])) +
  geom_hline(yintercept = sig_threshold, linetype = "dashed", color = "gray") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray") +
  xlab("Log2 Fold Change") + ylab("-log10 P-value") +
  ggtitle("Volcano Plot") +
  theme_classic()

# Add gene labels to the plot
proteomics.coreperiph.dotplot <- ggplot(data = proteomics.coreperiph, aes(x = Log2.Fold.Change, y = -log10(P.Value))) +
  geom_point(aes(color = ifelse(Log2.Fold.Change > 0, colors[2], colors[1])),
             size = 2) +
  scale_color_identity(guide = "legend", labels = c("Down-regulated", "Up-regulated"),
                       breaks = c(colors[1], colors[2])) +
  geom_hline(yintercept = sig_threshold, linetype = "dashed", color = "gray") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray") +
  xlab("Log2 Fold Change") + ylab("-log10 P-value") +
  ggtitle("Volcano Plot") +
  theme_classic() +
  geom_text_repel(aes(label = Gene), size = 3, point.padding = 0.5,
                  segment.alpha = 0, box.padding = 0.5) + NoLegend() 


setwd('/Users/judepops/Documents/DISSERTATION/Proteomics')

jpeg(filename="SMC_COREPERIPH_VOLCANO.jpg",units="cm",width=25,height=15, res=600)

proteomics.coreperiph.dotplot

dev.off()


legend <- get_legend(p)
(proteomics.coreperiph.dotplot)
library(cowplot)


  
  
  
  
  
