library(ggplot2)
library(reshape2)

Metric_Results <- read.csv('/Users/judepops/Documents/DISSERTATION/Benchmarking/Final_Results/Metric_Results.csv')


colnames(Metric_Results) <- NULL




my.names <- Metric_Results[1,]
colnames(Metric_Results) <- my.names
Metric_Results = Metric_Results[-1,]

Metric_Results = Metric_Results[-9,]
Metric_Results = Metric_Results[-9,]

## convert to long format with tidyr::pivot_longer
Metric_Results <- pivot_longer(Metric_Results, cols = colnames(Metric_Results))

Metric_Results %>% column_to_rownames(Metric)



ggplot(data=test_data_long_tidyr,
       aes(x=date, y=value, colour=name)) +
  geom_line() ## output not shown, it's equivalent to the below graph (with a tiny difference in the legend title)


library("tidyverse")
df <- Metric_Results %>%
  select(Metric, Harmony, Scanorama, rPCA) %>%
  gather(key = "variable", value = "value", -Metric)
head(df)


ggplot(df, aes(x = Metric, y = value, group = 1)) + 
  geom_line(aes(color = variable, linetype = variable)) + 
  scale_color_manual(values = c("red", "blue", 'yellow'))







Metric <- c('NMI[Louvain]', 'NMI[Leiden]', 'ARI[Louvain]', 'ARI[Leiden]', 'ASW', 'kBET', 'LISI/2')
Harmony <- as.numeric(c('0.6774049','0.6863792', '0.4646822', '0.4652630', '0.08549741', '0.996', '0.668813'))
rPCA <- as.numeric(c('0.6571805','0.6850007', '0.4519684', '0.4667722', '0.26104108', '0.9994', '0.693041'))
Scanorama <- as.numeric(c('0.6148498','0.6495001', '0.3892980', '0.4051145', '0.26058368', '0.9999', '0.7200455'))

Metric_Results <- data.frame(Metric, Harmony,Scanorama,rPCA)

# 
# df <- Metric_Results %>%
#   select(Metric, Harmony, Scanorama, rPCA) %>%
#   gather(key = "variable", value = "value", -Metric)
# head(df)
# 
# 
# ggplot(df, aes(x = Metric, y = value)) + 
#   geom_line(aes(color = variable, linetype = variable)) + 
#   scale_color_manual(values = c("red", "blue", 'yellow'))


library(ggplot2)
#> Warning: package 'ggplot2' was built under R version 4.0.5
ggplot(Metric_Results, aes(Metric)) +      
  geom_line(aes(y = Harmony), colour = "red", group = "Harmony") +
  geom_line(aes(y = Scanorama), colour = "blue", group = "Scanorama") +
  geom_line(aes(y = rPCA), colour = "grey", group = "rPCA") + theme_minimal() 

library(tidyr)
ggpLong <- Metric_Results %>% pivot_longer(cols = Harmony:rPCA, names_to = "Method", values_to = "Value")
pmetric <- ggplot(ggpLong, aes(x = Metric, y = Value, color = Method, group = Method)) + geom_line() + theme_classic() + 

my_colors <- c("Harmony" = "red", "Scanorama" = "blue", "rPCA" = "#00A300")
ggpLong <- Metric_Results %>% pivot_longer(cols = Harmony:rPCA, names_to = "Method", values_to = "Value")
pmetric <- ggplot(ggpLong, aes(x = Metric, y = Value, color = Method, group = Method)) + 
  geom_line(size = 1) + 
  theme_classic() + 
  scale_color_manual(values = my_colors) + theme(text = element_text(size=20)) 





setwd('/Users/judepops/Documents/DISSERTATION/Benchmarking/Final_Results/')

jpeg(filename="Metrics_Method_2.jpg",units="cm",width=40,height=20, res=600)

pmetric

dev.off()




# Batch_Results

Batch <- c('Biological', 'Clustering', 'Technical')
Harmony <- as.numeric(c('0.573432325',
                        '0.08549741',
                        '0.8324065'))
Scanorama <- as.numeric(c('0.5146906',
                     '0.26058368',
                     '0.85997275'))
rPCA <- as.numeric(c('0.56523045',
                          '0.26104108',
                          '0.8462205'))


Batch_Results <- data.frame(Batch, Harmony,Scanorama,rPCA)




# DO A BOX PLOT FOR THIS ONE
                                                                                                                                                                                

library(dplyr)
library(tidyr)
library(ggplot2)

library(tidyr)

df_long <- gather(Batch_Results, key = "Method", value = "Value", -Batch)

# Setting colors
colours <- c("Harmony" = "red", "Scanorama" = "blue", "rPCA" = "#00A300")


# plot the bar graph
Batch_Results_Bar <- ggplot(df_long, aes(x = Batch, y = Value, fill = Method)) +
  geom_bar(stat = "identity", position = "dodge") + scale_fill_manual(values = colours) +
  labs(x = "Batch", y = "Score", fill = "Method") + theme_classic() + theme(text = element_text(size=23))   





# Printing the image

setwd('/Users/judepops/Documents/DISSERTATION/Benchmarking/Final_Results/')

jpeg(filename="BioTechClust_Results_2.jpg",units="cm",width=25,height=15, res=600)

Batch_Results_Bar

dev.off()








































# Creating overall results plot
Integration <- c('Harmony','Scanorama','rPCA')
Overall <- as.numeric(c('0.497112078', '0.545082343', '0.557497343'))
Final_Results <- data.frame(Integration, Overall)

# Setting the levels so that rCPA is on top

Final_Results$Integration <- factor(Final_Results$Integration, levels = c("Harmony", "Scanorama", "rPCA"))

# Changing colours and plotting

ggplot(Final_Results, aes(x = Overall, y = Integration, fill = Integration)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("red", "blue", "#00A300")) +
  labs(title = "Score Comparison of Integration Methods",
       x = "Overall Score",
       y = "Integration Method") +
  theme(plot.title = element_text(hjust = 0.5)) +theme(text = element_text(size=23)) 





# Saving the plot

setwd('/Users/judepops/Documents/DISSERTATION/Benchmarking/Final_Results/')
jpeg(filename="Overall_Benchmark_Results_2.jpg",units="cm",width=25,height=15, res=600)

ggplot(Final_Results, aes(x = Overall, y = Integration, fill = Integration)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("red", "blue", "#00A300")) +
  labs(
       x = "Overall Score",
       y = "Integration Method") +
  theme(plot.title = element_text(hjust = 0.5)) + theme_classic() + theme(text = element_text(size=23)) 

dev.off()







# Creating a plot with biological axis and technical axis




Technical <- as.numeric(c('0.8324065','0.85997275', '0.8462205'))
Bio.Clust <- as.numeric(c('0.475845342','0.429988293', '0.463833993'))
Technique <- c('Harmony', 'Scanorana', 'rPCA')





Benchmark_Scatter <- data.frame(Technique, Technical,Bio.Clust)



ggplot(Benchmark_Scatter, aes(x = Technical, y = Bio.Clust, color = Technique)) +
  geom_point() +
  labs(x = "Technical", y = "Bio.Clust") +
  ggtitle("Scatter Plot") +
  theme_bw()



df <- data.frame(
  Technique = c("Harmony", "Scanorama", "rPCA"),
  Technical = c(0.8324065, 0.8599727, 0.8462205),
  Bio.Clust = c(0.4758453, 0.4299883, 0.4638340)
)

# Define colors for the plot
colours <- c("Harmony" = "red", "Scanorama" = "blue", "rPCA" = "#00A300")

# Create the scatter plot
Benchmark <- ggplot(df, aes(x = Technical, y = Bio.Clust, color = Technique)) +
  geom_point() +
  scale_x_continuous(limits = c(0, 1)) +
  scale_y_continuous(limits = c(0, 1)) +
  scale_color_manual(values = colours) +
  labs(x = "Technical", y = "Bio.Clust") +
  ggtitle("Scatter Plot") +
  theme_classic()


setwd('/Users/judepops/Documents/DISSERTATION/Benchmarking/Final_Results/')
jpeg(filename="Scatter_Benchmark_Results.jpg",units="cm",width=20,height=10, res=600)

Benchmark

dev.off()




