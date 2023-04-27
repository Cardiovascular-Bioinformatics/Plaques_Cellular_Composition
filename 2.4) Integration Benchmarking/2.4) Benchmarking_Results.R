# Results

# Working directory
setwd('/Users/judepops/Documents/DISSERTATION/Benchmarking/Final_Results')


# NMI/ARI results

## Louvain 0.8
saveRDS(ARI_RESULTS_LOUVAIN, 'ARI_RESULTS_LOUVAIN.rds')
saveRDS(NMI_RESULTS_LOUVAIN, 'NMI_RESULTS_LOUVAIN.rds')
## Leiden
saveRDS(NMI_RESULTS_LEIDEN, 'NMI_RESULTS_LEIDEN.rds')
saveRDS(ARI_RESULTS_LEIDEN, 'ARI_RESULTS_LEIDEN.rds')


# LISI results
LISI_RESULTS <- readRDS('LISI_RESULTS.rds')

# ASW score results
ASW_RESULTS <- readRDS('ASW_RESULTS.rds')

# kBET results
KBET_RESULTS <- readRDS('KBET_RESULTS.rds')


# Dataframe of all the results

Metric <- c('NMI[Louvain]','NMI[Leiden]', 'ARI[Louvain]', 'ARI[Leiden]', 'LISI', 'ASW', 'kBET')

Harmony.Benchmark <- c('0.6774049','0.6863792', '0.4646822', '0.4652630', '1.337626', '0.08549741', '0.996')
rPCA.Benchmark <- c('0.6571805','0.6850007', '0.4519684', '0.4667722', '1.386082', '0.26104108', '0.9994')
Scanorama.Benchmark <- c('0.6148498','0.6495001', '0.3892980', '0.4051145', '1.440091', '0.26058368', '0.9999')
Metric.Ranges <- c('0:1','0:1', '-1:1', '-1:1', '1:2', '-1:1', '0:1')


Metric_Results <- data.frame(Metric,Metric.Ranges, Harmony.Benchmark,Scanorama.Benchmark,rPCA.Benchmark)
Metric_Results$Harmony.Benchmark <- as.numeric(Metric_Results$Harmony.Benchmark)
Metric_Results$Scanorama.Benchmark <- as.numeric(Metric_Results$Scanorama.Benchmark)
Metric_Results$rPCA.Benchmark <- as.numeric(Metric_Results$rPCA.Benchmark)

mean(Metric_Results$Harmony.Benchmark)
mean(Metric_Results$Scanorama.Benchmark)
mean(Metric_Results$rPCA.Benchmark)

write.csv(Metric_Results, 'Metric_Results.csv')









