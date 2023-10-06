
##################################################################################################
# code provided by E Keats Shwab
#filtered cell type wilcoxon

source('path/to/config.R')
library(ggplot2)
library(dplyr)
library(ggrepel)
library(cowplot)
library(gtools)

rna <- readRDS('path/to/combined_SCT_for_DE_seq_sat.rds')
setwd('path/to/LOAD')


sample_percent=0.20
n_iterations = 30
proportion_df <- data.frame()
for(i in 1:n_iterations){
message(i)
cur_sample <- rna@meta.data[sample(rownames(rna@meta.data), round(sample_percent*ncol(rna))),]
meta_list <- cur_sample %>%
dplyr::group_split(sampID)
temp <- lapply(meta_list, function(meta){
df <- as.data.frame(meta$cell.type %>% table / nrow(meta))
colnames(df) <- c('cluster', 'proportion')
   df$sampID <- unique(meta$sampID)
    df$Diagnosis <- unique(meta$orig.ident)
    df
 })
 cur_df <- Reduce(rbind, temp)
 cur_df$iteration <- i
 proportion_df <- rbind(proportion_df, cur_df)
 }


clusters <- unique(rna@meta.data$cell.type)
clusters <- clusters[order(clusters)]

pvals <- sapply(clusters, function(cur_cluster){
result <- wilcox.test(
proportion_df %>% subset(Diagnosis=='LOAD' & cluster==cur_cluster) %>% .$proportion,
proportion_df %>% subset(Diagnosis=='Normal' & cluster==cur_cluster) %>% .$proportion
)

result$p.value
})

#result <- wilcox.test(
#proportion_df %>% subset(Diagnosis=='LOAD' & cluster== 'Exc1') %>% .$proportion,
#proportion_df %>% subset(Diagnosis=='Normal' & cluster== 'Exc1') %>% .$proportion
#)

result_table <- data.frame(
pval = as.numeric(pvals),
cluster = clusters,
signif = as.character(stars.pval(pvals))
)



setwd('/path/to/LOAD')

write.csv(result_table, 'filtered.cell.type.proportions.wilcoxon.30.iter.methods.csv')

#########################
#filtered cluster wilcoxon

sample_percent=0.20
n_iterations = 30
proportion_df <- data.frame()
for(i in 1:n_iterations){
message(i)
cur_sample <- rna@meta.data[sample(rownames(rna@meta.data), round(sample_percent*ncol(rna))),]
meta_list <- cur_sample %>%
dplyr::group_split(sampID)
temp <- lapply(meta_list, function(meta){
df <- as.data.frame(meta$cell.type.number %>% table / nrow(meta))
colnames(df) <- c('cluster', 'proportion')
   df$sampID <- unique(meta$sampID)
    df$Diagnosis <- unique(meta$orig.ident)
    df
 })
 cur_df <- Reduce(rbind, temp)
 cur_df$iteration <- i
 proportion_df <- rbind(proportion_df, cur_df)
 }


clusters <- unique(rna@meta.data$cell.type.number)
clusters <- clusters[order(clusters)]

pvals <- sapply(clusters, function(cur_cluster){
result <- wilcox.test(
proportion_df %>% subset(Diagnosis=='LOAD' & cluster==cur_cluster) %>% .$proportion,
proportion_df %>% subset(Diagnosis=='Normal' & cluster==cur_cluster) %>% .$proportion
)

result$p.value
})


result_table <- data.frame(
pval = as.numeric(pvals),
cluster = clusters,
signif = as.character(stars.pval(pvals))
)



setwd('/path/to/LOAD')

write.csv(result_table, 'filtered.cluster.proportions.wilcoxon.30.iter.methods.csv')


####################################
#cell type unfiltered wilcoxon

source('/path/to/config.R')
library(ggplot2)
library(dplyr)
library(ggrepel)
library(cowplot)
library(gtools)

rna <- readRDS('/path/to/combined_integrated_labeled_filt_noEndoVLMC.rds')
setwd('/path/to/LOAD')
rna



sample_percent=0.20
n_iterations = 30
proportion_df <- data.frame()
for(i in 1:n_iterations){
message(i)
cur_sample <- rna@meta.data[sample(rownames(rna@meta.data), round(sample_percent*ncol(rna))),]
meta_list <- cur_sample %>%
dplyr::group_split(sampID)
temp <- lapply(meta_list, function(meta){
df <- as.data.frame(meta$predicted.id %>% table / nrow(meta))
colnames(df) <- c('cluster', 'proportion')
   df$sampID <- unique(meta$sampID)
    df$Diagnosis <- unique(meta$orig.ident)
    df
 })
 cur_df <- Reduce(rbind, temp)
 cur_df$iteration <- i
 proportion_df <- rbind(proportion_df, cur_df)
 }


clusters <- unique(rna@meta.data$predicted.id)
clusters <- clusters[order(clusters)]

pvals <- sapply(clusters, function(cur_cluster){
result <- wilcox.test(
proportion_df %>% subset(Diagnosis=='LOAD' & cluster==cur_cluster) %>% .$proportion,
proportion_df %>% subset(Diagnosis=='Normal' & cluster==cur_cluster) %>% .$proportion
)

result$p.value
})

#result <- wilcox.test(
#proportion_df %>% subset(Diagnosis=='LOAD' & cluster== 'Exc1') %>% .$proportion,
#proportion_df %>% subset(Diagnosis=='Normal' & cluster== 'Exc1') %>% .$proportion
#)

result_table <- data.frame(
pval = as.numeric(pvals),
cluster = clusters,
signif = as.character(stars.pval(pvals))
)


setwd('/path/to/LOAD')

write.csv(result_table, 'unfiltered.cell.type.proportions.wilcoxon.30.iter.methods.csv')


##################################################################################################
#integrated atac cell type wilcoxon

source('/path/to/config.R')
library(ggplot2)
library(dplyr)
library(ggrepel)
library(cowplot)
library(gtools)


atac <- readRDS('/path/to/set4_srt_consensus.peaks.cluster-specific.-c=2.merged.rds')
meta <- readRDS('/path/to/set1_merged.cell.metadata.2.rds')

atac <- atac[,rownames(meta)]

atac <- AddMetaData(atac, metadata = meta)




atac@meta.data$diagnosis <- as.character(atac@meta.data$diagnosis)




sample_percent=0.20
n_iterations = 30
proportion_df <- data.frame()
for(i in 1:n_iterations){
message(i)
cur_sample <- atac@meta.data[sample(rownames(atac@meta.data), round(sample_percent*ncol(atac))),]
meta_list <- cur_sample %>%
dplyr::group_split(sample.id)
temp <- lapply(meta_list, function(meta){
df <- as.data.frame(meta$predicted.id %>% table / nrow(meta))
colnames(df) <- c('cluster', 'proportion')
   df$sampID <- unique(meta$sample.id)
    df$Diagnosis <- unique(meta$diagnosis)
    df
 })
 cur_df <- Reduce(rbind, temp)
 cur_df$iteration <- i
 proportion_df <- rbind(proportion_df, cur_df)
 }


clusters <- unique(atac@meta.data$predicted.id)
clusters <- clusters[order(clusters)]

pvals <- sapply(clusters, function(cur_cluster){
result <- wilcox.test(
proportion_df %>% subset(Diagnosis=='LOAD' & cluster==cur_cluster) %>% .$proportion,
proportion_df %>% subset(Diagnosis=='Normal' & cluster==cur_cluster) %>% .$proportion
)

result$p.value
})


result_table <- data.frame(
pval = as.numeric(pvals),
cluster = clusters,
signif = as.character(stars.pval(pvals))
)



setwd('/path/to/LOAD')

write.csv(result_table, 'atac.integrated.cell.type.proportions.wilcoxon.30.iter.methods.csv')


###############################################
### ATAC cluster proportions wilcoxon


sample_percent=0.20
n_iterations = 30
proportion_df <- data.frame()
for(i in 1:n_iterations){
message(i)
cur_sample <- atac@meta.data[sample(rownames(atac@meta.data), round(sample_percent*ncol(atac))),]
meta_list <- cur_sample %>%
dplyr::group_split(sample.id)
temp <- lapply(meta_list, function(meta){
df <- as.data.frame(meta$cluster.links %>% table / nrow(meta))
colnames(df) <- c('cluster', 'proportion')
   df$sampID <- unique(meta$sample.id)
    df$Diagnosis <- unique(meta$diagnosis)
    df
 })
 cur_df <- Reduce(rbind, temp)
 cur_df$iteration <- i
 proportion_df <- rbind(proportion_df, cur_df)
 }


clusters <- unique(atac@meta.data$cluster.links)
clusters <- clusters[order(clusters)]

pvals <- sapply(clusters, function(cur_cluster){
result <- wilcox.test(
proportion_df %>% subset(Diagnosis=='LOAD' & cluster==cur_cluster) %>% .$proportion,
proportion_df %>% subset(Diagnosis=='Normal' & cluster==cur_cluster) %>% .$proportion
)

result$p.value
})


result_table <- data.frame(
pval = as.numeric(pvals),
cluster = clusters,
signif = as.character(stars.pval(pvals))
)



setwd('/path/to/LOAD')

write.csv(result_table, 'atac.integrated.cluster.proportions.wilcoxon.30.methods.csv')

###