cond <- grepl('/hpc_path_identifier', getwd())
if (cond){
  source('path1/to/config.R')
} else {
  source('path2/to/config.R')
} ; rm(cond)
library(ggplot2)

setwd(rna_cell_type_annotation_objects)
combined.integrated <- readRDS('2.remove_hybrid_cells___integrated_rna_data_annotated_cell_types.rds')
DefaultAssay(combined.integrated) <- "integrated"

combined.integrated <- RunPCA(combined.integrated)
pca <- combined.integrated@reductions$pca
eigValues <- (pca@stdev)^2  ## EigenValues
varExplained <- eigValues / sum(eigValues)
varExplained
ggplot() + geom_point(mapping = aes(x = 1:50, y = varExplained))
cumsum(varExplained)

combined.integrated <- FindNeighbors(combined.integrated, reduction = "pca", dims = 1:30)
combined.integrated <- FindClusters(combined.integrated, resolution = 0.4, graph.name = 'integrated_snn')
combined.integrated <- RunUMAP(combined.integrated, metric = 'correlation', dims = 1:30)

DimPlot(combined.integrated, group.by = 'predicted.id')
DimPlot(combined.integrated, group.by = 'seurat_clusters')

# remove endos and vlmcs
ix <- !(combined.integrated$predicted.id %in% c('Endo', 'VLMC'))
cells_keep <- colnames(combined.integrated)[ix]
combined.integrated <- subset(combined.integrated, cells = cells_keep)

# assign cell types to clusters
cluster_names <- table(combined.integrated$seurat_clusters, combined.integrated$predicted.id)
cluster_names <- colnames(cluster_names)[ apply(cluster_names, 1, which.max) ] 
combined.integrated$seurat_clusters <- as.character(combined.integrated$seurat_clusters)
names(cluster_names) <- as.character( 0:( length(cluster_names)-1 ) )
combined.integrated$cell_type <- cluster_names[combined.integrated$seurat_clusters]

# rename clusters.  
# for each cell type, the bigest cluster gets 1, the next biggest gets 2, etc. 
# example: Exc1 = biggest excitatory cluster, Exc2 = second biggest excitatory cluster
df <- data.table(cell_type = combined.integrated$cell_type, cluster = combined.integrated$seurat_clusters)
df$id <- paste0(df$cell_type, df$cluster)
df <- split(df, f = df$cell_type)
df <- lapply(df, 
             function(x){
               cell_type <- x$cell_type[1]
               tb <- table(x$id)
               tb <- sort(tb, decreasing = T)
               new_names <- data.table(old_name = names(tb), new_name = paste0(cell_type, 1:length(tb)))
               return(new_names)
             })

df <- bind_rows(df)
new_ids <- df$new_name
names(new_ids) <- df$old_name
combined.integrated$cell_type_cluster <- paste0(combined.integrated$cell_type, combined.integrated$seurat_clusters)
combined.integrated$cell_type_cluster <- new_ids[combined.integrated$cell_type_cluster]
table(combined.integrated$cell_type_cluster, combined.integrated$seurat_clusters)


DimPlot(combined.integrated, group.by = 'cell_type_cluster')
DimPlot(combined.integrated, group.by = 'sampID')

setwd(rna_cell_type_annotation_objects)
saveRDS(combined.integrated, '3.pca_lovain_cluster_umap___annotated_integrated_processed_dataset.rds')
