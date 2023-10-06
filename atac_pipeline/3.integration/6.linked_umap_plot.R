cond <- grepl('/hpc_path_identifier', getwd())
if (cond){
  source('path1/to/config.R')
} else {
  source('path2/to/config.R')
} ; rm(cond)




# description -------------------------------------------------------------
# instructions on how to make umap plots of rna/atac datasets that have linked clusters.
# legend colors and labelling for each cluster will match for both umaps
setwd(atac_integration_objects)
atac <- readRDS('5.cluster_links___atac_integrated_and_linked.rds') 

setwd(rna_cell_type_annotation_objects)
rna <- readRDS('3.pca_lovain_cluster_umap___annotated_integrated_processed_dataset.rds')

rna$cell_type_cluster <- as.character(rna$cell_type_cluster)
clusters <- unique(rna$cell_type_cluster)
levels <- sort(clusters)
Idents(rna) <- 'cell_type_cluster'
levels(rna) <- levels
rna$cell_type_cluster <- factor(rna$cell_type_cluster, levels = levels)

atac.umap <- DimPlot(atac, reduction = 'umap', group.by = 'cluster_links', pt.size = 0.1)
rna.umap <- DimPlot(rna, reduction = 'umap', group.by = 'cell_type_cluster', pt.size = 0.1)

atac.levels <- levels(as.factor(atac$cluster_links))
rna.levels <- levels(rna$cell_type_cluster)

rna.cols <- as.data.frame(ggplot_build(rna.umap)$data)
atac.cols <- as.data.frame(ggplot_build(atac.umap)$data)

rna.cols <- data.frame("color"=rna.cols$colour, "group"=rna.cols$group)
rna.cols <- unique(rna.cols)
rna.cols <- rna.cols[order(rna.cols$group),]

atac.cols <- data.frame("color"=atac.cols$colour, "group"=atac.cols$group)
atac.cols <- unique(atac.cols)
atac.cols <- atac.cols[order(atac.cols$group),]

pos <- match(atac.levels, rna.levels)
atac.cols$color <- rna.cols$color[pos]

rm(rna.umap, atac.umap)

rna.umap <- DimPlot(rna, reduction = 'umap', group.by = 'cell_type_cluster', pt.size = 0.1, label = T, cols = rna.cols$color, repel = T)
atac.umap <- DimPlot(atac, reduction = 'umap', group.by = 'cluster_links', pt.size = 0.1, label = T, cols = atac.cols$color, repel = T)

rna.umap + atac.umap
