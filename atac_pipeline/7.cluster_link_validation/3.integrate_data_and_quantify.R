# ### have to install archive version of package
# packageurl <- 'https://cran.r-project.org/src/contrib/Archive/Matrix/Matrix_1.3-0.tar.gz'
# install.packages(packageurl, repos=NULL, type="source", lib = '/path/to/desired/install/location')
# library(Matrix, lib.loc = '/path/to/desired/install/location')

cond <- grepl('/hpc_path_identifier', getwd())
if (cond){
  source('path1/to/config.R')
} else {
  source('path2/to/config.R')
} ; rm(cond)
library(data.table)


setwd(cluster_link_validation_objects)
rna <- readRDS('set1_rna_multiome_preprocessed.rds')
atac <- readRDS('set1_atac_multiome_preprocessed.rds')
atac_gex <- readRDS('set2_gene_activities.rds')

atac[["gene_activity"]] <- CreateAssayObject(counts = atac_gex)

DefaultAssay(atac) <- 'gene_activity'
atac <- NormalizeData(object = atac, assay = 'gene_activity', method = 'LogNormalize')

DefaultAssay(atac) <- 'peaks'
atac <- RunTFIDF(atac, assay = 'peaks')


atac <- RunSVD(atac, assay = 'peaks', features = rownames(atac))
atac <- FindNeighbors(atac, reduction = 'lsi', dims = 2:30)
atac <- FindClusters(atac, algorithm = 3)

anchors <- FindTransferAnchors(reference = rna, query = atac, reference.assay = 'rna', 
                               query.assay = 'gene_activity', features = VariableFeatures(rna), 
                               reduction = 'cca')

setwd(cluster_link_validation_objects)
saveRDS(anchors, 'set4_atac_rna_anchors.rds')
anchors <- readRDS('set4_atac_rna_anchors.rds')


atac_cluster_links <- TransferData(anchorset = anchors, refdata = rna$seurat_clusters, weight.reduction = atac[['lsi']], dims = 2:30)


# quantify jac ------------------------------------------------------------
atac_cluster_links$barcode <- rownames(atac_cluster_links)
atac_cluster_links <- as.data.table(atac_cluster_links)
tmp <- data.table('barcode' = rownames(atac@meta.data), 'seurat_clusters' = atac$seurat_clusters)
atac_cluster_links <- full_join(atac_cluster_links, tmp, by = 'barcode')

atac_cluster_links <- as.data.frame(atac_cluster_links)
atac_cluster_links <- split(atac_cluster_links, f = atac_cluster_links$seurat_clusters)
atac_cluster_links <- lapply(atac_cluster_links, 
                             function(x){
                               prediction_scores <- colSums(x[, 2:16])
                               prediction_score_max <- max(prediction_scores)
                               column_names <- names(prediction_scores)
                               predicted_cluster <- column_names[which.max(prediction_scores)]
                               predicted_cluster <- gsub('prediction.score.', '', predicted_cluster)
                               
                               df <- matrix(prediction_scores, nrow = 1)
                               df <- as.data.frame(df)
                               colnames(df) <- column_names
                               df$seurat_clusters <- x$seurat_clusters[1]
                               df$prediction_score_max <- prediction_score_max
                               df$predicted_cluster <- predicted_cluster
                               return(df)
                             }) %>% bind_rows() %>% as.data.table()

tmp <- atac@meta.data
tmp$barcode <- rownames(tmp)
tmp <- as.data.table(tmp)
tmp <- tmp[, .(seurat_clusters, barcode)]
atac_cluster_links <- full_join(atac_cluster_links, tmp, by = 'seurat_clusters')

rna_clusters <- rna@meta.data
rna_clusters$barcode <- rownames(rna_clusters)
rna_clusters <- as.data.table(rna_clusters)

cells.use <- intersect(atac_cluster_links$barcode, rna_clusters$barcode)

atac_cluster_links <- atac_cluster_links[barcode %in% cells.use, ]
rna_clusters <- rna_clusters[barcode %in% cells.use, ]

table(atac_cluster_links$predicted_cluster %in% rna_clusters$seurat_clusters)

jaccardIndex <- function(v1, v2){
  v_int <- intersect(v1, v2) 
  v_unq <- unique(c(v1, v2))
  length(v_int) / length(v_unq)
}

rna_clusters$seurat_clusters <- as.character(rna_clusters$seurat_clusters)
rna_clusters$seurat_clusters <- as.numeric(rna_clusters$seurat_clusters)

atac_cluster_links$predicted_cluster <- as.numeric(atac_cluster_links$predicted_cluster)
table(atac_cluster_links$predicted_cluster %in% rna_clusters$seurat_clusters)

cluster_ix <- unique(atac_cluster_links$predicted_cluster)
cluster_ix <- sort(cluster_ix)
jac_df <- data.table('cluster_id' = cluster_ix, 'jaccard_ix' = 0)

for (i in 1:nrow(jac_df)){
  message(i)
  cluster_ix <- jac_df$cluster_id[i]
  atac_cells <- atac_cluster_links[predicted_cluster == cluster_ix, barcode]
  rna_cells <- rna_clusters[seurat_clusters == cluster_ix, barcode]
  jaccard_ix <- jaccardIndex(atac_cells, rna_cells)
  jac_df$jaccard_ix[i] <- jaccard_ix
}



# grubman filter ----------------------------------------------------------
tmp <- as.data.frame(atac_cluster_links)
tmp <- split(tmp, f = tmp$predicted_cluster)
tmp <- lapply(tmp,
              function(x){
                sums <- colSums( x[, grep('^prediction.score.[[:digit:]]+$', colnames(x))])
                sums[length(sums)+1] <- x$predicted_cluster[1]
                names(sums)[length(sums)] <- 'predicted_cluster'
                return(sums)
              }) %>% do.call(what = rbind) %>% as.data.frame()
            
grubmanFilter <- function(scores){
  c <- sort(scores, decreasing = T)[2:1]
  diff(c) / max(c)
}

grubman_score <- data.frame('cluster_id' = tmp$predicted_cluster, 'grubman_score' = 0)
for (i in 1:nrow(grubman_score)){
  cluster_ix <- grubman_score$cluster_id[i]
  scores <- tmp[tmp$predicted_cluster == cluster_ix, 
                grep('^prediction.score.[[:digit:]]+$', colnames(tmp))]
  scores <- as.numeric(scores)
  grubman_score$grubman_score[i] <- grubmanFilter(scores)
}

jac_df <- full_join(jac_df, grubman_score, by = 'cluster_id')

setwd(cluster_link_validation_objects)
saveRDS(jac_df, '3.integrate_data___jaccard_ix.rds')



