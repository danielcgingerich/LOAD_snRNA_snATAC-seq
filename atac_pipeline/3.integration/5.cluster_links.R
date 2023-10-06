cond <- grepl('/hpc_path_identifier', getwd())
if (cond){
  source('path1/to/config.R')
} else {
  source('path2/to/config.R')
} ; rm(cond)

library(stringr)

setwd(atac_integration_objects)
atac <- readRDS('4.run_harmony___integrated_atac_data.rds')



# find atac clusters ------------------------------------------------------
atac <- RunUMAP(atac, reduction = 'harmony.all.features', dims = 2:30)
atac <- FindNeighbors(atac, reduction = 'harmony.all.features', dims = 2:30)
atac <- FindClusters(atac, algorithm = 3, resolution = 1) # might want to tweak res to see if there is common structure across datasets



# link rna clusters -------------------------------------------------------
setwd(atac_cell_type_annotation_objects)
pred_clust <- list.files(pattern = '3.annotate_atac.*_cluster_predictions.rds')
pred_clust <- lapply(pred_clust, readRDS) %>% bind_rows
pred_clust <- pred_clust[colnames(atac),] # bad predictions filtered before harmony integration
pred_clust$atac_cluster <- atac$seurat_clusters

pred_clust <- split(pred_clust, f = pred_clust$atac_cluster)
tallies <- lapply(pred_clust, 
                  function(x){
                    scores <- grepl('prediction.score.*[[:digit:]]$', colnames(x))
                    sums <- colSums(x[,scores]) / nrow(x)
                    
                    predicted_id <- names(sums)[which.max(sums)] 
                    predicted_id <- gsub(pattern = 'prediction.score.', replacement = '', x = predicted_id)
                    
                    prediction_score_max <- max(sums) 
                    
                    grubman_filter <- sort(unname(sums), decreasing = T)[2]
                    grubman_filter <- (prediction_score_max - grubman_filter) / prediction_score_max
                    
                    data.frame('atac_cluster' = x$atac_cluster[1], 'cluster_links' = predicted_id,
                               'cluster_score_max' = prediction_score_max, 'cluster_score_filter' = grubman_filter)
                    }) %>% bind_rows %>% as.data.frame

atac$bcd <- rownames(atac@meta.data)
atac@meta.data <- full_join(atac@meta.data, tallies, 
                            by = c('seurat_clusters' = 'atac_cluster'))
rownames(atac@meta.data) <- atac$bcd

# save work ---------------------------------------------------------------
setwd(atac_integration_objects)
saveRDS(atac, '5.cluster_links___atac_integrated_and_linked.rds')

