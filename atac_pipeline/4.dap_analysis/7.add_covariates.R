cond <- grepl('/hpc_path_identifier', getwd())
if (cond){
  source('path1/to/config.R')
} else {
  source('path2/to/config.R')
} ; rm(cond)

setwd(atac_dap_analysis_objects)


setwd(atac_cell_type_annotation_objects)
rds_list <- list.files(pattern = '_cell_type_predictions.rds')
sample_id <- gsub('.*___|_cell_.*', '', rds_list)
rds_list2 <- paste0('3.annotate_atac_cell_types___', sample_id, '_cluster_predictions.rds')

#load objects\~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
setwd(atac_integration_objects)
harmony <- readRDS('5.cluster_links___atac_integrated_and_linked.rds')
setwd(atac_preprocessing_objects)
metadata <- readRDS('3.merge_metadata___merged_cell_metadata.rds')
pred_celltype <- list()
pred_cluster <- list()
for (i in 1:length(rds_list)){
  setwd(atac_cell_type_annotation_objects)
  pred_celltype[[i]] <- readRDS(rds_list[i])
  pred_celltype[[i]]$barcode <- rownames(pred_celltype[[i]])
  pred_cluster[[i]] <- readRDS(rds_list2[i])
  pred_cluster[[i]]$barcode <- rownames(pred_cluster[[i]])
}



#merge metadata\~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
pred_celltype <- bind_rows(pred_celltype)
pred_cluster <- bind_rows(pred_cluster)
 
colnames(pred_cluster) <- paste0(colnames(pred_cluster), '_cluster')
colnames(pred_celltype) <- paste0(colnames(pred_celltype), '_celltype')

metadata <- metadata[rownames(pred_cluster),]

metadata <- bind_cols(metadata, pred_celltype, pred_cluster)

#add cluster_links to metadata\~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
metadata <- metadata[colnames(harmony),]
metadata$cluster_links <- harmony$cluster_links
metadata$seurat_clusters <- harmony$seurat_clusters

#add nuclei proportion to metadata\~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
levels <- levels(as.factor(metadata$cluster_links))
metadata <- split(metadata, f = metadata$sample_id)
metadata <- lapply(metadata, 
                   function(x){
                     x$cluster_links <- factor(x$cluster_links, levels = levels)
                     p <- table(x$cluster_links) / nrow(x)
                     for (i in 1:length(levels)){
                       x[,paste0(levels[i], '.proportion')] <- p[i]
                     }
                     return(x)
                   })
metadata <- bind_rows(metadata)

#save\~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
saveRDS(metadata, '7.add_covarates___metadata_for_dap_analysis.rds')
