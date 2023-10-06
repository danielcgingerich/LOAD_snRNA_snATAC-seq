cond <- grepl('/hpc_path_identifier', getwd())
if (cond){
  source('path1/to/config.R')
} else {
  source('path2/to/config.R')
} ; rm(cond)


setwd(atac_preprocessing_objects)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
atac <- list.files(pattern = '_seurat_object_qc_filtered.rds')
atac <- lapply(atac, readRDS)

metadata <- lapply(atac, FUN = function(x){x@meta.data})
metadata <- bind_rows(metadata)

saveRDS(metadata, '3.merge_metadata___merged_cell_metadata.rds')
