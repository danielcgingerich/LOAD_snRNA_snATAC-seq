cond <- grepl('/hpc_path_identifier', getwd())
if (cond){
  source('path1/to/config.R')
} else {
  source('path2/to/config.R')
} ; rm(cond)

# set paths ---------------------------------------------------------------

setwd(atac_dap_analysis_objects)
rds_list <- list.files(pattern = 'consensus_peak_object')
sample_id <- gsub('.*___|_consensus.*', '', rds_list)

setwd(atac_preprocessing_objects)
rds_list2 <- paste0('4.dim_reduction_clustering___', sample_id, '_seurat_qc_dim_reduced.rds')

# prepare objects for merge -----------------------------------------------

message('\nloading objects \n')
atac.list <- list()
for (i in 1:length(rds_list)){
  message(i, ": ", sample_id[i])
  
  # load objects
  setwd(atac_dap_analysis_objects)
  atac.list[[i]] <- readRDS( rds_list[i] )
  
  setwd(atac_preprocessing_objects)
  atac.orig <- readRDS( rds_list2[i] )
  
  # if you run into this issue: https://github.com/stuart-lab/signac/issues/872 ,
  # manually set the object key. this fixed the issue.
  atac.list[[i]]@assays$recalled.peaks@key <- 'recalledpeaks_'
  
  Fragments(atac.list[[i]]) <- Fragments(atac.orig)
  
  # process data
  message('running tf-idf \n')
  atac.list[[i]] <- RunTFIDF(atac.list[[i]])
  message('finding variable features \n')
  atac.list[[i]] <- FindTopFeatures(atac.list[[i]], min.cutoff = 'q0')
  message('running svd \n')
  atac.list[[i]] <- RunSVD(atac.list[[i]])
}



# merge objects -----------------------------------------------------------

message('\nmerging objects \n')
atac.merged <- merge(x = atac.list[[1]], y = atac.list[2:length(atac.list)])



message('\nsaving... \n')
setwd(atac_dap_analysis_objects)
saveRDS(atac.merged, '6.process_new_peaks___srt_consensus_processed_merged.rds')
message('done!')
