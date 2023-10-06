cond <- grepl('/hpc_path_identifier', getwd())
if (cond){
  source('path1/to/config.R')
} else {
  source('path2/to/config.R')
} ; rm(cond)

setwd(atac_integration_objects)
atac_list <- list.files(pattern = '2.create_combined_peak_assay___')
sample_id <- gsub('.*___|.rds', '', atac_list)

lapply(1:length(atac_list),
       function(x){
         message('processing...')
         obj <- readRDS(atac_list[x])
         obj_id <- sample_id[x]
         obj@assays$combined_peaks@key <- 'combinedpeaks_'
         obj <- RunTFIDF(obj)
         obj <- FindTopFeatures(obj, min.cutoff = 'q0')
         obj <- RunSVD(obj)
         message('saving...')
         saveRDS(obj, paste0('3.re-process_and_merge___', obj_id, '_combined_peak_assay.rds'))
       })

atac_list <- list.files(pattern = '3.re-process_and_merge___')
sample_id <- gsub('.*___|_combined_.*', '', atac_list)
atac_list <- lapply(atac_list, readRDS)
atac_merged <- merge(x = atac_list[[1]], y = atac_list[2:length(atac_list)])


setwd(atac_integration_objects)
saveRDS(atac_merged, '3.re-process_and_merge___merged_dataset.rds')

