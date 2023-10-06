cond <- grepl('/hpc_path_identifier', getwd())
if (cond){
  source('path1/to/config.R')
} else {
  source('path2/to/config.R')
} ; rm(cond)

# obtain merged peaks -----------------------------------------------------
setwd(atac_preprocessing_objects)
atac.list <- list.files(pattern = '^2.apply_qc_filter___')
atac.list <- lapply(atac.list, readRDS)
peaks <- lapply(atac.list, rownames)

peaks <- unlist(peaks)
peaks <- StringToGRanges(peaks)
combined_peaks <- reduce(peaks)

peakwidths <- width(combined_peaks)
combined_peaks <- combined_peaks[peakwidths < 10000 & peakwidths > 20]

setwd(atac_integration_objects)
saveRDS(combined_peaks, '1.create_common_peak_set___preak_granges.rds')


