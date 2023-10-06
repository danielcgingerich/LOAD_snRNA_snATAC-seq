cond <- grepl('/hpc_path_identifier', getwd())
if (cond){
  source('path1/to/config.R')
} else {
  source('path2/to/config.R')
} ; rm(cond)


for (i in 1:8){
  message(i, '\nloading data')
  setwd(atac_preprocessing_objects)
  atac <- list.files(pattern = '1.data_qc_.*_srt_preprocessed_unfiltered.rds')
  sample_id <- gsub(pattern = '1.data_qc___|_srt.*', replacement = '', atac)
  atac <- readRDS(atac[i])
  
  message('filtering cells...')
  max.frag <- quantile(atac$peak_region_fragments, probs = 0.99)
  cells.use <- atac@meta.data
  cells.use <- cells.use[cells.use$peak_region_fragments > 1000 & 
                           cells.use$peak_region_fragments < max.frag &
                           cells.use$pct_reads_in_peaks > 15 &
                           cells.use$blacklist_ratio < 0.05 &
                           cells.use$nucleosome_signal < 4 &
                           cells.use$TSS.enrichment > 2,]
  

  atac <- subset(atac, cells = rownames(cells.use))
  message("done! saving...")
  saveRDS(atac, paste0('2.apply_qc_filter___', sample_id[i], '_seurat_object_qc_filtered.rds'))
  rm(atac, max.frag, cells.use)
}


