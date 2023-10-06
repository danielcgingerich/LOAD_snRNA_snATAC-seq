cond <- grepl('/hpc_path_identifier', getwd())
if (cond){
  source('path1/to/config.R')
} else {
  source('path2/to/config.R')
} ; rm(cond)

library(Matrix.utils)


setwd(atac_dap_analysis_objects)
atac <- readRDS('6.process_new_peaks___srt_consensus_processed_merged.rds')
mtx <- atac@assays$recalled.peaks@data
metadata <- readRDS('7.add_covarates___metadata_for_dap_analysis.rds')

metadata$cluster_links_diagnosis <- paste0(metadata$cluster_links, '_', metadata$diagnosis)
n = length(unique(metadata$cluster_links))
for (i in 1:n){
  message(i, ' of ', n)
  message('Preping data for:')
  cluster <- levels(as.factor(metadata$cluster_links))[i]
  message(cluster, '\n')
  
  load_cells <- rownames(metadata[metadata$cluster_links_diagnosis == paste0(cluster, '_LOAD'), ])
  norm_cells <- rownames(metadata[metadata$cluster_links_diagnosis == paste0(cluster, '_Normal'), ])
  message('total LOAD: ', length(load_cells))
  message('total Normal: ', length(norm_cells))
  
  message('calculating percent exprs...')
  pct_load <- round(rowSums(mtx[, load_cells] > 0)/length(load_cells), digits = 3)
  pct_norm <- round(rowSums(mtx[, norm_cells] > 0)/length(norm_cells), digits = 3)
  
  message('calculating log2FC...')
  u_load <- log(x = rowMeans(mtx[, load_cells]) + 1, base = 2)
  u_norm <- log(x = rowMeans(mtx[, norm_cells]) + 1, base = 2)
  logfc <- u_load - u_norm
  
  fc_results <- data.frame('pct.1' = pct_load, 'pct.2' = pct_norm, 'avg_log2FC' = logfc)
  
  message('saving...')
  setwd(atac_dap_analysis_objects)
  saveRDS(fc_results, paste0('8.calculate_fold_change___', cluster, '.fc_results.rds'))
  message('done! \n')
}
