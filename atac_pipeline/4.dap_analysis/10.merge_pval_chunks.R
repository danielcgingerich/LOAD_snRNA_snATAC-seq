cond <- grepl('/hpc_path_identifier', getwd())
if (cond){
  source('path1/to/config.R')
} else {
  source('path2/to/config.R')
} ; rm(cond)

slurm.array.id <- Sys.getenv('SLURM_ARRAY_TASK_ID')
i <- as.numeric(slurm.array.id)

setwd(atac_dap_analysis_objects)

metadata <- readRDS('7.add_covarates___metadata_for_dap_analysis.rds')

cluster <- levels(as.factor(metadata$cluster_links))[i]
peaks <- readRDS('4.combine_consensus_peaks___cluster_specific_consensus_peaks_c=2.rds')
peaks <- peaks[grep(cluster, peaks$peak.called.in),]

n <- list.files(pattern = paste0(cluster, '_pvals')) %>% length
message('number of partitions: ', n)
p_val <- paste0('9.likelihood_ratio___', cluster, '_pvals', c(1:n), '.rds') 
p_val <- unlist(lapply(p_val, readRDS))
names(p_val) <- GRangesToString(peaks)

fc <- readRDS(paste0('8.calculate_fold_change___', cluster, '.fc.results.rds'))
fc <- fc[names(p_val),]
daps <- cbind(fc, as.data.frame(p_val))
daps$p_val_adj <- p.adjust(p = daps$p_val, method = 'fdr', n = nrow(daps))
daps <- daps[order(daps$p_val),]


# save data ---------------------------------------------------------------
setwd(atac_dap_analysis_objects)
saveRDS(daps, paste0('10.merge_pval_chunks___', cluster, '_fast_daps.rds'))





