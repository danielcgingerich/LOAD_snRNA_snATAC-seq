cond <- grepl('/hpc_path_identifier', getwd())
if (cond){
  source('path1/to/config.R')
} else {
  source('path2/to/config.R')
} ; rm(cond)

library(GenomicRanges)
library(SNPlocs.Hsapiens.dbSNP151.GRCh38)
library(MafDb.gnomAD.r2.1.GRCh38)
library(data.table)
# use maf cutoff of 5%


# -------------------------------------------------------------------------
message('loading motif info')
setwd(atac_motif_enrichment_objects)
motif_info <- readRDS('2.pwm_info___motif.info.rds')

setwd(atac_cicero_objects)
kunkle_ccans <- readRDS('3.annotate_ccans___annotated.ccans.rds')
gwas_ccan_names <- kunkle_ccans[grep('kunkle', kunkle_ccans$snp_overlap), 'ccan_cluster']
kunkle_ccans <- kunkle_ccans[kunkle_ccans$ccan_cluster %in% gwas_ccan_names, ]

hg38 <- rtracklayer::import(paste0(ref_genome_path, 'gencode.v32.primary_assembly.annotation.gtf'))

# find highly coaccessible dap-degs -----------------------------------------

# function extracts daps and degs from each ccan and finds dap-deg pairs that
# are highly coaccessible and have same sign log FC
findPeakPairs <- function(directionality, coaccess_cutoff){
  atsnp_daps <- list()
  ccans0 <- kunkle_ccans[kunkle_ccans$directionality == directionality, ]
  cluster <- unique(ccans0$cluster)
  message(length(cluster), ' total clusters')
  
  for (i in cluster){
    message(i)
    ccan_subset <- ccans0[ccans0$cluster == i, ]
    
    daps <- ccan_subset$Peak[ccan_subset$p_val_adj <= 0.05]
    degs <- ccan_subset$Peak[ccan_subset$deg_overlap != '']
    
    setwd(atac_cicero_objects)
    conns <- readRDS(paste0('2.run_cicero___', i, '_cicero_conns.rds'))
    conns <- conns[!is.na(conns$coaccess), ]
    conns$Peak2 <- as.character(conns$Peak2)
    
    # cicero returns all reverse pairs of peaks, so only one search is sufficient
    # dont need to look for daps in peak2 and degs in peak1
    ix <- conns$Peak1 %in% daps & conns$Peak2 %in% degs
    target_conns <- conns[ix, ]
    
    target_conns <- target_conns[target_conns$coaccess >= coaccess_cutoff, ]
    
    colnames(target_conns) <- c('dap_peak', 'deg_peak', 'coaccess')
    
    target_peaks <- c(target_conns$dap_peak, target_conns$deg_peak)
    tmp <- ccan_subset[ccan_subset$Peak %in% target_peaks, ]
    
    target_conns <- left_join(x = target_conns, 
                              y = tmp[, c('Peak', 'avg_log2FC', 'CCAN')],
                              by = c('dap_peak' = 'Peak'))
    target_conns <- left_join(x = target_conns, 
                              y = tmp[, c('Peak', 'deg_logfc', 'CCAN')], 
                              by = c('deg_peak' = 'Peak'))
    ix <- target_conns$avg_log2FC*target_conns$deg_logfc >= 0 & target_conns$CCAN.x == target_conns$CCAN.y
    
    target_conns <- target_conns[ix, ]
    atsnp_daps[[i]] <- target_conns
  }
  return(atsnp_daps)
}

atsnp_daps_uni <- findPeakPairs(directionality = 'unidirectional', coaccess_cutoff = 0.2)
atsnp_daps_mixed <- findPeakPairs(directionality = 'mixed', coaccess_cutoff = 0.2)

atsnp_daps_uni <- atsnp_daps_uni[sapply(atsnp_daps_uni, nrow) > 0]
atsnp_daps_mixed <- atsnp_daps_mixed[sapply(atsnp_daps_mixed, nrow) > 0]

setwd(atac_motif_enrichment_objects)
saveRDS(atsnp_daps_uni, '1.atsnp_target_dap_deg_pairs___unidirectional.rds')
saveRDS(atsnp_daps_mixed, '1.atsnp_target_dap_deg_pairs___mixed.rds')
