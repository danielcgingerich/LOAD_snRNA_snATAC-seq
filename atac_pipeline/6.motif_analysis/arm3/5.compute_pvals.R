cond <- grepl('/hpc_path_identifier', getwd())
if (cond){
  source('path1/to/config.R')
} else {
  source('path2/to/config.R')
} ; rm(cond)

library(SNPlocs.Hsapiens.dbSNP151.GRCh38)
library(data.table)
library(parallel)
library(BSgenome.Hsapiens.UCSC.hg38)
library(JASPAR2020)
library(atSNP)

i <- Sys.getenv('SLURM_ARRAY_TASK_ID') %>% as.numeric()

# lets use cutoff of 5%. can easily get cutoff recomendations on google.

setwd(atac_motif_enrichment_objects)
snp_info <- readRDS('3.compute_motif_score___snp_info.rds')
pwms <- readRDS('3.compute_motif_score___atsnp_ppms.rds')
atsnp.scores <- readRDS('3.compute_motif_score___atsnp_motif_scores.rds')

generateResults <- function(dummy_argument){ # dont actually use any argument. lapply needs one though. 
  message('computing p values')
  atsnp_results <- ComputePValues(motif.lib = pwms, snp.info = snp_info,
                                  motif.scores = atsnp.scores$motif.scores,
                                  testing.mc=FALSE)
  
  message('done. converting format')
  atsnp_results <- as.data.table(atsnp_results)
  atsnp_results <- atsnp_results[, .(pval_ref, pval_snp, pval_diff, pval_rank)]
  return(atsnp_results)
}

atsnp_sum <- generateResults('asdf')

setwd(atac_motif_enrichment_objects)
saveRDS(atsnp_sum, paste0('5.compute_pvals___atsnp', i, '.rds'))
message('done running atsnp')

