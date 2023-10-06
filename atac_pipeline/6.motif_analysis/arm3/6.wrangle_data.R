cond <- grepl('/hpc_path_identifier', getwd())
if (cond){
  source('path1/to/config.R')
} else {
  source('path2/to/config.R')
} ; rm(cond)
library(data.table)
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg38)
library(Biostrings)
library(SNPlocs.Hsapiens.dbSNP151.GRCh38)


setwd(atac_motif_enrichment_objects)

results <- list.files(pattern = '5.compute_pvals___atsnp.*.rds$')

results <- lapply(results, 
                  readRDS) %>% bind_cols()

results <- as.matrix(results)

p <- p_sigma <- p_min <- p_max <- list()

for (i in 1:4){
  ix <- seq(i, ncol(results), 4)
  df <- results[, ix]
  sigma <- apply(df, MARGIN = 1, var)
  mu <- apply(df, MARGIN = 1, mean)
  q0 <- apply(df, MARGIN = 1, min)
  q1 <- apply(df, MARGIN = 1, max)
  
  p[[i]] <- mu
  p_sigma[[i]] <- sigma
  p_min[[i]] <- q0
  p_max[[i]] <- q1
}

p <- bind_cols(p)
p_sigma <- bind_cols(p_sigma)
p_min <- bind_cols(p_min)
p_max <- bind_cols(p_max)

colnames(p) <- c('pval_ref', 'pval_snp', 'pval_diff', 'pval_rank')
colnames(p_sigma) <- c('sigma_ref', 'sigma_snp', 'sigma_diff', 'sigma_rank')
colnames(p_min) <- c('min_ref', 'min_snp', 'min_diff', 'min_rank')
colnames(p_max) <- c('max_ref', 'max_snp', 'max_diff', 'max_rank')

pvals <- cbind(p, p_sigma, p_min, p_max)
pvals <- as.data.table(pvals)

pvals <- pvals[, .(pval_ref, sigma_ref, min_ref, max_ref, 
                   pval_snp, sigma_snp, min_snp, max_snp, 
                   pval_rank, sigma_rank, min_rank, max_rank,
                   pval_diff, sigma_diff, min_diff, max_diff)]

setwd(atac_motif_enrichment_objects)
saveRDS(pvals, '6.wrangle_data___pval_stats.rds')