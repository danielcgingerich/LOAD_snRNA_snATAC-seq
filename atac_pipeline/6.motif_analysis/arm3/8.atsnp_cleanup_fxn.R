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


# set up data ------------------------------------------------------------
###
maf = 0.01




###
getResults <- function(maf = NULL){
  setwd(atac_motif_enrichment_objects)
  motifs <- readRDS('7.atsnp_candidate_motifs___all_motifs_passed_detection_threshold.rds')
  
  message('retrieving atSNP data')
  atsnp_scores <- readRDS('3.compute_motif_score___atsnp_motif_scores.rds')
  atsnp_scores <- as.data.table(atsnp_scores$motif.scores)
  atsnp_results <- cbind(atsnp_scores, readRDS('6.wrangle_data___pval_stats.rds'))
  
 
  
  
  atsnp_results$id <- atsnp_results[, paste(motif, snpid, snpbase, sep = '_')]
  motifs$id <- motifs[, paste(motif, rsid, snp_allele, sep = '_')]
  atsnp_results <- atsnp_results[id %in% motifs$id, ] 
  # (above) motifmatchr uses exact p values. this is better for matching 
  # than atsnps approximate p values. there are slight misalignments, but
  # this occurs very little and is due to different background
  # models and PWM scoring methods.  (motifbreakr assumes each base is iid multinomial) 
  # (atsnp assumes background is first order markov)

  message('calculating FDRs')
  atsnp_results$fdr_diff <- atsnp_results[, p.adjust(pval_diff, method = 'fdr')]
  atsnp_results$fdr_rank <- atsnp_results[, p.adjust(pval_rank, method = 'fdr')]

  if ( !is.null(maf) ){
    message('filtering out SNPs with MAF < ', maf*100, '%')
    snps <- readRDS('2.atsnp_target_snps_and_maf___atsnp_target_snps.rds')
    snps <- mcols(snps) %>% as.data.table()
    keep <- snps[gnomad_maf >= maf | bravo_maf >= maf, RefSNP_id]
    message(length(unique(keep)), ' SNPs retained')
    atsnp_results <- atsnp_results[snpid %in% keep]
  }
  
  # regions of interest
  message('finding regions of interest')
  setwd(atac_motif_enrichment_objects)
  loadRegions <- function(file){
    directionality <- gsub('.*___|.rds', '', file)
    regions <- readRDS(file)
    regions <- lapply(1:length(regions), 
                      function(x){
                        regions[[x]]$cluster <- names(regions)[x]
                        regions[[x]]$directionality <- directionality
                        return(regions[[x]])
                      })
    regions <- bind_rows(regions)
    return(regions)
  }
  regions <- list.files(pattern = '1.atsnp_target_dap_deg_pairs.*rds$')
  regions <- lapply(regions, loadRegions)
  regions <- bind_rows(regions)
  tmp <- regions
  regions <- StringToGRanges(regions$dap_peak)
  mcols(regions) <- tmp ; rm(tmp)
  
  
  
  
  
  
  
  
  
  snp_locs <- readRDS('2.atsnp_target_snps_and_maf___atsnp_target_snps.rds')
  tmp <- paste0('chr', as.character(seqnames(snp_locs)), 
                '-', start(snp_locs), '-', start(snp_locs))
  names(tmp) <- snp_locs$RefSNP_id
  atsnp_results$snp_position <- tmp[atsnp_results$snpid]
  
  atsnp_positions <- StringToGRanges(atsnp_results$snp_position)
  atsnp_overlaps <- findOverlaps(query = regions, subject = atsnp_positions)
  dap_hits <- queryHits(atsnp_overlaps)
  atsnp_hits <- subjectHits(atsnp_overlaps)
  atsnp_results_filtered <- list(atsnp_results[atsnp_hits, ], 
                                 as.data.table(mcols(regions)[dap_hits, ]))
  atsnp_results_filtered <- bind_cols(atsnp_results_filtered)
  
  
  
  message('filtering out motifs not expressed')
  tf_exprs <- readRDS('6.expressed_enriched_tfs___tfs.expressed.in.clusters.rds')
  tf_exprs <- lapply(1:length(tf_exprs), 
                     function(x){
                       df <- tf_exprs[[x]]
                       df$cluster <- names(tf_exprs)[x]
                       df
                     }) %>% bind_rows() %>% as.data.table()
  tf_exprs$id <- tf_exprs[, paste0(id, '_', cluster)]
  tf_exprs <- tf_exprs[!(duplicated(id)), ]
  atsnp_results_filtered$id <- atsnp_results_filtered[, paste0(motif, '_', cluster)]
  atsnp_results_filtered <- left_join(x = atsnp_results_filtered, 
                                      y = tf_exprs[, .(id, tf_name = gene)], 
                                      by = 'id')
  atsnp_results_filtered$tf_expressed <- atsnp_results_filtered[, ifelse(is.na(tf_name), 0, 1)]
  
  message('filtering out motifs not enriched')
  homer <- list.files(pattern = '6.expressed_enriched_tfs___homer.*enriched.tfs.rds$')
  homer <- lapply(homer, 
                  function(x){
                    df <- readRDS(x)
                    df <- lapply(df, as.data.table)
                    d <- gsub('.*___homer.|.enriched.tfs.rds', '', x)
                    df <- lapply(names(df), 
                                 function(y){
                                   
                                   if (nrow(df[[y]]) > 0){
                                     df[[y]]$cluster <- y
                                     df[[y]]
                                   }
                                   
                                 })
                    df <- bind_rows(df)
                    
                    df$directionality <- d
                    return(df)
                  })
  homer <- bind_rows(homer)
  
  atsnp_results_filtered$id <- atsnp_results_filtered[, paste(motif, cluster, directionality)]
  homer$id <- homer[, paste(motif, cluster, directionality)]
  homer <- homer[!(duplicated(id)), ]
  atsnp_results_filtered <- left_join(x = atsnp_results_filtered, 
                                      y = homer[, .(fold_enrichment, id)], by = 'id')
  
  ix <- is.na(atsnp_results_filtered$fold_enrichment)
  atsnp_results_filtered$fold_enrichment[ix] <- 0
  
  
  
  
  
  
  # setwd(motif_enrichment_objects)
  # message('filtering out degs < ', deg_filter, ' logFC')
  # ccans <- readRDS('set2_annotated.ccans.rds')
  # ccans <- as.data.table(ccans)
  # ccans <- ccans[!is.na(deg.logfc)]
  # ccans <- ccans[abs(deg.logfc) >= deg_filter, ]
  # ccans$id <- ccans[, paste0(Peak, '_', cluster)]
  # 
  # atsnp_results_filtered$id <- atsnp_results_filtered[, paste0(deg_peak, '_', cluster)]
  # atsnp_results_filtered <- atsnp_results_filtered[id %in% ccans$id, ]
  
  
  message('adding annotations')
  setwd(atac_cicero_objects)
  annotated_ccans <- readRDS('3.annotate_ccans___annotated.ccans.rds')
  annotated_ccans <- as.data.table(annotated_ccans)
  
  atsnp_results_filtered$id <- atsnp_results_filtered[, paste(deg_peak, cluster, sep = '_')]
  annotated_ccans$id <- annotated_ccans[, paste(Peak, cluster, sep = '_')]
  
  atsnp_results_filtered <- left_join(x = atsnp_results_filtered, y = annotated_ccans[, .(deg_overlap, deg_overlap_type, id)], 
                                      by = 'id')
  
  gwas_snps <- annotated_ccans[grep('kunkle', snp_overlap), .(CCAN, cluster, snp_overlap)]
  gwas_snps$snp_overlap <- gsub('kunkle:|,.*', '', gwas_snps$snp_overlap)
  gwas_snps <- unique(gwas_snps)
  
  gwas_snps$id <- gwas_snps[, paste(CCAN, cluster, sep = '_')]
  atsnp_results_filtered$id <- atsnp_results_filtered[, paste(CCAN.x, cluster, sep = '_')]
  
  atsnp_results_filtered <- left_join(x = atsnp_results_filtered, y = gwas_snps[, .(snp_overlap, id)], by = 'id')
  
  motifs$id <- motifs[, paste(motif, rsid, snp_allele, sep = '_')]
  atsnp_results_filtered$id <- atsnp_results_filtered[, paste(motif, snpid, snpbase, sep = '_')]
  
  atsnp_results_filtered <- left_join(x = atsnp_results_filtered, 
                                      y = motifs[, .(id, tfbs_position_ref, 
                                                     tfbs_position_snp, 
                                                     pwm_score_ref, 
                                                     pwm_score_snp)], 
                                      by = 'id')
  
  setwd(atac_motif_enrichment_objects)
  df <- readRDS("2.pwm_info___motif.info.rds")
  tf_names <- df$tf_name
  names(tf_names) <- df$id
  atsnp_results_filtered$tf_name <- tf_names[atsnp_results_filtered$motif]
  
  atsnp_results_filtered <- atsnp_results_filtered[, .(cluster, gwas_locus = snp_overlap, ccan = CCAN.x, directionality,
                                                       motif, tf_name, 
                                                       rsid = snpid, snp_position,
                                                       snp_base = snpbase, tfbs_position_ref,
                                                       tfbs_position_snp, affinity_change = -1*log_lik_ratio,
                                                       dap = dap_peak, deg_peak,
                                                       deg = paste0(deg_overlap, '(', deg_overlap_type, ')'),
                                                       coaccess, dap_logfc = avg_log2FC,
                                                       deg_logfc = deg_logfc, tf_expressed, fold_enrichment,
                                                       pval_ref, pval_snp, 
                                                       pval_rank, fdr_rank, sigma_rank, min_rank, max_rank,
                                                       pval_diff, fdr_diff, sigma_diff, min_diff, max_diff)]
  
  return(atsnp_results_filtered)
}
