cond <- grepl('/hpc_path_identifier', getwd())
if (cond){
  source('path1/to/config.R')
} else {
  source('path2/to/config.R')
} ; rm(cond)

# import needed files -----------------------------------------------------

# load rna gene exprs
setwd(rna_cell_type_annotation_objects)
rna <- readRDS('3.pca_lovain_cluster_umap___annotated_integrated_processed_dataset.rds')
metadata <- rna$cell_type_cluster
counts <- rna@assays$RNA@counts

# motif info 
setwd(atac_motif_enrichment_objects)
motif_info <- readRDS('2.pwm_info___motif.info.rds')

cleanHomer <- function(directionality){
  message('loading homer file')
  setwd(atac_motif_enrichment_objects)
  setwd('homer_outs')
  files <- list.files(pattern = directionality)
  cluster <- gsub('_.*', '', files)
  for (i in 1:length(files)){
    message(i)
    # load motif results
    setwd(atac_motif_enrichment_objects)
    setwd('homer_outs')
    
    setwd(files[i])
    homer <- readLines('knownResults.txt')
    homer[1] <- tolower(homer[1])
    homer <- strsplit(homer, split = '\t')
    homer <- lapply(homer, unlist)
    header <- homer[[1]]
    homer <- homer[-1]
    homer <- lapply(homer, as.matrix)
    homer <- lapply(homer, t)
    homer <- lapply(homer, as.data.frame)
    homer <- bind_rows(homer)
    colnames(homer) <- c('motif', 'consensus', 'p_val', 'log_p', 'fdr',
                         'n_tgt', 'pct_tgt', 'n_bg', 'pct_bg')
    homer$p_val <- as.numeric(homer$p_val)
    homer$log_p <- as.numeric(homer$log_p)
    homer$fdr <- as.numeric(homer$fdr)
    homer$pct_tgt <- gsub('%', '', homer$pct_tgt) %>% as.numeric()
    homer$pct_bg <- gsub('%', '', homer$pct_bg) %>% as.numeric()
    
    ### homer finding 0 matches for some motifs. 
    # i checked this with motifmatchr and it is accurate
    # remove unmatched motifs
    homer <- homer[homer$pct_bg != 0, ] 
    homer$fold_enrichment <- homer$pct_tgt/homer$pct_bg
    
    # motif name
    motif_info_tmp <- unique(motif_info[, c('id', 'tf_name')])
    homer <- left_join(x = homer, y = motif_info_tmp, by = c('motif' = 'id'))
    
    # motif is deg 
    message('loading degs')
    setwd(rna_deg_analysis_objects)
    degs <- paste0('1.mast___', cluster[i], '_model_results.csv')
    degs <- read.csv(degs, header = T, row.names = 1)
    degs <- degs[degs$fdr <= 0.05, ]
    
    message('finding degs')
    enriched_motif_info <- motif_info[motif_info$id %in% homer$motif,]
    degs <- degs[, c('primerid', 'coef', 'fdr')]
    enriched_motif_info <- left_join(x = enriched_motif_info, y = degs, by = c('gene' = 'primerid'))
    de_motif <- enriched_motif_info[!is.na(enriched_motif_info$coef), ]
    colnames(de_motif) <- c('id', 'deg','tf_name', 'coef', 'fdr_deg')
    ix <- which(colnames(de_motif) == 'tf_name')
    de_motif <- de_motif[, -ix] # homer already has this column
    
    if (nrow(de_motif) > 0){
      homer <- left_join(x = homer, y = de_motif, by = c('motif' = 'id'))
      colnames(homer)[5] <- 'fdr_motif'
    } else {
      homer$deg <- NA
      homer$coef <- NA
      homer$fdr_deg <- NA
      colnames(homer)[5] <- 'fdr_motif'
    }
    
    # isolate sig results 
    message('filtering p values')
    homer <- homer[pmax(homer$p_val, homer$fdr_motif)<=0.05, ]
    
    # motif expressed at least 10% of cells 
    message('filtering pct exprs')
    enriched_motif_info <- motif_info[motif_info$id %in% homer$motif,]
    cells.use <- names(metadata)[metadata == cluster[i]]
    genes.use <- unique(enriched_motif_info$gene)
    mtx <- counts[genes.use, cells.use, drop = F]
    pct_exprs <- rowSums(mtx > 0)/length(cells.use)
    expressed_motifs <- names(pct_exprs[pct_exprs >= 0.1])
    expressed_motifs <- enriched_motif_info[enriched_motif_info$gene %in% expressed_motifs,]
    homer <- homer[homer$motif %in% expressed_motifs$id, ]
    
    # save results 
    message('saving')
    setwd(atac_motif_enrichment_objects)
    saveRDS(homer, paste0('5.clean_homer_results___', cluster[i], '_homer.results.cleaned_', directionality, '.rds'))
    message('done!\n')
  }
}

cleanHomer(directionality = 'unidirectional')
cleanHomer(directionality = 'mixed')



