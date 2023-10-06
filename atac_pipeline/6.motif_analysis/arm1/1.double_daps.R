cond <- grepl('/hpc_path_identifier', getwd())
if (cond){
  source('path1/to/config.R')
} else {
  source('path2/to/config.R')
} ; rm(cond)

library(data.table)
# import needed files -----------------------------------------------------




# -------------------------------------------------------------------------
findCCRES <- function(uni_or_mixed, fc_cut){
  setwd(atac_cicero_objects)
  
  regulatory_ccans <- readRDS('3.annotate_ccans___annotated.ccans.rds') %>% as.data.table() # only includes ccans with dap/deg overlaps
  regulatory_ccans <- regulatory_ccans[directionality == uni_or_mixed, ]
  
  regulatory_ccans$peak_ccan <- regulatory_ccans[, paste0(Peak, '_', ccan_cluster)]
  ccan_daps <- regulatory_ccans[p_val_adj <= 0.05, ] # all daps in any ccan

  pro_daps <- regulatory_ccans[p_val_adj <= 0.05, ] 
  pro_daps <- pro_daps[deg_overlap != '', ]
  pro_daps <- pro_daps[abs(deg_logfc) >= fc_cut,] 
  
  # overlapped deg and dap have same sign logFC
  pro_daps <- pro_daps[deg_logfc*avg_log2FC > 0, ]
  
  cluster <- unique(ccan_daps$cluster)
  conns <- lapply(cluster, 
                  function(x){
                    setwd(atac_cicero_objects)
                    message(x)
                    # both peaks are daps
                    tmp <- ccan_daps[cluster == x, ]
                    conns_x <- paste0('2.run_cicero___', x, '_cicero_conns.rds')
                    conns_x <- readRDS(conns_x) %>% as.data.table()
                    conns_x$Peak2 <- as.character(conns_x$Peak2)
                    conns_x <- conns_x[!is.na(coaccess), ]
                    conns_x <- conns_x[Peak1 %in% tmp$Peak, ]
                    conns_x <- conns_x[Peak2 %in% tmp$Peak, ]
                    
                    tmp <- tmp[, .(Peak, peak_ccan, ccan_cluster)]
                    conns_x <- left_join(x = conns_x, y = tmp, by = c('Peak1' = 'Peak'))
                    conns_x <- left_join(x = conns_x, y = tmp, by = c('Peak2' = 'Peak'))
                    colnames(conns_x) <- gsub('.x', '1', colnames(conns_x))
                    colnames(conns_x) <- gsub('.y', '2', colnames(conns_x))
                    conns_x <- conns_x[ccan_cluster1 == ccan_cluster2, ]
                    return(conns_x)
                  })
  conns <- conns[sapply(conns, nrow) > 0] %>% bind_rows()
  conns <- conns[coaccess >= 0.2,]

  # one of the peaks overlaps a deg
  conns <- conns[peak_ccan1 %in% pro_daps$peak_ccan | peak_ccan2 %in% pro_daps$peak_ccan, ]
  
  # merge with rest of data
  conns <- left_join(x = conns, y = regulatory_ccans[, .(peak_ccan, avg_log2FC)], by = c('peak_ccan1' = 'peak_ccan'))
  conns <- left_join(x = conns, y = regulatory_ccans[, .(peak_ccan, avg_log2FC)], by = c('peak_ccan2' = 'peak_ccan'))

  # both daps go same direction
  conns <- conns[avg_log2FC.x*avg_log2FC.y > 0, ]
  
  cherrypicker <- unique(c(conns$peak_ccan1, conns$peak_ccan2))
  ccres <- regulatory_ccans[peak_ccan %in% cherrypicker, ]
  
  # n degs 
  n_degs <- ccres[deg_overlap != '', unique(deg_overlap), by = cluster][, .N, by = cluster]
  # names
  deg_names <- ccres[deg_overlap != '', paste(unique(deg_overlap), collapse = ', '), by = cluster]
  # n ccans
  n_ccans <- ccres[deg_overlap != '', unique(ccan_cluster), by = cluster][, .N, by = cluster]
  
  ccres <- full_join(x = n_ccans, y = n_degs, by = 'cluster')
  ccres <- full_join(x = ccres, y = deg_names, by = 'cluster')
  colnames(ccres) <- c('cluster', 'n_ccans', 'n_degs', 'deg_names')
  ccres$degs <- paste0(ccres$n_degs, '\n', ccres$deg_names)
  ccres <- ccres[, .(cluster, n_ccans, degs)]
  return(ccres)
}


fc <- 0.15

uni <- findCCRES(uni_or_mixed = 'unidirectional', fc_cut = fc)
mix <- findCCRES(uni_or_mixed = 'mixed', fc_cut = fc)

setwd(atac_motif_enrichment_objects)
write.csv(uni, '1.double_daps___unidirectional_ccans.csv')
write.csv(mix, '1.double_daps___mixed_ccans.csv')
