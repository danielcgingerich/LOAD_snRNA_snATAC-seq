cond <- grepl('/hpc_path_identifier', getwd())
if (cond){
  source('path1/to/config.R')
} else {
  source('path2/to/config.R')
} ; rm(cond)

library(GenomicRanges)

setwd(atac_cicero_objects)
all_ccans <- readRDS('3.annotate_ccans___annotated.ccans.rds')

setwd(atac_motif_enrichment_objects)
motif_info <- readRDS('2.pwm_info___motif.info.rds')

setwd(rna_cell_type_annotation_objects)
rna <- readRDS('3.pca_lovain_cluster_umap___annotated_integrated_processed_dataset.rds')
metadata <- rna$cell_type_cluster
counts <- rna@assays$RNA@counts


# expression filter -------------------------------------------------------

cluster <- unique(all_ccans$cluster)
tf_exp <- list()
for (i in 1:length(cluster)){
  message(i)
  cells.use <- names(metadata)[metadata == cluster[i] ]
  genes.use <- unique(motif_info$gene)
  mtx <- counts[genes.use, cells.use]
  pct_exp <- rowSums(mtx>0)/length(cells.use)
  pct_exp <- pct_exp[pct_exp>=0.1]
  gene <- names(pct_exp)
  ix <- motif_info$gene %in% gene
  tf_exp[[i]] <- motif_info[ix, ]
}
names(tf_exp) <- cluster

setwd(atac_motif_enrichment_objects)
saveRDS(tf_exp, '6.expressed_enriched_tfs___tfs.expressed.in.clusters.rds')



# enrichment filter -------------------------------------------------------
setwd(atac_motif_enrichment_objects)
homer_uni <- list.files(pattern = '5.clean_homer_results___.*_unidirectional.rds')
homer_mixed <- list.files(pattern = '5.clean_homer_results___.*_mixed.rds')

uni_clusters <- gsub('5.clean_homer_results___|_homer.*', '', homer_uni)
mixed_clusters <- gsub('5.clean_homer_results___|_homer.*', '', homer_mixed)

names(homer_uni) <- uni_clusters
names(homer_mixed) <- mixed_clusters

fe <- 1.2 # how much fold enrichment should you use?

homer_uni <- lapply(1:length(homer_uni),
                    function(x){
                      results <- homer_uni[x]
                      results <- readRDS(results)
                      cluster <- names(homer_uni)[x]
                      exp <- tf_exp[[cluster]]
                      ix <- results$motif %in% exp$id  & results$fold_enrichment >= fe
                      results[ix, ]
                    })
names(homer_uni) <- uni_clusters

homer_mixed <- lapply(1:length(homer_mixed),
                    function(x){
                      results <- homer_mixed[x]
                      results <- readRDS(results)
                      cluster <- names(homer_mixed)[x]
                      exp <- tf_exp[[cluster]]
                      ix <- results$motif %in% exp$id & results$fold_enrichment >= fe
                      results[ix, ]
                    })
names(homer_mixed) <- mixed_clusters

setwd(atac_motif_enrichment_objects)
saveRDS(homer_uni, '6.expressed_enriched_tfs___homer.uni.enriched.tfs.rds')
saveRDS(homer_mixed, '6.expressed_enriched_tfs___homer.mixed.enriched.tfs.rds')
