cond <- grepl('/hpc_path_identifier', getwd())
if (cond){
  source('path1/to/config.R')
} else {
  source('path2/to/config.R')
} ; rm(cond)

# i <- Sys.getenv('SLURM_ARRAY_TASK_ID') %>% as.numeric()

message('loading...')
setwd(atac_cluster_link_validation_objects)
atac <- readRDS('1.preprocessing_rna_atac___atac_preprocessed.rds')


# rna <- readRDS('1.preprocessing_rna_atac___rna_preprocessed.rds')
# features <- rownames(rna)
# setwd(atac_cluster_link_validation_objects)
# saveRDS(features, '2.gene_activity___gene_names_for_gene_activity.rds')

features <- readRDS('2.gene_activity___gene_names_for_gene_activity.rds')

Annotation(atac)$gene_biotype <- 'protein_coding'
activity_mtx <- GeneActivity(atac, features = features, max.width = Inf)

setwd(atac_cluster_link_validation_objects)
saveRDS(activity_mtx, paste0('2.gene_activity__gene_activity_matrix.rds'))

