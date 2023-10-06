cond <- grepl('/hpc_path_identifier', getwd())
if (cond){
  source('path1/to/config.R')
} else {
  source('path2/to/config.R')
} ; rm(cond)


slurm.array.id <- Sys.getenv('SLURM_ARRAY_TASK_ID')
i <- as.numeric(slurm.array.id)


setwd(atac_cell_type_annotation_objects)
atac <- list.files(pattern = '1.gene.activity___.*rds$')[i]
sample_id <- gsub('.*___|_srt_.*', '', atac)
atac <- readRDS(atac)
rna <- readRDS('2.snrna_lognormalized_integration___data.rds')



# set default assays ------------------------------------------------------
DefaultAssay(atac) <- 'gene.activity'
DefaultAssay(rna) <- 'integrated'



# normalize gene activities -----------------------------------------------
atac <- NormalizeData(atac, assay = 'gene.activity', 
                      normalization.method = 'LogNormalize', 
                      scale.factor = median(atac$nCount_gene.activity))
atac <- FindVariableFeatures(atac, nfeatures = 5000)



# find transfer anchors ---------------------------------------------------
message('\nfinding anchors \n')
transfer.anchors <- FindTransferAnchors(reference = rna, 
                                        query = atac,
                                        reduction = 'cca',
                                        features = rownames(rna))



# transfer data -----------------------------------------------------------
message('\ntransferring metadata \n')
cell_type_predictions <- TransferData(anchorset = transfer.anchors,
                                 refdata = rna$cell_type,
                                 weight.reduction = atac[['lsi']],
                                 dims = 2:30)
cluster_predictions <- TransferData(anchorset = transfer.anchors,
                                 refdata = rna$cell_type_cluster,
                                 weight.reduction = atac[['lsi']],
                                 dims = 2:30)

message('\ncomputing weight mtx \n')
weight.mtx <-  TransferData(anchorset = transfer.anchors,
                            refdata = rna$cell_type_cluster,
                            query = atac,
                            weight.reduction = atac[['lsi']],
                            dims = 2:30, store.weights = T)
weight.mtx <- weight.mtx@tools$TransferData$weights.matrix
colnames(weight.mtx) <- gsub('_query', '', transfer.anchors@query.cells)
rna.index <- transfer.anchors@anchors[,1]
rownames(weight.mtx) <- gsub('_reference', '', transfer.anchors@reference.cells[rna.index])




# save files -------------------------------------------------------------
setwd(atac_cell_type_annotation_objects)
saveRDS(cell_type_predictions, paste0('3.annotate_atac_cell_types___', sample_id, '_cell_type_predictions.rds'))
saveRDS(cluster_predictions, paste0('3.annotate_atac_cell_types___', sample_id, '_cluster_predictions.rds'))
saveRDS(weight.mtx, paste0('3.annotate_atac_cell_types___', sample_id, '_anchor_weight_matrix.rds'))


