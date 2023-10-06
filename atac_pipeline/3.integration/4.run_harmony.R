cond <- grepl('/hpc_path_identifier', getwd())
if (cond){
  source('path1/to/config.R')
} else {
  source('path2/to/config.R')
} ; rm(cond)


# load objects ------------------------------------------------------------
setwd(atac_integration_objects)
atac_merged <- readRDS('3.re-process_and_merge___merged_dataset.rds')
setwd(atac_cell_type_annotation_objects)
pred_id <- list.files(pattern = '3.ann.*cell_type_predictions.rds$')
sample_id <- gsub(pattern = '.*___|_cell_.*', replacement = '', x = pred_id)
pred_id <- lapply(pred_id, 
                  function(x){
                    df <- readRDS(x)
                    df$barcode <- rownames(df)
                    return(df)
                  })


# filter low pred scores --------------------------------------------------

message('\nfiltering low scores \n')
metadata <- bind_rows(pred_id)
metadata$sample_id <- gsub('_.*', '', metadata$barcode)

#Grubman et.al 2019

for (i in 1:nrow(metadata)){
  x1 <- metadata$prediction.score.max[[i]]
  x2 <- sort(metadata[i,c(2:9)], decreasing = T)[[2]] # '2:7' changed to '2:9', (2 additional cell types) 16feb22
  metadata$x1[i] <- x1 ; metadata$x2[i] <- x2
  metadata$score.filter[i] <- (x1 - x2)/x1
}

###
mtx <- metadata[, grep('prediction.score.*', colnames(metadata))]
mtx <- mtx[, -grep('prediction.score.max', colnames(mtx))]
mtx <- as.matrix(mtx)
mtx <- apply(mtx, 1, sort, decreasing = T)
mtx <- t(mtx)
grubman_filter <- (mtx[, 1] - mtx[, 2]) / mtx[, 1]
metadata$grubman_filter <- grubman_filter
###

metadata <- metadata[metadata$grubman_filter > 0.20 & metadata$prediction.score.max > 0.50,] 
metadata <- metadata[!(metadata$predicted.id %in% c('Endo', 'VLMC')), ]

atac_merged <- subset(atac_merged, cells = rownames(metadata)) ; gc()
atac_merged <- AddMetaData(atac_merged, metadata = metadata)




# run svd on merged data --------------------------------------------------

message('\nrunning lsi... \n')
atac_merged <- FindTopFeatures(atac_merged, min.cutoff = 'q0')
atac_merged <- RunSVD(atac_merged, features = rownames(atac_merged), reduction.name = 'lsi.all.features')
gc()



# run harmony -------------------------------------------------------------

message('\nrunning harmony... \n')
library(harmony)
atac_merged <- RunHarmony(atac_merged, group.by.vars = 'sample_id', theta = 10, lambda = 0.05,sigma = 0.01,
                          dims.use = 2:30, reduction = 'lsi.all.features', assay.use = 'combined_peaks',
                          project.dim = FALSE, 
                          max.iter.harmony = 10, 
                          epsilon.harmony = -Inf,  
                          epsilon.cluster = 1e-6, 
                          reduction.save = 'harmony.all.features')

DimPlot(atac_merged, reduction = 'harmony.all.features', group.by = 'sample_id')


setwd(atac_integration_objects)
saveRDS(atac_merged, '4.run_harmony___integrated_atac_data.rds')

