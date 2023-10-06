cond <- grepl('/hpc_path_identifier', getwd())
if (cond){
  source('path1/to/config.R')
} else {
  source('path2/to/config.R')
} ; rm(cond)

setwd(rna_cell_type_annotation_objects)
combined.integrated <- readRDS('1.annotate_cell_types___integrated_rna_data_annotated_cell_types.rds')


# find difference between 1st & 2nd pred.scores, scaled by max prediction score per each cell.
# values close to 1 are accurate

# these cells represent multiplets or poor quality cells.

mtx <- combined.integrated@meta.data[, grep(pattern = 'prediction.score', x = colnames(combined.integrated@meta.data))]
mtx <- mtx[, !(colnames(mtx) %in% 'prediction.score.max')]
mtx <- as.matrix(mtx)
sorted_scores <- apply(mtx, 1, sort, decreasing = T)
sorted_scores <- t(sorted_scores)
second_highest_score <- sorted_scores[, 2]

grubman_filter <- (combined.integrated$prediction.score.max - second_highest_score) / combined.integrated$prediction.score.max

combined.integrated$grubman_filter <- grubman_filter


combined.integrated <- subset(combined.integrated, subset = grubman_filter > 0.65)
combined.integrated

setwd(rna_cell_type_annotation_objects)
saveRDS(combined.integrated, '2.remove_hybrid_cells___integrated_rna_data_annotated_cell_types.rds')

