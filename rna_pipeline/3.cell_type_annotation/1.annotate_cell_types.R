cond <- grepl('/hpc_path_identifier', getwd())
if (cond){
  source('path1/to/config.R')
} else {
  source('path2/to/config.R')
} ; rm(cond)

setwd(rna_cell_type_annotation_objects)
reference_dataset <- readRDS('m1_reference_dataset.rds')

setwd(rna_integration_objects)
combined.integrated <- readRDS('1.integrate_data___integrated_rna.rds')

DefaultAssay(reference_dataset) <- "integrated"
DefaultAssay(combined.integrated) <- "integrated"

transfer.anchors <- FindTransferAnchors(reference = reference_dataset, query = combined.integrated,
                                        normalization.method = "SCT",
                                        reference.assay = "integrated", query.assay = "integrated", 
                                        project.query = FALSE, features = intersect(rownames(reference_dataset), rownames(combined.integrated)))

predictions <- TransferData(anchorset = transfer.anchors, refdata = reference_dataset$cell_type, dims = 1:30)
combined.integrated <- AddMetaData(object = combined.integrated, metadata = predictions)

##Filter based on prediction score - remove cells with max score of less than 0.5
combined.integrated <- subset(combined.integrated, subset = prediction.score.max > 0.5)
combined.integrated

setwd(rna_cell_type_annotation_objects)
saveRDS(combined.integrated, '1.annotate_cell_types___integrated_rna_data_annotated_cell_types.rds')
