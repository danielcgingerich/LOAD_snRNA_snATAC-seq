cond <- grepl('/hpc_path_identifier', getwd())
if (cond){
  source('path1/to/config.R')
} else {
  source('path2/to/config.R')
} ; rm(cond)





# laod data ---------------------------------------------------------------

# rna was already integrated, but that is sct normalized data. 
# is recommended to use log normalized for cca with gene activity matrix
# so, must re-integrate rna data. 

#load srt obj's:
srt.obj <- list()
for (i in 1:8){
  message(i)
  setwd(rna_preprocessing_objects)
  files <- list.files(pattern = '1.data_qc___')
  
  srt.obj[[i]] <- readRDS(files[i])
  
  # add sample id prefix to cell names before integration
  sample_id <- srt.obj[[i]]$sampID[1]
  sample_id <- gsub('^rna', '', sample_id)
  srt.obj[[i]] <- RenameCells(srt.obj[[i]], 
                              new.names = paste0(sample_id, '_', colnames(srt.obj[[i]]) ))
  
  DefaultAssay(srt.obj[[i]]) <- 'RNA'
  srt.obj[[i]][['SCT']] <- NULL
  all.genes <- rownames(srt.obj[[i]])
  mt <- grep('^MT-', all.genes)
  unlocalized.contigs <- paste(c("BX004987.1", "AC145212.1",
                           "MAFIP",      "AC011043.1", "AC011043.2", "AC011841.1", "BX072566.1",
                           "AL354822.1", "AL592183.1", "AC240274.1", "AC213203.2", "AC213203.1",
                           "AC004556.3", "AC233755.2", "AC233755.1", "AC136352.3", "AC136352.2",
                           "AC171558.3", "AC171558.1", "AC133551.1", "AC136612.1", "AC136616.1",
                           "AC136616.3", "AC136616.2", "AC141272.1", "AC023491.2", "AC007325.1",
                           "AC007325.4", "AC007325.2"), collapse = '|')
  unlocalized.contigs <- grep(unlocalized.contigs, all.genes)
  genes_remove <- c(mt, unlocalized.contigs)
  genes.use <- all.genes[-genes_remove]
  srt.obj[[i]] <- subset(srt.obj[[i]], features = genes.use)
  srt.obj[[i]] <- NormalizeData(srt.obj[[i]])
  srt.obj[[i]] <- FindVariableFeatures(srt.obj[[i]], selection.method = 'vst', nfeatures = 5000)
  srt.obj[[i]] <- ScaleData(srt.obj[[i]])
  srt.obj[[i]] <- RunPCA(srt.obj[[i]])
  rm(genes.use)
}

message("\n selecting integration features \n")
int.features <- SelectIntegrationFeatures(object.list = srt.obj, nfeatures = 500)


reference_dataset <- which.max(sapply(srt.obj, ncol)) # choose dataset with most cells as reference. 
ref.anchors <- FindIntegrationAnchors(object.list = srt.obj, normalization.method = "LogNormalize", 
                                      anchor.features = int.features, dims = 1:30, reference = reference_dataset)
rna <- IntegrateData(anchorset = ref.anchors, normalization.method = "LogNormalize", dims = 1:30)





# add metadata ------------------------------------------------------------
setwd(rna_cell_type_annotation_objects)
metadata <- readRDS('3.pca_lovain_cluster_umap___annotated_integrated_processed_dataset.rds')
metadata <- metadata@meta.data
ix <- rownames(metadata) %in% colnames(rna)
metadata <- metadata[ix, ]

rna <- subset(rna, cells = rownames(metadata))
rna <- AddMetaData(rna, metadata = metadata)

rna$cell_type <- as.character(rna$cell_type)
rna$cell_type_cluster <- as.character(rna$cell_type_cluster)
na.index <- is.na(rna$cell_type_cluster)
rna$cell_type[na.index] <- rna$predicted.id[na.index]
rna$cell_type_cluster[na.index] <- rna$predicted.id[na.index]

setwd(atac_cell_type_annotation_objects)
saveRDS(rna, '2.snrna_lognormalized_integration___data.rds')
