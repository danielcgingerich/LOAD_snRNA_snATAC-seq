cond <- grepl('/hpc_path_identifier', getwd())
if (cond){
  source('path1/to/config.R')
} else {
  source('path2/to/config.R')
} ; rm(cond)
library(glmGamPoi)
library(sctransform)

# integration -------------------------------------------------------------


setwd(rna_preprocessing_objects)
rna <- list.files(pattern = '1.data_qc.*.rds')
sample_id <- gsub('.*___rna|_qc_.*', '', rna)
rna <- lapply(rna, readRDS)
rna <- lapply(rna, 
              function(x){
                old_names <- colnames(x)
                prefix <- gsub('rna', '', x$sampID)
                new_names <- paste0(prefix, '_', old_names)
                RenameCells(x, new.names = new_names)
              })


int_features <- SelectIntegrationFeatures(object.list = rna, nfeatures = 10000)
rna <- PrepSCTIntegration(object.list = rna, anchor.features = int_features)
anchors <- FindIntegrationAnchors(object.list = rna, normalization.method = 'SCT', anchor.features = int_features)
rna_int <- IntegrateData(anchorset = anchors, normalization.method = 'SCT')

setwd(rna_integration_objects)
saveRDS(rna_int, '1.integrate_data___integrated_rna.rds')
