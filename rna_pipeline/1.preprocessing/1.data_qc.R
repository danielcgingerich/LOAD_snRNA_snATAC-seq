cond <- grepl('/hpc_path_identifier', getwd())
if (cond){
  source('path1/to/config.R')
} else {
  source('path2/to/config.R')
} ; rm(cond)
library(glmGamPoi)
library(sctransform)

setwd(rna_preprocessing_objects)
rna <- list.files(pattern = 'set0_')
sample_id <- gsub('set0_rna|.rds', '', rna)
rna <- lapply(rna, readRDS)

# qc ----------------------------------------------------------------------

# nFeatures between 200-10000, % mito <= 17.5%
rna <- lapply(rna, 
              function(x){
                x[['percent.mt']] <- PercentageFeatureSet(x, pattern = '^MT')
                return(x)
              })
meta <- lapply(rna, function(x){x@meta.data})
meta <- bind_rows(meta)
meta$bcd <- rownames(meta)
meta <- as.data.table(meta)
meta <- meta[nFeature_RNA >= 200 & nFeature_RNA <= 10000 & percent.mt <= 17.5, ]

rna <- lapply(rna, 
              function(x){
                ix <- meta$bcd %in% colnames(x)
                cells_keep <- meta$bcd[ix]
                subset(x, cells = cells_keep)
              })


# preprocessing -----------------------------------------------------------
#remove MT genes
for (i in 1:length(list)) {
  data_list <- GetAssayData(rna[[i]], assay = "RNA")
  data_list <- data_list[-(which(rownames(data_list) %in% c('MT-ND1','MT-ND2','MT-CO1','MT-CO2','MT-ATP8', 'MT-ATP6','MT-CO3', 'MT-ND3', 'MT-ND4L', 'MT-ND4', 'MT-ND5', 'MT-ND6', 'MT-CYB'))),]
  rna[[i]] <- subset(list[[i]], features = rownames(data_list))
}

for (i in 1:length(rna)) {
  rna[[i]] <- SCTransform(rna[[i]], method = 'glmGamPoi', return.only.var.genes = FALSE, verbose = FALSE)
}

setwd(rna_preprocessing_objects)
for (i in 1:length(rna)){
  message('saving')
  saveRDS(object = rna[[i]], file = paste0('1.data_qc___rna', sample_id[i], '_qc_normalized_dimreduct.rds'))
}
        



