cond <- grepl('/hpc_path_identifier', getwd())
if (cond){
  source('path1/to/config.R')
} else {
  source('path2/to/config.R')
} ; rm(cond)


library(GenomeInfoDb)
library(patchwork)
library(Matrix.utils)



# normalize + lsi + clustering --------------------------------------------
# objects processed simultaneously via slurm array
slurm.array.id <- Sys.getenv('SLURM_ARRAY_TASK_ID')                                               
i <- as.numeric(slurm.array.id)      

message(i)
atac <- list.files(pattern = 'seurat_object_qc_filtered.rds')             
sample.id <- gsub(pattern = '.*___|_seurat.*', replacement = '', atac)
atac <- readRDS(atac[i])    



# normalize ---------------------------------------------------------------
atac <- RunTFIDF(atac)   



# singular value decomp ---------------------------------------------------
atac <- FindTopFeatures(atac, min.cutoff = 'q0')                                                  
atac <- RunSVD(atac)



# nonlinear dim reduction (for visualization only) -------------------------
atac <- RunUMAP(atac, reduction = 'lsi', dims = 2:30)



# find communities --------------------------------------------------------
atac <- FindNeighbors(atac, reduction = 'lsi', dims = 2:30)
atac <- FindClusters(atac, algorithm = 3)



# save --------------------------------------------------------------------
message('saving...')
saveRDS(atac, paste0('4.dim_reduction_clustering___', sample.id[i], '_seurat_qc_dim_reduced.rds'))
