cond <- grepl('/hpc_path_identifier', getwd())
if (cond){
  source('path1/to/config.R')
} else {
  source('path2/to/config.R')
} ; rm(cond)

slurm.array.id <- Sys.getenv('SLURM_ARRAY_TASK_ID')
i <- as.numeric(slurm.array.id)

# set paths ---------------------------------------------------------------


setwd(atac_preprocessing_objects)
rds_list <- list.files(pattern = '4.dim_red')
sample_id <- gsub('.*___|_seurat_qc.*', '', rds_list)
setwd(atac_dap_analysis_objects)
peaks <- list.files(pattern = 'consensus_peaks')


message(paste0('beginning script.  dataset = ', sample_id[i]))


# load objects ------------------------------------------------------------

message('\nloading objects... \n')
setwd(atac_preprocessing_objects)
atac <- readRDS(rds_list[i])
setwd(atac_dap_analysis_objects)
peaks <- readRDS(peaks)


# build mtx ---------------------------------------------------------------

message('\nbuilding feature matrix... \n')
message(sample_id[i])
mtx <- FeatureMatrix(fragments = Fragments(atac), features = peaks,
                     cells = colnames(atac))


# create assay ------------------------------------------------------------

message('\ncreating assay... \n')
atac_recall <- CreateChromatinAssay(counts = mtx,
                                    sep = c(":", "-"), 
                                    min.cells = -1, min.features = -1)
message('\ncreating seurat obj... \n')
atac_recall <- CreateSeuratObject(atac_recall, assay = 'recalled.peaks')


message(paste0('saving object.  dataset = ', sample_id[i]))

setwd(atac_dap_analysis_objects)
saveRDS(atac_recall, paste0('5.quantify_peaks___', sample_id[i], '_consensus_peak_object.rds'))

message('\ndone!')


