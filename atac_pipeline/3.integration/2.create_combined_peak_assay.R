cond <- grepl('/hpc_path_identifier', getwd())
if (cond){
  source('path1/to/config.R')
} else {
  source('path2/to/config.R')
} ; rm(cond)

slurm.array.id <- Sys.getenv('SLURM_ARRAY_TASK_ID')
i <- as.numeric(slurm.array.id)

# load data ---------------------------------------------------------------

setwd(atac_integration_objects)
combined_peaks <- readRDS('1.create_common_peak_set___preak_granges.rds')

setwd(atac_preprocessing_objects)
frags <- list.files(pattern = 'fragment_object.rds')[i]
sample_id <- gsub('.*___|_frag.*', '', frags)
frags <- readRDS(frags)

setwd(atac_cell_type_annotation_objects)
cells_use <- readRDS(paste0('3.annotate_atac_cell_types___', sample_id, '_cell_type_predictions.rds'))
cells_use <- rownames(cells_use)

###
frags[[1]]@path <- paste0('\\path\\to\\fragment_files\\', sample_id, '\\fragments.tsv.gz')
###
# quantify fragments ------------------------------------------------------
mtx <- FeatureMatrix(fragments = frags,
                     features = combined_peaks,
                     cells = cells_use)

assay.obj <- CreateChromatinAssay(counts = mtx, sep = c(":", "-"), 
                                  fragments = frags, min.cells = -1)

srt.obj <- CreateSeuratObject(assay.obj, assay = 'combined_peaks')

setwd(atac_integration_objects)
saveRDS(srt.obj, paste0('2.create_combined_peak_assay___', sample_id, '.rds'))
