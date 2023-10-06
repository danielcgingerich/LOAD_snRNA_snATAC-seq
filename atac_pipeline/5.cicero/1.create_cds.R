cond <- grepl('/hpc_path_identifier', getwd())
if (cond){
  source('path1/to/config.R')
} else {
  source('path2/to/config.R')
} ; rm(cond)

library(harmony)
library(cicero)
library(SingleCellExperiment)
slurm.array.id <- Sys.getenv('SLURM_ARRAY_TASK_ID')
i <- as.numeric(slurm.array.id)

setwd(atac_dap_analysis_objects)
metadata <- readRDS('7.add_covarates___metadata_for_dap_analysis.rds')
peaks <- readRDS('4.combine_consensus_peaks___cluster_specific_consensus_peaks_c=2.rds')
# need both harmony object and unintegrated object
atac <- readRDS('6.process_new_peaks___srt_consensus_processed_merged.rds')
setwd(atac_integration_objects)
atac_int <- readRDS('5.cluster_links___atac_integrated_and_linked.rds')



# -------------------------------------------------------------------------
cluster <- sort(unique(metadata$cluster_links))[i] %>% as.character

# subset necessary peaks/cells
cells_use <- rownames(metadata[metadata$cluster_links == cluster & metadata$diagnosis == 'LOAD',])
peaks_use <- peaks[grep(cluster, peaks$peak.called.in),] %>% GRangesToString
message('total cells in ', cluster, ': ', length(cells_use))
message('total peaks in ', cluster, ': ', length(peaks_use))

counts <- atac[['recalled.peaks']]@counts
data <- atac[['recalled.peaks']]@data

counts <- counts[peaks_use, cells_use]
data <- data[peaks_use, cells_use]

atac_subset <- CreateChromatinAssay(counts = counts, sep = c('-','-'), min.cells = -1, min.features = -1)
atac_subset <- CreateSeuratObject(atac_subset, assay = 'peaks')
atac_subset[['peaks']]@data <- data
atac_subset <- AddMetaData(atac_subset, metadata = metadata[colnames(atac_subset), ])



# add harmony embeddings --------------------------------------------------
emb <- subset(atac_int, cells = colnames(atac_subset))
atac_subset[['harmony']] <- CreateDimReducObject(embeddings = emb[['harmony.all.features']]@cell.embeddings, assay = 'peaks')

cds <- as.CellDataSet(atac_subset, reduction = 'harmony')
# we used bin size of 50 because we had enough data, 
# but if your data is smaller, consider smaller bins
bin_size <- min(50, round(ncol(cds)/4)) 
cds <- make_cicero_cds(cds, reduced_coordinates = atac_subset[['harmony']]@cell.embeddings, k = bin_size)
warnings()

setwd(atac_cicero_objects)
saveRDS(cds, paste0('1.create_cds___', cluster, '_cicero_cds.rds'))
