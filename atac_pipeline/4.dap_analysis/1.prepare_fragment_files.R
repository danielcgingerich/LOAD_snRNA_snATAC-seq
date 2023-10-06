cond <- grepl('/hpc_path_identifier', getwd())
if (cond){
  source('path1/to/config.R')
} else {
  source('path2/to/config.R')
} ; rm(cond)


slurm.array.task.id <- Sys.getenv('SLURM_ARRAY_TASK_ID')
i <- as.numeric(slurm.array.task.id)

setwd(atac_preprocessing_objects)
atac <- list.files(pattern = '4.dim_red')[i]
sample_id <- gsub('.*___|_seurat_.*', '', atac)
atac <- readRDS(atac)
setwd(atac_integration_objects)
atac_integrated <- readRDS('5.cluster_links___atac_integrated_and_linked.rds')


# macs2 works on mac os or *nix systems only. so you may need to move files to another OS.  
# if moving files from windows to macOS/unix, need to change the @path slot in the fragment object.
frags <- Fragments(atac)
new_path <- paste0('/path/to/fragment_files/', 
                   sample_id, '/fragments.tsv.gz')
frags[[1]]@path <- new_path
Fragments(atac) <- NULL
Fragments(atac) <- frags

# subset atac -------------------------------------------------------------

message('subsetting data to split...')
cells_use <- colnames(atac_integrated)
cells_use <- cells_use[ grep(pattern = paste0('^', sample_id, '_'), x = cells_use)]

atac <- subset(atac, cells = cells_use)
atac_meta_subset <- atac_integrated@meta.data[colnames(atac),]
atac@meta.data <- data.frame(row.names = colnames(atac))
atac <- AddMetaData(atac, metadata = atac_meta_subset)


# split fragments ---------------------------------------------------------

message('splitting fragments...')
frag_path <- Fragments(atac)[[1]]@path
frag_path <- gsub('fragments.tsv.gz', '', frag_path)
SplitFragments(atac, group.by = 'cluster_links', outdir = frag_path, 
               file.suffix = paste0('_', sample_id))

