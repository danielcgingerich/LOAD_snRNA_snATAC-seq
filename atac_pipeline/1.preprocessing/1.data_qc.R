cond <- grepl('/hpc_path_identifier', getwd())
if (cond){
  source('path1/to/config.R')
} else {
  source('path2/to/config.R')
} ; rm(cond)



# load objects ------------------------------------------------------------
setwd(atac_preprocessing_objects)
obj.list <- list.files(pattern = 'set0_atac.*v2.rds')
sample_id <- gsub('.*atac|_v2.rds', '', obj.list) # sample id is contained in file name
obj.list <- lapply(obj.list, readRDS)



# prepare annotations -----------------------------------------------------
message("loading annotations")
setwd(ref_genome_path)
hg38 <- rtracklayer::import('gencode.v32.primary_assembly.annotation.gtf')
genome(hg38) <- "hg38"
seqlevelsStyle(hg38) <- "UCSC"
hg38$gene_biotype <- hg38$gene_type



# add sample id ------------------------------------------------------------
for (i in 1:length(obj.list)){
  obj.list[[i]] <- RenameCells(obj.list[[i]], add.cell.id = sample_id[[i]])
}



# nucleosome signal -------------------------------------------------------
for (i in 1:length(obj.list)){
  message(i)
  obj.list[[i]] <- NucleosomeSignal(obj.list[[i]])
}



# tss enrichment ----------------------------------------------------------
for (i in 1:length(obj.list)){
  message(i)
  Annotation(obj.list[[i]]) <- hg38
  obj.list[[i]] <- TSSEnrichment(obj.list[[i]], fast = FALSE)
}

for (i in 1:length(obj.list)){
  message(i)
  setwd(atac_preprocessing_objects)
  saveRDS(obj.list[[i]], paste0('1.data_qc___', sample_id[i], '_srt_nuc.signal-tss.enrichment.rds'))
}

setwd(atac_preprocessing_objects)
obj.list <- list.files(pattern = 'srt_nuc.signal-tss.enrichment.rds')
obj.list <- lapply(obj.list, readRDS)

# pct reads + blacklist ---------------------------------------------------
for (i in 1:length(obj.list)){
  message(i)
  obj.list[[i]]$pct_reads_in_peaks <- obj.list[[i]]$peak_region_fragments / obj.list[[i]]$passed_filters * 100
  obj.list[[i]]$blacklist_ratio <- obj.list[[i]]$blacklist_region_fragments / obj.list[[i]]$peak_region_fragments
}



# add covariates for diff exprs -------------------------------------------
cohort <- readRDS('merged.cell.metadata.2.rds')
cohort <- cohort[, c(26:72)]
cohort <- unique(cohort)

for (i in 1:length(obj.list)){
  message(i)
  df <- obj.list[[i]]@meta.data
  df <- df[, !(colnames(df) %in% colnames(cohort))]
  df$barcode <- rownames(df)
  df$sample_id <- as.integer(sample_id[i])
  df <- left_join(x = df, cohort, by = c('sample_id' = 'sample.id'))
  rownames(df) <- df$barcode
  
  obj.list[[i]]@meta.data <- df
}



# save work ---------------------------------------------------------------
for (i in 1:length(obj.list)){
  message(i)
  setwd(atac_preprocessing_objects)
  saveRDS(obj.list[[i]], paste0('1.data_qc___', sample_id[i], '_srt_preprocessed_unfiltered.rds'))
}

frag.list <- lapply(obj.list, Fragments)
rds.list <- paste0('rds', c(25:(25+23)))
for (i in 1:length(obj.list)){
  message(i)
  saveRDS(frag.list[[i]], paste0('1.data_qc___', sample_id[i], '_fragment_object.rds'))
}
