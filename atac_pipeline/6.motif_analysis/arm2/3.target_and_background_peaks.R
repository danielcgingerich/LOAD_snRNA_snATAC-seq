cond <- grepl('/hpc_path_identifier', getwd())
if (cond){
  source('path1/to/config.R')
} else {
  source('path2/to/config.R')
} ; rm(cond)

library(GenomicRanges)
library(rtracklayer)


getTargetPeaks <- function(directionality){
  setwd(atac_cicero_objects)
  ccans <- readRDS('3.annotate_ccans___annotated.ccans.rds')
  ccans <- ccans[ccans$directionality == directionality, ]
  cluster <- unique(ccans$cluster)
  
  for (i in 1:length(cluster)){
    message(cluster[i])
    poi <- ccans$Peak[ccans$cluster == cluster[i]]
    poi <- StringToGRanges(poi)
    strand(poi) <- '+'
    setwd(atac_motif_enrichment_objects)
    export.bed(poi, paste0('3.target_background_peaks___', cluster[i], '_', directionality, '_target.peaks.bed'))
  }
}

getBackgroundPeaks <- function(directionality){
  setwd(atac_cicero_objects)
  ccans <- readRDS('3.annotate_ccans___annotated.ccans.rds')
  ccans <- ccans[ccans$directionality == directionality, ]
  cluster <- unique(ccans$cluster)
  
  setwd(atac_dap_analysis_objects)
  peaks <- readRDS('4.combine_consensus_peaks___cluster_specific_consensus_peaks_c=2.rds')
  
  setwd(motif_enrichment_objects)
  for (i in 1:length(cluster)){
    message(cluster[i])
    poi <- ccans$Peak[ccans$cluster == cluster[i]]
    bg <- peaks[grep(cluster[i], peaks$peak.called.in)]
    bg <- GRangesToString(bg)
    bg <- bg[ !(bg %in% poi) ]
    bg <- StringToGRanges(bg) 
    strand(bg) <- '+'
    setwd(atac_motif_enrichment_objects)
    export.bed(bg, paste0('3.target_background_peaks___', cluster[i], '_', directionality, '_background.peaks.bed'))
  }
}

# uni ccans
getTargetPeaks(directionality = 'unidirectional')
getBackgroundPeaks(directionality = 'unidirectional')
# mixed ccans
getTargetPeaks(directionality = 'mixed')
getBackgroundPeaks(directionality = 'mixed')
