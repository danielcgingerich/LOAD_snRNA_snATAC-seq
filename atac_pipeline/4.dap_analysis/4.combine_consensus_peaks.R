cond <- grepl('/hpc_path_identifier', getwd())
if (cond){
  source('path1/to/config.R')
} else {
  source('path2/to/config.R')
} ; rm(cond)

setwd(atac_dap_analysis_objects)
setwd('recalled_peaks/consensus_peaks')

path_list <- list.files()
path_list <- paste0(path_list, '/ConsensusPeaks.bed')
cluster <- gsub('_.*', '', path_list)

names(path_list) <- cluster
cpeaks <- lapply(path_list, read.table, header = T)
cpeaks <- lapply(cpeaks, makeGRangesFromDataFrame)
cpeaks <- lapply(1:length(cpeaks), 
                 function(x){
                   ranges <- cpeaks[[x]]
                   ranges$peak.called.in <- names(cpeaks)[x]
                   return(ranges)
                 })
gr.combined <- do.call("c", cpeaks)
gr <- reduce(gr.combined, with.revmap = T, min.gapwidth = 0)



ident.vec <- vector(mode = 'character', length = length(gr))
for (i in 1:length(cluster)){
  message(cluster[i])
  subject <- gr.combined[gr.combined$peak.called.in == cluster[i],]
  ovrlp <- findOverlaps(query = gr, subject = subject)
  ident.vec[queryHits(ovrlp)] <- paste(ident.vec[queryHits(ovrlp)], cluster[i], sep = ',')
}
ident.vec <- gsub('^,', '', ident.vec)
gr$peak.called.in <- ident.vec



#save data\~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
setwd(atac_dap_analysis_objects)
saveRDS(gr, '4.combine_consensus_peaks___cluster_specific_consensus_peaks_c=2.rds')







