cond <- grepl('/hpc_path_identifier', getwd())
if (cond){
  source('path1/to/config.R')
} else {
  source('path2/to/config.R')
} ; rm(cond)
library(data.table)
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg38)
library(Biostrings)
library(SNPlocs.Hsapiens.dbSNP151.GRCh38)
library(motifmatchr)
library(JASPAR2020)
library(TFBSTools)
library(TFMPvalue)



# calculate nucleotide distribution ---------------------------------------
setwd(atac_motif_enrichment_objects)
target_peaks <- list.files(pattern = 'atsnp_target_dap_deg_pairs.*rds$')
target_peaks <- lapply(target_peaks, function(x){
  x <- readRDS(x) ; x <- bind_rows(x) ; x$dap_peak
})
target_peaks <- unique( unlist(target_peaks) )
target_peaks <- StringToGRanges(target_peaks)
target_peaks <- reduce(target_peaks)

nuc_dist <- getSeq(x = BSgenome.Hsapiens.UCSC.hg38, names = target_peaks, strand = '+')
nuc_dist <- c(nuc_dist, reverseComplement(nuc_dist))
nuc_dist <- as.character(nuc_dist)
nuc_dist <- strsplit(nuc_dist, split = '')
nuc_dist <- unlist(nuc_dist)
nuc_dist <- table(nuc_dist)
n <- sum(nuc_dist)
nuc_dist <- nuc_dist / n
acgt <- names(nuc_dist)
nuc_dist <- as.numeric(nuc_dist)
names(nuc_dist) <- acgt

# make empirical pwms -----------------------------------------------------
pfm <- getMatrixSet(JASPAR2020, opts = list(species = 9606, all_versions = F))
pfm <- lapply(pfm, 
              function(x){
                x@profileMatrix
              })

ppm <- lapply(pfm, function(x){
  ppm <- x + 0.25
  ppm <- ppm %*% diag(1/(colSums(ppm)))
  return(ppm)
})

# have to use log2 odds to be compatible with tfmpvalue
log2pwm <- lapply(ppm, 
                  function(x){
                    x['A', ] <- x['A', ] / nuc_dist['A']
                    x['C', ] <- x['C', ] / nuc_dist['C']
                    x['G', ] <- x['G', ] / nuc_dist['G']
                    x['T', ] <- x['T', ] / nuc_dist['T']
                    log2(x)
                  })



# find empirical detection threshold --------------------------------------

i <- Sys.getenv('SLURM_ARRAY_TASK_ID') %>% as.numeric()
log2pwm_i <- log2pwm[[i]]
pfm_i <- pfm[[i]]

t1 <- TFMpv2sc(mat = log2pwm_i, pvalue = 5e-5, type = 'PWM', bg = nuc_dist) ###
t2 <- TFMpv2sc(mat = log2pwm_i, pvalue = 5e-6, type = 'PWM', bg = nuc_dist) ###
t3 <- TFMpv2sc(mat = log2pwm_i, pvalue = 5e-7, type = 'PWM', bg = nuc_dist) ###
t4 <- TFMpv2sc(mat = log2pwm_i, pvalue = 5e-8, type = 'PWM', bg = nuc_dist) ###
det_thresh <- c(t1, t2, t3, t4)

setwd(atac_motif_enrichment_objects)
if (i == 1){
  saveRDS(log2pwm, '4.detection_thresholds___empirical_pwms.rds')
}

saveRDS(det_thresh, paste0('4.detection_thresholds___', names(log2pwm)[i] ,'_empirical_detection_threshold.rds'))



# merge output ------------------------------------------------------------
setwd(atac_motif_enrichment_objects)
det_thresh <- list.files(pattern = 'empirical_detection_threshold.rds$')
motifs <- gsub('.*___|_empirical.*', '', det_thresh)

det_thresh <- lapply(det_thresh, readRDS)
det_thresh <- do.call(what = rbind, args = det_thresh)
det_thresh <- as.data.table(det_thresh)
det_thresh$motif <- motifs
det_thresh <- det_thresh[, .(motif = motif, `p<5e-5` = V1, 
                             `p<5e-6` = V2, `p<5e-7` = V3, 
                             `p<5e-8` = V4)]


system('rm *_empirical_detection_threshold.rds')

saveRDS(det_thresh, '4.detection_thresholds___all_empirical_detection_thresholds.rds')
