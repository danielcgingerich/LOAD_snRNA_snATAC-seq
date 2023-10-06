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

setwd(atac_motif_enrichment_objects)


# get snp info ----------------------------------------------------------
pwms <- readRDS('4.detection_thresholds___empirical_pwms.rds')
snps <- readRDS('2.atsnp_target_snps_and_maf___atsnp_target_snps.rds')
snp_info <- readRDS('3.compute_motif_score___snp_info.rds')

base_table <- data.table(rsid = snp_info$snpids, snp_base = snp_info$snp_base, ref_base = snp_info$ref_base)
names(snps) <- snps$RefSNP_id
snps <- snps[base_table$rsid]

snp_position <- GRanges(snps)
snp_position <- GRangesToString(snps)
snp_position <- paste0('chr',snp_position)
base_table$snp_position <- snp_position

base_table <- cbind(base_table, as.data.table(mcols(snps)))
base_table <- base_table[, .(rsid, snp_base, ref_base, snp_position, gnomad_maf, bravo_maf, rare_or_common)]

ref_letter <- getSeq(x = BSgenome.Hsapiens.UCSC.hg38, names = StringToGRanges(snp_position), strand = '+')
base_table$ref_allele <- as.character(ref_letter)
number_to_base <- c('1' = 'A', '2' = 'C', '3' = 'G', '4' = 'T')
base_table$snp_allele <- number_to_base[ as.character(base_table$snp_base)]


# create sequence matrixes ------------------------------------------------
i <- Sys.getenv('SLURM_ARRAY_TASK_ID') %>% as.numeric()

message('array ID: ', i)
pwm_i <- pwms[[i]]

regions <- StringToGRanges(base_table$snp_position)

regions <- suppressWarnings(
  Extend(regions, ncol(pwm_i)-1, ncol(pwm_i)-1))
ref_seqs <- getSeq(x = BSgenome.Hsapiens.UCSC.hg38, names = regions, strand = '+')
ref_seqs_rev <- reverseComplement(ref_seqs)

middle <- rep(x = ncol(pwm_i), length(ref_seqs))
middle <- as.list(middle)
replace_at <- IRangesList(start = middle, end = middle)
snp_seqs <- replaceAt(x = ref_seqs, at = replace_at, 
                      value = as.list(base_table$snp_allele))
snp_seqs_rev <- reverseComplement(snp_seqs)

seqs <- list(ref_seqs = ref_seqs, ref_seqs_rev = ref_seqs_rev, 
             snp_seqs = snp_seqs, snp_seqs_rev = snp_seqs_rev)

seqs <- lapply(seqs, 
               function(x){
                 x <- as.character(x)
                 x <- strsplit(x, split = '')
                 do.call(what = rbind, args = x)
               })


# scan sequences for motifs -----------------------------------------------
pwmScan <- function(pwm, sequence_matrix){
  
  frames <- lapply(1:ncol(pwm), 
                   function(p){ 
                     sequence_matrix[, p:(p+ncol(pwm)-1) ] 
                   } 
  )
  
  scores <- lapply(X = frames, pwm = pwm,
                   function(frame_i, pwm){
                     
                     suppressWarnings(
                       ix <- melt(frame_i)
                     )
                     ix <- as.matrix(ix)
                     ix <- ix[, 3:2] # pwm index
                     ix <- gsub(' ', '', ix)
                     
                     colnames(pwm) <- as.character(1:ncol(pwm))
                     values <- pwm[ ix ]
                     values <- matrix(values, ncol = ncol(pwm))
                     rowSums(values)
                     
                   })
  scores <- do.call(what = cbind, args = scores)
  max_ix <- apply(scores, MARGIN = 1, FUN = which.max)
  max_ix <- cbind(1:length(max_ix), max_ix)
  best_score <- scores[max_ix]
  motif_start <- max_ix[, 2]
  motif_end <- motif_start + ncol(pwm) - 1
  data.table(pwm_score = best_score, motif_start = motif_start, motif_end = motif_end)
}

ref_scores <- pwmScan(pwm = pwm_i, sequence_matrix = seqs$ref_seqs)
rev_ref_scores <- pwmScan(pwm = pwm_i, sequence_matrix = seqs$ref_seqs_rev)
snp_scores <- pwmScan(pwm = pwm_i, sequence_matrix = seqs$snp_seqs)
rev_snp_scores <- pwmScan(pwm = pwm_i, sequence_matrix = seqs$snp_seqs_rev)

ref_scores$strand <- '+'
snp_scores$strand <- '+'

rev_ref_scores$strand <- '-'
rev_snp_scores$strand <- '-'




# take best match per sequence --------------------------------------------
combineFwdAndRevScores <- function(fwd, rev){
  
  score_mtx <- cbind(fwd$pwm_score, rev$pwm_score)
  score_mtx <- as.matrix(score_mtx)
  ix <- apply(score_mtx, 1, which.max)
  ix <- cbind(1:length(ix), ix)
  motifs <- data.table(id = 1:nrow(fwd))
  
  for (column_i in colnames(fwd)){
    mtx <- cbind(fwd[, ..column_i], rev[, ..column_i])
    mtx <- as.matrix(mtx)
    motifs[, column_i] <- mtx[ix]
  }
  return(motifs)
}

ref_motifs <- combineFwdAndRevScores(fwd = ref_scores, rev = rev_ref_scores)
snp_motifs <- combineFwdAndRevScores(fwd = snp_scores, rev = rev_snp_scores)

colnames(ref_motifs) <- paste0(colnames(ref_motifs), '_ref')
colnames(snp_motifs) <- paste0(colnames(snp_motifs), '_snp')

unfiltered_motifs <- bind_cols(list(base_table, ref_motifs, snp_motifs))




# convert positions to genomic coordinates --------------------------------
unfiltered_motifs$width <- unfiltered_motifs$motif_end_ref - unfiltered_motifs$motif_start_ref + 1

convertCoords <- function(pos_, width_, start_, end_, strand_){
  
  pos_ <- StringToGRanges(pos_)
  chr <- as.character(seqnames(pos_))
  
  plus_new_start <- start(pos_) - width_ + start_
  plus_new_end <- plus_new_start + width_ - 1
  
  minus_new_end <- start(pos_) + width_ - start_
  minus_new_start <- minus_new_end - width_ + 1
  
  ix <- ifelse(strand_ == '+', yes = 1, no = 2)
  ix <- cbind(1:length(ix), ix)
  
  new_start <- cbind(plus_new_start, minus_new_start)
  new_end <- cbind(plus_new_end, minus_new_end)
  
  new_start <- new_start[ix] %>% as.numeric()
  new_end <- new_end[ix] %>% as.numeric()
  
  paste0(chr, '-', new_start, '-', new_end)
}

unfiltered_motifs$tfbs_position_ref <- 
  convertCoords(pos_ = unfiltered_motifs$snp_position, 
                width_ = unfiltered_motifs$width, 
                start_ = unfiltered_motifs$motif_start_ref, 
                end_ = unfiltered_motifs$motif_end_ref,
                strand_ = unfiltered_motifs$strand_ref)

unfiltered_motifs$tfbs_position_snp <- 
  convertCoords(pos_ = unfiltered_motifs$snp_position, 
                width_ = unfiltered_motifs$width, 
                start_ = unfiltered_motifs$motif_start_snp, 
                end_ = unfiltered_motifs$motif_end_snp,
                strand_ = unfiltered_motifs$strand_snp)

# filter out motifs -------------------------------------------------------
setwd(atac_motif_enrichment_objects)
threshold <- readRDS('4.detection_thresholds___all_empirical_detection_thresholds.rds')
motif_id <- names(pwms)[i]
threshold <- threshold[motif == motif_id, `p<5e-5`]

cond1 <- unfiltered_motifs$pwm_score_ref >= threshold
cond2 <- unfiltered_motifs$pwm_score_snp >= threshold

unfiltered_motifs$motif <- motif_id

all_motifs <- unfiltered_motifs[cond1 | cond2, ]
ref_motifs <- unfiltered_motifs[cond1, ]
snp_motifs <- unfiltered_motifs[cond2, ]




# save --------------------------------------------------------------------
message('saving')
setwd(atac_motif_enrichment_objects)
saveRDS(all_motifs, paste0('7.atsnp_candidate_motifs___', motif_id, '_passed_detection_threshold.rds'))
saveRDS(ref_motifs, paste0('7.atsnp_candidate_motifs___', motif_id, '_passed_detection_threshold_major_allele.rds'))
saveRDS(snp_motifs, paste0('7.atsnp_candidate_motifs___', motif_id, '_passed_detection_threshold_minor_allele.rds'))


# 50 gb per node
# 30 task per node
# 100 at once 



# merge motifs ------------------------------------------------------------
setwd(atac_motif_enrichment_objects)
threshold <- readRDS('4.detection_thresholds___all_empirical_detection_thresholds.rds')
motif_id <- threshold$motif

all_motifs <- paste0('7.atsnp_candidate_motifs___', motif_id, '_passed_detection_threshold.rds')
ref_motifs <- paste0('7.atsnp_candidate_motifs___', motif_id, '_passed_detection_threshold_major_allele.rds')
snp_motifs <- paste0('7.atsnp_candidate_motifs___', motif_id, '_passed_detection_threshold_minor_allele.rds')

all_motifs <- lapply(all_motifs, readRDS)
ref_motifs <- lapply(ref_motifs, readRDS)
snp_motifs <- lapply(snp_motifs, readRDS)

all_motifs <- bind_rows(all_motifs)
ref_motifs <- bind_rows(ref_motifs)
snp_motifs <- bind_rows(snp_motifs)

system('rm 7.atsnp_candidate_motifs___*passed_detection_threshold*')

saveRDS(all_motifs, '7.atsnp_candidate_motifs___all_motifs_passed_detection_threshold.rds')
saveRDS(ref_motifs, '7.atsnp_candidate_motifs___ref_motifs_passed_detection_threshold.rds')
saveRDS(snp_motifs, '7.atsnp_candidate_motifs___snp_motifs_passed_detection_threshold.rds')

