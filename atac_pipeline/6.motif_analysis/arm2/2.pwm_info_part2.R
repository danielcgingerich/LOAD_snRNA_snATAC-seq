cond <- grepl('/hpc_path_identifier', getwd())
if (cond){
  source('path1/to/config.R')
} else {
  source('path2/to/config.R')
} ; rm(cond)


i <- Sys.getenv('SLURM_ARRAY_TASK_ID') %>% as.numeric()

library(motifmatchr)
library(Rsamtools)
library(TFBSTools)
library(JASPAR2020)
library(TFMPvalue)

setwd(atac_motif_enrichment_objects)

files <- list.files(pattern = '2.pwm_info___.*homer.*.motifs')
files <- lapply(files, readLines)
files <- unlist(files)
write(files, '2.pwm_info___homer.formatted.all.motifs')

system('rm set1_MA*.formatted.motifs')

files <- list.files(pattern = 'detection.thresholds.rds')
n <- gsub('2.pwm_info___|_det.*', '', files)
files <- lapply(files, readRDS)
df <- unlist(files) %>% matrix(ncol = 6, byrow = T) %>% as.data.frame()
rownames(df) <- n
colnames(df) <- c('e_neg5', 'e_neg6', 'e_neg7', 'e_neg8', 'e_neg9', 'e_neg10')
saveRDS(df, '2.pwm_info___all.motif.detection.thresholds.rds')

system('rm set1_MA*thresholds.rds')
