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

pfm0 <- getMatrixSet(JASPAR2020, opts = list(species = 9606, all_versions = F))
pfm <- lapply(pfm0, 
              function(x){
                x@profileMatrix
              })

ppm <- lapply(pfm, function(x){
  ppm <- x + 0.25
  ppm <- ppm %*% diag(1/(colSums(ppm)))
  return(ppm)
})

pwm <- lapply(ppm, function(x){
  odds_ratio <- x/0.25
  return(log(odds_ratio))
})

pwm <- lapply(1:length(pwm), function(x){
  mtx_info <- pfm0[[x]]
  mtx <- pwm[[x]]
  pwmtx <- PWMatrix(ID = mtx_info@ID, name = mtx_info@name, 
                    matrixClass = mtx_info@matrixClass, 
                    strand = mtx_info@strand, bg = mtx_info@bg, 
                    tags = mtx_info@tags, profileMatrix = mtx,
                    pseudocounts = 1)
  return(pwmtx)
})
names(pwm) <- lapply(pwm, function(x){x@ID})

pwm <- do.call(c(pwm, list(use.names = TRUE)), what = PWMatrixList)

if (i == 1){
  message('creating pwms for motifmatchr')
  setwd(atac_motif_enrichment_objects)
  saveRDS(pwm, '2.pwm_info___signac.formatted.pwms.rds')
}


# get motif info ----------------------------------------------------------

motif_info <- lapply(pwm, 
                     function(x){
                       mtx_id <- x@ID
                       gene <- x@name
                       gene <- gsub('\\(.*\\)', '', gene)
                       gene <- strsplit(gene, split = '::') %>% unlist
                       data.frame('id' = mtx_id, 'gene' = gene)
                     }) %>% bind_rows
ix <- grep('EWSR1', motif_info$gene)
motif_info[ix, 'gene'] <- 'EWSR1'
tmp <- data.frame('id' = 'MA0149.1', 'gene' = 'FLI1')
motif_info <- bind_rows(motif_info, tmp)

motif_names <- lapply(pwm, 
                      function(x){
                        data.frame('id' = x@ID, 'tf_name' = x@name)
                      }) %>% bind_rows
motif_info <- left_join(x = motif_info, y = motif_names, by = 'id')
setwd(motif_enrichment_objects)
if (i == 1){
  message('creating pwms for motifmatchr')
  setwd(atac_motif_enrichment_objects)
  saveRDS(motif_info, '2.pwm_info___motif.info.rds')
}


# odds ratio threshold ---------------------------------------------------------------
message('calculating pwm detection threshold: ', names(pwm)[i])


message('calculating thresholds')
mtx <- pwm[[i]]@profileMatrix

p_cutoff <- c(5e-5, 5e-6, 5e-7, 5e-8, 5e-9, 5e-10)
threshold <- numeric()
for (j in seq_along(p_cutoff)){
  threshold[j] <- TFMpv2sc(mat = mtx, pvalue = p_cutoff[j], type = 'PWM')
}
setwd(atac_motif_enrichment_objects)
filename <- paste0('2.pwm_info___', names(pwm)[i], '_detection.thresholds.rds')
saveRDS(threshold, filename)

# homer formatted pwms ----------------------------------------------------
message('creating .motif file')

motif <- names(pfm0)[i]
mtx <- pfm0[[i]]@profileMatrix
mtx <- mtx %*% diag(1/colSums(mtx)) # we need ppms
consensus_sequence <- apply(mtx, 2, 
                            function(x){
                              ix <- which.max(x)
                              rownames(mtx)[ix]
                            }) %>% paste(collapse = '')
first_line <- paste0('>', 
                     consensus_sequence, '\t',
                     motif, '\t', 
                     threshold[1])
body_lines <- sapply(1:ncol(mtx), 
                     function(x){
                       tmp <- mtx[, x]
                       tmp <- paste(tmp, collapse = '\t')
                       return(tmp)
                     })
file <- c(first_line, body_lines)

write(file, paste0('2.pwm_info___', motif, '_homer.formatted.motifs'))

message('done!')
