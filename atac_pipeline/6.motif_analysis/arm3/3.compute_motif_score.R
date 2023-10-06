cond <- grepl('/hpc_path_identifier', getwd())
if (cond){
  source('path1/to/config.R')
} else {
  source('path2/to/config.R')
} ; rm(cond)

library(SNPlocs.Hsapiens.dbSNP151.GRCh38)
library(BSgenome.Hsapiens.UCSC.hg38)
library(data.table)
library(TFBSTools)
library(JASPAR2020)
library(atSNP)

# lets use cutoff of 5%. can easily get cutoff recomendations on google.
setwd(atac_motif_enrichment_objects)
snps <- readRDS('2.atsnp_target_snps_and_maf___atsnp_target_snps.rds')
snps <- mcols(snps) %>% as.data.table()
rsid <- snps[gnomad_maf >= 0.01 | bravo_maf >= 0.01, RefSNP_id]

snp_info <- LoadSNPData(snpids = rsid,
                        genome.lib ="BSgenome.Hsapiens.UCSC.hg38",
                        snp.lib = "SNPlocs.Hsapiens.dbSNP151.GRCh38",
                        half.window.size = 25, default.par = TRUE, mutation = FALSE)
setwd(atac_motif_enrichment_objects)
saveRDS(snp_info, '3.compute_motif_score___snp_info.rds')

snp_info <- readRDS('3.compute_motif_score___snp_info.rds')


message('getting ppms')
pfm <- getMatrixSet(JASPAR2020, opts = list(species = 9606, all_versions = F))

ppm <- toPWM(pfm, type = 'prob', pseudocounts = 1)

ppm <- lapply(ppm, 
              function(x){
                mtx <- t(x@profileMatrix)
                mtx <- mtx[, c('A', 'C', 'G', 'T')]
                rownames(mtx) <- NULL
                colnames(mtx) <- NULL
                mtx
              })

setwd(atac_motif_enrichment_objects)
saveRDS(ppm, '3.compute_motif_score___atsnp_ppms.rds')

# run atSNP 
atsnp_scores <- ComputeMotifScore(motif.lib = ppm, snp.info = snp_info, ncores = 1)



setwd(atac_motif_enrichment_objects)
saveRDS(atsnp_scores, '3.compute_motif_score___atsnp_motif_scores.rds')
