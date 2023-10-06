cond <- grepl('/hpc_path_identifier', getwd())
if (cond){
  source('path1/to/config.R')
} else {
  source('path2/to/config.R')
} ; rm(cond)

library(GenomicRanges)
library(GenomicScores)
library(SNPlocs.Hsapiens.dbSNP151.GRCh38)
library(MafDb.gnomAD.r2.1.GRCh38)
library(MafDb.TOPMed.freeze5.hg38)
library(BSgenome.Hsapiens.UCSC.hg38)
library(data.table)
# lets use cutoff of 5%. can easily get cutoff recomendations on google.


# -------------------------------------------------------------------------
message('loading motif info')
setwd(atac_motif_enrichment_objects)

uni_regions <- readRDS('1.atsnp_target_dap_deg_pairs___unidirectional.rds') %>% bind_rows()
mix_regions <- readRDS('1.atsnp_target_dap_deg_pairs___mixed.rds') %>% bind_rows()

daps <- unique(c(uni_regions$dap_peak, mix_regions$dap_peak))
daps <- StringToGRanges(daps)

chr <- as.character(seqnames(daps)) %>% gsub(pattern = 'chr', replacement = '')
start <- start(daps)
end <- end(daps)
daps <- GRanges(seqnames = chr, ranges = IRanges(start = start, end = end))

snps <- snpsByOverlaps(x = SNPlocs.Hsapiens.dbSNP151.GRCh38, ranges = daps)
tmp <- mcols(snps)
chr <- as.character(seqnames(snps))
chr <- paste0('chr', chr)
snps <- GPos(seqnames = chr, pos = pos(snps))
mcols(snps) <- tmp ; rm(tmp)

###
alleles <- inferRefAndAltAlleles(gpos = snps, genome = BSgenome.Hsapiens.UCSC.hg38)
alleles <- as.data.table(alleles)

mcols(snps) <- cbind(mcols(snps), alleles)
snps$multiallelic <- sapply(snps$alt_alleles, length)
snps$multiallelic <- ifelse(snps$multiallelic > 1, 1, 0)
snps_multi <- snps[snps$multiallelic == 1]
snps_di <- snps[snps$multiallelic == 0]
snps_multi <- split(snps_multi, f = snps_multi$RefSNP_id)
snps_multi <- lapply(snps_multi, 
                     function(x){
                       alts <- unlist(x$alt_alleles)
                       gr <- rep(x, length(alts))
                       gr$alt_alleles <- alts
                       gr
                     })
snps_multi <- GRangesList(snps_multi)
snps_multi <- unlist(snps_multi)

snps <- c(snps_di, snps_multi)
snps$alt_alleles <- unlist(snps$alt_alleles)
setwd(atac_motif_enrichment_objects)
saveRDS(snps, '2.atsnp_target_snps_and_maf___snp_alleles.rds')

# reports global frequency (all ancestry combined)
snps <- snpsByOverlaps(x = SNPlocs.Hsapiens.dbSNP151.GRCh38, ranges = daps)

gnomad <- gscores(x = MafDb.gnomAD.r2.1.GRCh38, ranges = snps, ref = snps$ref_allele, alt = snps$alt_alleles)

snps_bravo <- paste0('chr', as.character(seqnames(snps)), '-', start(snps), '-', start(snps))
snps_bravo <- StringToGRanges(snps_bravo)
mcols(snps_bravo) <- mcols(snps)
bravo <- gscores(x = MafDb.TOPMed.freeze5.hg38, ranges = snps_bravo)

snps$gnomad_maf <- gnomad$AF
snps$bravo_maf <- bravo$AF
snps$gnomad_maf[is.na(snps$gnomad_maf)] <- 0
snps$bravo_maf[is.na(snps$bravo_maf)] <- 0

snps <- snps[!(snps$gnomad_maf == 0 & snps$bravo_maf == 0)]

snps$rare_or_common <- ifelse( snps$gnomad_maf > 0.001 | snps$bravo_maf > 0.001, 'common', 'rare')
table(snps$rare_or_common)


setwd(atac_motif_enrichment_objects)
saveRDS(snps, '2.atsnp_target_snps_and_maf___atsnp_target_snps.rds')


