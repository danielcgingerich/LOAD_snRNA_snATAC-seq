cond <- grepl('/hpc_path_identifier', getwd())
if (cond){
  source('path1/to/config.R')
} else {
  source('path2/to/config.R')
} ; rm(cond)

library(GenomicRanges)

# load objects ------------------------------------------------------------
setwd(atac_cicero_objects)
ccans <- list.files(pattern = '.cicero_ccans.rds')
ccans <- lapply(ccans, 
                function(x){
                  cluster <- gsub('2.run_cicero___|.cicero.*', '', x)
                  tmp <- readRDS(x)
                  tmp$cluster <- cluster
                  rownames(tmp) <- NULL
                  return(tmp)
                }) %>% bind_rows

setwd(atac_dap_analysis_objects)
daps <- list.files(pattern = '10.merge_pval.*rds$')
daps <- lapply(daps, 
               function(x){
                 cluster <- gsub('.*___|.fast.*', '', x)
                 tmp <- readRDS(x)
                 tmp$peak <- rownames(tmp)
                 rownames(tmp) <- NULL
                 tmp$cluster <- cluster
                 return(tmp)
               }) %>% bind_rows

setwd(rna_deg_analysis_objects)
degs <- list.files(pattern = '_model_results.csv')
degs <- lapply(degs, 
               function(x){
                 cluster <- gsub('.*___|_model_results.csv', '', x)
                 tmp <- read.csv(x)
                 tmp$cluster <- cluster
                 return(tmp)
               }) %>% bind_rows

# ccan overlaps >=1 dap ----------------------------------------------------------
ccan_peaks <- unique(ccans$Peak)
ccans$peak_cluster <- paste0(ccans$Peak, '_', ccans$cluster)
ccans$ccan_cluster <- paste0(ccans$CCAN, '_', ccans$cluster)
daps$peak_cluster <- paste0(daps$peak, '_', daps$cluster)

ccans <- left_join(ccans, daps[, c('avg_log2FC', 'p_val_adj', 'peak_cluster')], by = 'peak_cluster')

# ccan overlaps >=1 pro/int1 of deg ---------------------------------------------

degs <- degs[degs$cluster %in% unique(ccans$cluster), ]

# get introns

setwd(ref_genome_path)
hg38 <- rtracklayer::import('gencode.v32.primary_assembly.annotation.gtf')
deg_coords <- hg38[hg38$gene_name %in% degs$primerid]
deg_exons <- deg_coords[deg_coords$type == 'exon',]
deg_exons <- split(deg_exons, f = deg_exons$gene_name)
deg_exons <- lapply(deg_exons, reduce)
deg_exons <- deg_exons[sapply(deg_exons, length) > 1]
deg_introns <- lapply(deg_exons, 
                      function(x) {
                        gr = GRanges(seqnames = seqnames(x)[1], 
                                     ranges = IRanges(start=min(start(x)),end=max(end(x))), 
                                     strand = strand(x)[1])
                        db = disjoin(c(x, gr))
                        ints = db[countOverlaps(db, x) == 0]
                        if(as.character(strand(ints)[1]) == "-") {
                          ints$intron_id = c(length(ints):1)
                        } else {
                          ints$intron_id = c(1:length(ints))
                        }
                        ints
                      })
intron1 <- lapply(1:length(deg_introns), 
                  function(x){
                    deg_introns[[x]]$gene_name <- names(deg_introns)[x]
                    deg_introns[[x]] <- deg_introns[[x]][deg_introns[[x]]$intron_id == 1,]
                    return(deg_introns[[x]])
                  })
intron1 <- GRangesList(intron1)
intron1 <- unlist(intron1)
intron1$type <- 'intron'

deg_promoters <- promoters(deg_coords[deg_coords$type == 'gene',])
deg_promoters$type <- 'promoter'
deg_cre <- sort(c(deg_promoters, intron1), ignore.strand = F)

setwd(atac_cicero_objects)
saveRDS(deg_cre, '3.annotate_ccans___deg_promoters_and_intron1.rds')


setwd(atac_cicero_objects)
deg_cre <- readRDS('3.annotate_ccans___deg_promoters_and_intron1.rds')
deg_cre$coords <- GRangesToString(deg_cre)

ccans <- split(ccans, f = ccans$cluster)
degs <- split(degs, f = degs$cluster)

ccans <- ccans[names(degs)] # there are less deg clusters
nrow(bind_rows(ccans))
# [1] 180407 <- some were removed bc of removing clusters

for (i in 1:length(ccans)){
  cluster <- names(ccans)[i]
  
  ccans_tmp <- ccans[[cluster]]
  degs_tmp <- degs[[cluster]]
  
  message(cluster)
  genes <- degs_tmp$primerid
  deg_ranges <- deg_cre[deg_cre$gene_name %in% genes ]
  ccan_ranges <- ccans_tmp$Peak
  ccan_ranges <- StringToGRanges(ccan_ranges)
  ovrlp <- findOverlaps(query = ccan_ranges, subject = deg_ranges)
  
  deg_hits <- subjectHits(ovrlp)
  ccan_hits <- queryHits(ovrlp)
  
  ccans_tmp$deg_overlap <- ''
  ccans_tmp$deg_logfc <- NA
  ccans_tmp$deg_overlap_type <- ''
  ccans_tmp$deg_overlap_coords <- ''
  
  ccans_tmp$deg_overlap[ccan_hits] <- deg_ranges[deg_hits]$gene_name
  ccans_tmp$deg_overlap_type[ccan_hits] <- deg_ranges[deg_hits]$type
  ccans_tmp$deg_overlap_coords[ccan_hits] <- deg_ranges[deg_hits]$coords
  
  rownames(degs_tmp) <- degs_tmp$primerid
  ccans_tmp$deg_logfc[ccan_hits] <- degs_tmp[deg_ranges$gene_name[deg_hits], 'coef']
  
  ccans[[cluster]] <- ccans_tmp
}

ccans <- bind_rows(ccans)

# ccan directionality -----------------------------------------------------

dap_ccans <- ccans[ccans$p_val_adj <= 0.05, 'ccan_cluster']
deg_ccans <- ccans[ccans$deg_overlap != '', 'ccan_cluster']
dap_deg_ccans <- intersect(dap_ccans, deg_ccans)

dap_deg_ccans <- ccans[ccans$ccan_cluster %in% dap_deg_ccans, ]
dap_deg_ccans <- split(dap_deg_ccans, f = dap_deg_ccans$ccan_cluster)

dap_deg_ccans <- lapply(dap_deg_ccans, 
                        function(x){
                          dap_fc <- x[x$p_val_adj <= 0.05, 'avg_log2FC']
                          deg_fc <- x[!is.na(x$deg_logfc), 'deg_logfc']
                          direction <- expand.grid(dap_fc, deg_fc)
                          cond <- direction[,1]*direction[,2] > 0
                          
                          uni_cond <- all(cond) # if all are TRUE
                          bi_cond <- all(!cond) # if all are FALSE
                          mix_cond <- length(unique(cond)) > 1 # if mixed TRUE and FALSE
                          
                          if (uni_cond){
                            directionality <- 'unidirectional'
                          }
                          if (bi_cond){
                            directionality <- 'bidirectional'
                          }
                          if (mix_cond){
                            directionality <- 'mixed'
                          }
                          x$directionality <- directionality
                          return(x)
                        })
dap_deg_ccans <- bind_rows(dap_deg_ccans)


# gwas snp overlaps ------------------------------------------------------------

setwd(atac_motif_enrichment_objects)
kunkle_snps <- read.csv('kunkle.snps.csv', header = T)
download.file(url = 'https://hgdownload.cse.ucsc.edu/goldenpath/hg19/liftOver/hg19ToHg38.over.chain.gz',
              destfile = paste0(ref_genome_path, 'hg19ToHg38.over.chain.gz'))
# decompress file with windows 7-zip or equivalent software
hg19_to_hg38 <- rtracklayer::import.chain(paste0(ref_genome_path, 'hg19ToHg38.over.chain/hg19ToHg38.over.chain'))
kunkle_snps_hg19 <- paste0('chr', kunkle_snps$chr, '-', kunkle_snps$pos, '-', kunkle_snps$pos) %>% StringToGRanges()
kunkle_snps_hg19$rsid <- kunkle_snps$rsid
kunkle_snps_hg19$gene <- kunkle_snps$gene
kunkle_snps_hg38 <- rtracklayer::liftOver(x = kunkle_snps_hg19, chain = hg19_to_hg38) %>% unlist()

target_peaks <- unique(dap_deg_ccans$Peak)
target_peaks <- StringToGRanges(target_peaks)

annotateGWAS <- function(peaks, gwas_snps, max_distance){
  message('snp ranges must have column called gene')
  ovrlp <- distanceToNearest(x = peaks, subject = gwas_snps, 
                             ignore.strand = T)
  gwas_hits <- subjectHits(ovrlp)
  peak_hits <- queryHits(ovrlp)
  distances <- mcols(ovrlp)$distance
  df <- data.frame('peak' = GRangesToString(peaks)[peak_hits],
                   'snp' = GRangesToString(gwas_snps)[gwas_hits],
                   'snp_name' = gwas_snps$gene[gwas_hits],
                   'distance' = distances
  )
  df <- df[df$distance <= max_distance, ]
  return(df)
}

# kunkle 
kunkle_ovrlps <- annotateGWAS(peaks = target_peaks, 
                              gwas_snps = kunkle_snps_hg38,
                              max_distance = 5e5)
kunkle_ovrlps$kunkle_summary <- paste0('kunkle:', kunkle_ovrlps[,3], ',',
                                       kunkle_ovrlps[,4])
all(kunkle_ovrlps$distance <= 5e5)

# bellenguez 
setwd(atac_motif_enrichment_objects)
bellenguez <- read.csv('Bellenguez_LOAD_GWAS_APOE_MEF2C_85_loci.csv', header = T)
bellenguez_loci <- paste0('chr', bellenguez$chromosome, '-', 
                          bellenguez$list.of.loci, '-', 
                          bellenguez$list.of.loci)
bellenguez_loci <- StringToGRanges(bellenguez_loci)
bellenguez_loci$gene <- bellenguez$name
bellenguez_ovrlps <- annotateGWAS(peaks = target_peaks, 
                                  gwas_snps = bellenguez_loci,
                                  max_distance = 5e5)
bellenguez_ovrlps$bellenguez_summary <- paste0('bellenguez:', 
                                               bellenguez_ovrlps[,3], ',',
                                               bellenguez_ovrlps[,4])
all(bellenguez_ovrlps$distance <= 5e5)

# join all information
ovrlp_summary <- full_join(kunkle_ovrlps, bellenguez_ovrlps, by = 'peak')
ovrlp_summary[is.na(ovrlp_summary)] <- ""
ovrlp_summary$snp_overlap <- paste(ovrlp_summary$kunkle_summary, 
                                   ovrlp_summary$bellenguez_summary, sep = '_')
ovrlp_summary$snp_overlap <- gsub('^_|_$', '', ovrlp_summary$snp_overlap)
ovrlp_summary <- ovrlp_summary[, c('peak', 'snp_overlap')]

dap_deg_ccans <- full_join(x = dap_deg_ccans, y = ovrlp_summary, by = c('Peak'='peak'))
ix <- is.na(dap_deg_ccans$snp_overlap)
dap_deg_ccans$snp_overlap[ix] <- ""

# save work ---------------------------------------------------------------

setwd(atac_cicero_objects)
saveRDS(ccans, '3.annotate_ccans___all.ccans.rds')
saveRDS(dap_deg_ccans, '3.annotate_ccans___annotated.ccans.rds')

