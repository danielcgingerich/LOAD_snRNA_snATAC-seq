cond <- grepl('/hpc_path_identifier', getwd())
if (cond){
  source('path1/to/config.R')
} else {
  source('path2/to/config.R')
} ; rm(cond)

library(SingleCellExperiment)
library(cicero)
i <- Sys.getenv('SLURM_ARRAY_TASK_ID') %>% as.numeric

message('loading data...')
setwd(atac_cicero_objects)
cds <- list.files(pattern = '1.create_cds___')[i]
cluster <- gsub('.*___|.cicero.*', '', cds)
cds <- readRDS(cds)
hg38 <- read.table('http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes')
hg38 <- hg38[1:24, ]

message('running cicero...')
conns <- run_cicero(cds,
                    genomic_coords = hg38, 
                    window = 5e5)

message('saving conns...')
setwd(atac_cicero_objects)
saveRDS(conns, paste0('2.run_cicero___', cluster, '_cicero_conns.rds'))

ccans <- generate_ccans(conns, coaccess_cutoff_override = 0.2)

message('saving ccans...')
saveRDS(ccans, paste0('2.run_cicero___', cluster, '.cicero_ccans.rds'))
