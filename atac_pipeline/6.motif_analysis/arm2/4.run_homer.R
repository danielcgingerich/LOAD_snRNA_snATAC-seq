# export PATH="/path/to/homer/bin/:$PATH"

cond <- grepl('/hpc_path_identifier', getwd())
if (cond){
  source('path1/to/config.R')
} else {
  source('path2/to/config.R')
} ; rm(cond)
i <- Sys.getenv('SLURM_ARRAY_TASK_ID') %>% as.numeric

runHomer <- function(directionality){
  # define input
  message('defining input')
  setwd(atac_motif_enrichment_objects)
  
  # target regions
  peakfile <- list.files(pattern = paste0(directionality, '_target.peaks.bed'))[i] 
  peakfile <- paste0(atac_motif_enrichment_objects, peakfile)
  
  cluster <- gsub('.*___|_uni.*|_mixed.*', '', peakfile)
  
  # bg regions
  bgfile <- paste0('3.target_background_peaks___', cluster, '_', directionality, '_background.peaks.bed') 
  bgfile <- paste0(atac_motif_enrichment_objects, bgfile)
  
  # sequences
  fa.path <- '/path/to/GRCh38.primary_assembly.genome.fa' 
  
  # pwms
  all.motifs <- paste0(atac_motif_enrichment_objects, '2.pwm_info___homer.formatted.all.motifs') # pwms to use
  
  # output directory name
  outdir <- paste0(cluster, '_', directionality)
  
  message('peak file: ', peakfile)
  message('cluster: ', cluster)
  message('background peaks: ', bgfile)
  message('fast a files: ', fa.path)
  message('input pwms: ', all.motifs)
  message('output directory: ', outdir)
  
  # define output path, build cmd, run cmd 
  
  message('creating output directory')
  setwd(atac_motif_enrichment_objects)
  setwd('homer_outs')
  system(paste0('mkdir', ' ', outdir))
  
  cmd <- paste0('findMotifsGenome.pl', ' ', # call the script
                peakfile, ' ', # target peaks
                'grch38p13', ' ', # genome I configured
                outdir,  '/ ', # output directory name
                '-size given', ' ', # region sizes
                '-mknown', ' ', all.motifs, ' ', # specify known motifs to search for
                '-bg', ' ', bgfile, ' ', # background peaks if specified
                '-mask', ' ',# should sequences be masked?
                '-nomotif') # skip de novo analysis

  message('running command:\n', cmd)
  system(cmd)
}


runHomer(directionality = 'unidirectional')
runHomer(directionality = 'mixed')

