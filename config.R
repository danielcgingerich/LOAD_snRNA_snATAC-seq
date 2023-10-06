# we worked on both HPC cluster and our local machines.  config file detects whether we on HPC or 
# local machine, and sets paths accordingly. 

# if hpc:
if (grepl('/hpc_path_identifier', getwd())){ 
  # on university hpc, must specify personal directory for installing and loading r packages.
  .libPaths(c('path/to/your/r/package/directory/on/hpc/cluster', .libPaths()))
  ref_genome_path <- '/path/to/reference/genome/gtf/file'
  r_fxns <- '/path/to/scripts/containing/r/functions'
  
  atac_preprocessing_objects <- '/path/to/atac_pipeline/1.preprocessing/'
  atac_cell_type_annotation_objects <- '/path/to/atac_pipeline/2.cell_type_annotation/'
  atac_integration_objects <- '/path/to/atac_pipeline/3.integration/'
  atac_dap_analysis_objects <- '/path/to/atac_pipeline/4.dap_analysis/'
  atac_cicero_objects <- '/path/to/atac_pipeline/5.cicero/'
  atac_motif_enrichment_objects <- '/path/to/atac_pipeline/6.motif_enrichment/'
  atac_cluster_link_validation_objects <- '/path/to/atac_pipeline/7.cluster_link_validation/'
  
  rna_preprocessing_objects <- '/path/to/rna_pipeline/1.preprocessing/'
  rna_cell_type_annotation_objects <- '/path/to/rna_pipeline/3.cell_type_annotation/'
  rna_integration_objects <- '/path/to/rna_pipeline/2.integration/'
  rna_deg_analysis_objects <- '/path/to/rna_pipeline/4.deg_analysis/'
  
  
  # if local machine:
} else { 
  
  ref_genome_path <- '/path/to/reference/genome/gtf/file'
  r_fxns <- '/path/to/scripts/containing/r/functions'
  
  atac_preprocessing_objects <- '/path/to/atac_pipeline/1.preprocessing/'
  atac_cell_type_annotation_objects <- '/path/to/atac_pipeline/2.cell_type_annotation/'
  atac_integration_objects <- '/path/to/atac_pipeline/3.integration/'
  atac_dap_analysis_objects <- '/path/to/atac_pipeline/4.dap_analysis/'
  atac_cicero_objects <- '/path/to/atac_pipeline/5.cicero/'
  atac_motif_enrichment_objects <- '/path/to/atac_pipeline/6.motif_enrichment/'
  atac_cluster_link_validation_objects <- '/path/to/atac_pipeline/7.cluster_link_validation/'
  
  rna_preprocessing_objects <- '/path/to/rna_pipeline/1.preprocessing/'
  rna_cell_type_annotation_objects <- '/path/to/rna_pipeline/3.cell_type_annotation/'
  rna_integration_objects <- '/path/to/rna_pipeline/2.integration/'
  rna_deg_analysis_objects <- '/path/to/rna_pipeline/4.deg_analysis/'
  
  }



# load packages
suppressWarnings(suppressMessages( library(Seurat) ))
suppressWarnings(suppressMessages( library(Signac) ))
suppressWarnings(suppressMessages( library(Signac) ))
suppressWarnings(suppressMessages( library(IRanges) ))
suppressWarnings(suppressMessages( library(GenomicRanges) ))
suppressWarnings(suppressMessages( library(dplyr) ))
suppressWarnings(suppressMessages( library(Matrix.utils) ))
suppressWarnings(suppressMessages( library(stringr) ))
suppressWarnings(suppressMessages( library(data.table) ))

message('\n[packages loaded]\n')



# test your configuration works
# loop sets working directory for all paths. if error, correct the path. 
# paths <- ls()
# paths <- paths[grepl('_objects$|_scripts', paths)]
# for (ix in 1:length(paths)){
#   d <- eval(parse(text = paths[ix]))
#   message(d)
#   setwd(d)
# }
