# Integrative analysis of parallel single nucleus ATAC and RNA seq data

This repository contains the analysis pipeline used in Gamache, Gingerich & Shwab et al. “Integrative single-nucleus multi-omics analysis prioritizes candidate cis and trans regulatory networks and their target genes in Alzheimer’s disease brains”.  The pipeline was created for examining the genomic and epigenomic landscape of Alzheimer’s disease in the temporal cortex, but can be used for a myriad of tissue regions and pathologies.  In addition to cross modality integration at the level of cell type, this pipeline uses Seurat's implementation of diagonalized CCA to match snATAC subclusters to snRNA subclusters.  The method of subcluster matching was validated using 10X genomics PBMC multiome data and found to be ~70% accurate (see `/atac_pipeline/7.cluster_link_validation/`).  Be sure to keep this in mind when interpreting downstream results.  

# config.R file

This file is sourced at the beginning of each script.  If working on multiple computers or HPCs, the file detects what machine R is running on based on the path to the working directory (`grepl(‘/hpc_path_identifier’, getwd()`).  Then, directory paths and library paths are set accordingly.  The config file also silently loads commonly used packages to avoid repetitive use of `library()` at the beginning of each script.

# Directory layout

R scripts in each subdirectory of the pipeline save all output files to a respective user-specified directory.  These directories can be specified in the config.R file. <br />
R scripts in this repository are located as follows: 

./config.R (this sets all user-specified paths from which to read/write files)

/rna_pipeline/ <br />
&emsp;1.preprocessing/ <br />
&emsp;2.integration/ <br />
&emsp;3.cell_type_annotation/ <br />
&emsp;4.deg_analysis/ <br />
<br />
  
/atac_pipeline/ <br />
&emsp;1.preprocessing/ <br />
&emsp;2.cell_type_annotation/ <br />
&emsp;3.integration/ <br />
&emsp;4.dap_analysis/ <br />
&emsp;5.cicero/ <br />
&emsp;6.motif_enrichment/ <br />
&emsp;7.cluster_link_validation/ <br />
