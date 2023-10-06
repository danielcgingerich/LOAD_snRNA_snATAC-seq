cond <- grepl('/hpc_path_identifier', getwd())
if (cond){
  source('path1/to/config.R')
} else {
  source('path2/to/config.R')
} ; rm(cond)

library(GenomeInfoDb)
library(patchwork)
library(Matrix.utils)

slurm.array.id <- Sys.getenv('SLURM_ARRAY_TASK_ID')                                               
i <- as.numeric(slurm.array.id)                                                                   

setwd(atac_preprocessing_objects)
atac <- list.files(pattern = '_seurat_qc_dim_reduced.rds')             
sample.id <- gsub(pattern = '.*___|_seurat.*', replacement = '', atac)
atac <- readRDS(atac[i])

# prepare ref genome ------------------------------------------------------


# gene names in hg38 are not unique. the seurat count matrix appends decimals to gene names in order
# to make them unique. this is so the count matrix can be indexed properly - you dont want to pull 
# the wrong gene out.  features.tsv contains the original gene name with its corresponding ensembl 
# id. below, i create a data frame: col1=appended gene name, col2=ensembl id, col3=original gene name.
# data.frame$col3 finds the hg38 gene with the same id. then, appended gene name is replaced.
# this makes the names in hg38 the same as rna object, so that the gene activity matrix will align
# with the rna data.

message("loading genome")
rna <- readRDS('set0_rna688.rds')                                                                               
features.tsv <- read.table('features.tsv')                                                             
hg38 <- rtracklayer::import(paste0(ref_genome_path, 'gencode.v32.primary_assembly.annotation.gtf'))                                                            
genome(hg38) <- "hg38"
seqlevelsStyle(hg38) <- "UCSC"
hg38 <- keepStandardChromosomes(hg38, pruning.mode = 'coarse')                                    
hg38$gene_biotype <- 'protein_coding'                                                             
hg38$cherry <- substr(hg38$gene_id, 1, 15)                                                        
hg38 <- hg38[hg38$cherry %in% features.tsv$V1,]                                                   
df <- data.frame("rownames"=rownames(rna), "feature"=features.tsv$V2, "ensg" = features.tsv$V1)   
df2 <- df[df[,1] != df[,2],]                                                                                                             

for (z in 1:nrow(df2)){
  hg38$gene_name[which(hg38$cherry == df2[z,3])] <- df2[z,1]                                      
}
hg38$gene_name[which(hg38$tag == 'PAR')] <- paste0(hg38$gene_name[which(hg38$tag == 'PAR')], 
                                                   '_PAR')                                        
length(unique(hg38$gene_name))                                                                    

features <- rownames(rna)
hg38 <- hg38[hg38$gene_name %in% features]                                                        

###
for (i in 1:8){
  message(i)
  setwd(atac_preprocessing_objects)
  atac <- list.files(pattern = '_seurat_qc_dim_reduced.rds')             
  sample.id <- gsub(pattern = '.*___|_seurat.*', replacement = '', atac)
  atac <- readRDS(atac[i])
###

Annotation(atac) <- hg38                                                                          
                                                                                                  

# construct gene activity matrix ------------------------------------------

mtx <- GeneActivity(atac, assay = 'peaks', features = features, max.width = Inf)                  
message(paste0("features present: ", nrow(mtx)))                                                  
atac[['gene.activity']] <- CreateAssayObject(counts = mtx)

message("saving...")

setwd(atac_cell_type_annotation_objects)
saveRDS(atac, paste0('1.gene_activity___', sample.id[i],  '_srt_qc_gene-activity.rds')) 
message('done! checking data...')

###
}
###



