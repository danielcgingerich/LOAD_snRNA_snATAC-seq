# module load macs2
cond <- grepl('/hpc_path_identifier', getwd())
if (cond){
  source('path1/to/config.R')
} else {
  source('path2/to/config.R')
} ; rm(cond)


slurm.array.task.id <- Sys.getenv('SLURM_ARRAY_TASK_ID')
i <- as.numeric(slurm.array.task.id)

setwd(atac_integration_objects)
atac_integrated <- readRDS('5.cluster_links___atac_integrated_and_linked.rds')
metadata <- atac_integrated@meta.data

sample_id <- unique(metadata$sample_id)
metadata <- metadata[metadata$sample_id == sample_id[i], ]
cluster <- metadata[, 'cluster_links']
cluster <- as.character(sort(unique(cluster)))

message(paste0('dataset: ', sample_id[i]))
message('\nbeginning peak call')
message('\n', length(cluster), ' total clusters')
message('\n', nrow(metadata), ' total cells')

frag_path <- 'path/to/fragment_files/directory/'
outs <- paste0(atac_dap_analysis_objects, 'recalled_peaks/')  
# run macs2 algorithm
for (j in 1:length(cluster)){
  message('\n', j, ' of ', length(cluster), ': ', cluster[j]) 
  cmd <- paste0('macs2 callpeak -t ', frag_path,'/', sample_id[i], '/', cluster[j], '_', sample_id[i], '.bed ',
                '-g 2.7e+09 -f BED --nomodel --extsize 200 --shift -100 -n ', cluster[j], '.peaks.',sample_id[i],' --outdir ',
                outs)
  message('call: \n',cmd)
  system(cmd)
  message('\ndone!')
  rm(cmd)
}



