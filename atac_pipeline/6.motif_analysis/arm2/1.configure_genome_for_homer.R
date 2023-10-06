# this was done on *nix based compute cluster:
# must have homer installed: http://homer.ucsd.edu/homer/motif/
# in BASH, run this line before openning R session:
#     export PATH="path/to/homer/bin/:$PATH"
# 
# open up an R session


cond <- grepl('/hpc_path_identifier', getwd())
if (cond){
  source('path1/to/config.R')
} else {
  source('path2/to/config.R')
} ; rm(cond)

fasta <- paste0(ref_genome_path, 'GRCh38.primary_assembly.genome.fa')
gtf <- paste0(ref_genome_path, 'gencode.v32.primary_assembly.annotation.gtf')

cmd <- paste0('loadGenome.pl', ' ',
              '-name', ' ', 'grch38p13', ' ',
              '-org', ' ', 'human', ' ', 
              '-fasta', ' ', fasta, ' ', 
              '-gtf', ' ', gtf)
system(cmd)


