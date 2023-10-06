cond <- grepl('/hpc_path_identifier', getwd())
if (cond){
  source('path1/to/config.R')
} else {
  source('path2/to/config.R')
} ; rm(cond)

peak.path <- '\\path\\to\\recalled_peaks\\'
mspc.path <- '\\path\\to\\mspc_folder\\mspc.dll' 

setwd( paste0(atac_dap_analysis_objects, 'recalled_peaks/') )
cluster <- list.files(pattern = '_peaks.narrowPeak') 
cluster <- gsub(pattern = '.peaks.*_peaks.*', replacement = '', cluster) %>% unique



# access your OS command line shell through R
# windows:
system('cmd.exe', input = 'dotnet')
system('cmd.exe', input = paste0('mkdir ', peak.path, 'consensus_peaks\\', cluster, '_consensus_peaks\n'))

# mac/nix:
system('dotnet')
system(paste0('mkdir ', peak.path, 'consensus_peaks\\', cluster, '_consensus_peaks\n'))



for (i in 1:length(cluster)){
  message('\n',cluster[i])
  setwd(paste0(atac_dap_analysis_objects, 'recalled_peaks'))
  peak.list <- list.files(pattern = paste0(cluster[i], '.peaks.*.narrowPeak'))  #0 KB files removed 
  peak.list <- paste0('-i', ' ', peak.path, peak.list)
  message('found in ', length(peak.list), ' total datasets')
  
  out.dir <- paste0(peak.path, 'consensus_peaks\\', cluster[i], '_consensus_peaks')
  cmd <- paste('dotnet', mspc.path, paste(peak.list, collapse = ' '), 
               '-r bio -w 1e4 -s 1e-8','-c 2', '-o', out.dir, sep = ' ')
  
  # execute on windows
  system('cmd.exe', input = cmd)
  
  # execute on mac/nix
  # system(cmd)
  rm(cmd)
}








