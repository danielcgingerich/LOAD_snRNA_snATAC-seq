cond <- grepl('/hpc_path_identifier', getwd())
if (cond){
  source('path1/to/config.R')
} else {
  source('path2/to/config.R')
} ; rm(cond)
library(lmtest)
library(pbapply)
library(future)
library(Matrix.utils)

# count matrix is split into chunks and parallel computing used to speed up the process. 
# i = cluster index
# j = chunk index
i=1
slurm.array.id <- Sys.getenv('SLURM_ARRAY_TASK_ID')
j <- as.numeric(slurm.array.id)


message('loading objects')
setwd(atac_dap_analysis_objects)

# mtx <- readRDS('6.process_new_peaks___srt_consensus_processed_merged.rds')
# mtx <- mtx[['recalled.peaks']]@data
# saveRDS(mtx, '9.likelihood_ratio___normalized_count_matrix.rds')
data_use <- readRDS('9.likelihood_ratio___normalized_count_matrix.rds')

peaks <- readRDS('4.combine_consensus_peaks___cluster_specific_consensus_peaks_c=2.rds')
metadata <- readRDS('7.add_covarates___metadata_for_dap_analysis.rds')
data_use <- data_use[,rownames(metadata)]




#prep data\~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
message('\nprepping data')
cluster <- levels(as.factor(metadata$cluster_links))[i]
celltype <- gsub('[[:digit:]]', '', cluster)

metadata$cluster_links_diagnosis <- paste0(metadata$cluster_links, '_', metadata$diagnosis)
cells.1 <- rownames(metadata[metadata$cluster_links_diagnosis == paste0(cluster, '_LOAD'),])
cells.2 <- rownames(metadata[metadata$cluster_links_diagnosis == paste0(cluster, '_Normal'),])

group_info <- data.frame(row.names = c(cells.1, cells.2))
group_info$group <- ifelse(rownames(group_info) %in% cells.1, 0, 1)

# latent variables chosen based on correlation with principal components. 
# see "choose_covariates.R" in rna_pipeline/4.deg_analysis
latent_vars <- c('peak_region_fragments',
                 'percent.fragments.overlapping.any.targeted.region', 
                 paste0(cluster, '.proportion'), 
                 'age', 'sex', 'pmi')
latent_vars <- metadata[c(cells.1, cells.2), latent_vars]
latent_vars$sex <- ifelse(as.character(latent_vars$sex) == 'M', 0, 1)
latent_vars$percent.fragments.overlapping.any.targeted.region <- gsub('%', '', latent_vars$percent.fragments.overlapping.any.targeted.region)

if(all(rownames(latent_vars) == rownames(group_info))){
  message('rownames match. merging.')
  latent_vars <- apply(latent_vars, MARGIN = 2, as.numeric)
  metadata <- lapply(list(group_info, latent_vars), function(x){as(as.matrix(x), 'sparseMatrix')})
  metadata <- cbind(metadata[[1]], metadata[[2]])
}

peaks <- GRangesToString(peaks)[grep(cluster, peaks$peak.called.in)]
data_use <- data_use[peaks, rownames(group_info), drop = FALSE]
npeaks <- length(peaks)





#split mtx\~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
message('\nsplitting data into chunks')
split_rows <- split(rownames(data_use), ceiling(seq_along(rownames(data_use))/(npeaks/100)))
model_data <- pblapply(split_rows, function(x){data_use[x,,drop = F]})
model_data <- pblapply(model_data, FUN = t)
model_data <- pblapply(model_data, 
                       function(x){
                         x <- cbind(x, metadata)
                         return(x)
                       })
rm(npeaks)



#calculate p.val\~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
message('\ncalculating pvals')

npeaks <- ncol(model_data[[j]])-ncol(metadata)
model_data_subset <- model_data[[j]]
colnames(model_data_subset) <- gsub('-', '_', colnames(model_data_subset)) # or else glm() thinks '-' is subtraction
model_data_subset <- as.data.frame(model_data_subset)

#null model
group <- model_data_subset[,'group']
var1 <- model_data_subset[,colnames(latent_vars)[1]]
var2 <- model_data_subset[,colnames(latent_vars)[2]]
var3 <- model_data_subset[,colnames(latent_vars)[3]]
var4 <- model_data_subset[,colnames(latent_vars)[4]]
var5 <- model_data_subset[,colnames(latent_vars)[5]]
var6 <- model_data_subset[,colnames(latent_vars)[6]]
fmla2 <- group ~ var1 + var2 + var3 + var4 + var5 + var6
model2 <- glm(formula = fmla2, family = "binomial")

p_val <- pbsapply(
  X = 1:npeaks,
  FUN = function(x) {
    message(x)
    peak <- model_data_subset[,x]
    
    # experimental model
    fmla <- group ~ peak + var1 + var2 + var3 + var4 + var5 + var6
    model1 <- glm(formula = fmla, family = "binomial")
    
    # model comparison (LR test)
    lrtest <- lrtest(model1, model2)
    return(lrtest$Pr[2])
  }
)

names(p_val) <- colnames(model_data_subset)[1:npeaks]

#save data\~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
message('\ndone! saving...')
setwd(atac_dap_analysis_objects)
saveRDS(p_val, paste0('9.likelihood_ratio___', cluster, '_pvals',j, '.rds'))
