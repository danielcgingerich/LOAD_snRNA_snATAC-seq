# code provided by Zhaohui Man
cond <- grepl('/hpc_path_identifier', getwd())
if (cond){
  source('path1/to/config.R')
} else {
  source('path2/to/config.R')
} ; rm(cond)

library(devtools)
library(rlang)
library(ggplot2)
library(scater)
library(Seurat)
library(cowplot)
library(Matrix.utils)
library(edgeR)
library(dplyr)
library(magrittr)
library(Matrix)
library(purrr)
library(reshape2)
library(S4Vectors)
library(tibble)
library(SingleCellExperiment)
library(pheatmap)
library(apeglm)
library(png)
library(RColorBrewer)
library(fitdistrplus) 
library(MASS) 
library(tidyr) 
library(gdata)
library(data.table)
library(EnvStats)
library(sn)
library(matrixStats)
library(fmsb)
library(lme4)
library(nebula)

obj.path <- '~/AD/'
out.path <- '~/results/'
setwd (obj.path)
obj_list <- list.files(pattern = '_for_DE.rds')
obj_list
#i=1
for (i in 1:length(obj_list)) {
  tryCatch( {rna.obj <- obj_list[i]
  cluster.id <- gsub('_for_DE.rds', '', rna.obj)
  
  print (cluster.id)
  print (i)
  
  #####################################
  ### Format and transform the data ###
  #####################################
  setwd (obj.path)
  rna.obj <- readRDS(rna.obj)
  
  DefaultAssay(rna.obj) <- "RNA"
  rna.obj
  
  print('removing genes that are not expressed')
  #Remove genes that are not expressed
  data <- rna.obj@assays$RNA@data
  nrow(data) 
  genes.use <- rowSums(data) > 0
  genes.use <- as.data.frame(genes.use)
  head(genes.use)
  genes.use <- subset(genes.use, genes.use=='TRUE')
  genes.use <- as.vector(rownames(genes.use))
  head(genes.use)
  length(genes.use) 
  gc()
  
  rna.obj <- subset(rna.obj, features=genes.use)
  
  data.subset <-rna.obj@assays$RNA@data
  nrow(data.subset) 
  allgenes <- t(as.data.frame(rna.obj@assays$RNA@data))
  allgenes <- as.data.frame(allgenes)
  
  ngenes <- ncol(allgenes)
  allgenes$orig.ident <- rna.obj@meta.data$orig.ident
  allgenes$sampID <- rna.obj@meta.data$sampID
  allgenes$Sex <- rna.obj@meta.data$Sex
  allgenes$age <- rna.obj@meta.data$age
  allgenes$PMI <- rna.obj@meta.data$PMI
  allgenes$seq.sat <- rna.obj@meta.data$seq.sat
  
  allgenes$wellKey <- rownames(rna.obj@meta.data)
  rownames(allgenes) <- allgenes$wellKey
  allgenes$cell.type.number.proportion <- rna.obj@meta.data$cell.type.number.proportion
  allgenes <- allgenes[,c(ngenes+7,ngenes+2,ngenes+1,ngenes+3,ngenes+4,ngenes+5,ngenes+6,ngenes+8,1:ngenes)]
  
  if (nrow(allgenes) == 0 ) {next}
  
  genedata <- as.matrix(t(allgenes[,c(-1,-2,-3,-4,-5,-6,-7,-8)]))
  
  nrow(genedata) 
  
  
  coldata <- allgenes[,1:8]
  coldata$orig.ident <- as.factor(coldata$orig.ident)
  
  coldata$sampID <- as.factor(coldata$sampID)
  
  for (j in 5:8) {coldata[ ,j] <- scale (coldata[, j])}
  
  coldata[1:3, ]
  
  rm(rna.obj)
  #####################################
  ###        Filter the data        ###
  #####################################
  
  print('filtering the data...')
  
  #Calculate %cells expressing each gene in each group
  LOAD.cells <- as.matrix(t(allgenes[allgenes$orig.ident == 'LOAD', c(-1,-2,-3,-4,-5,-6,-7,-8)]))
  ncol(LOAD.cells) #count LOAD cells
  Normal.cells <- as.matrix(t(allgenes[allgenes$orig.ident == 'Normal', c(-1,-2,-3,-4,-5,-6,-7,-8)]))
  ncol(Normal.cells) #count Normal cells
  
  PercentAbove <- function(x, threshold){
    return(length(x = x[x > threshold]) / length(x = x))
  }
  
  pct.exp.LOAD <- apply(X = LOAD.cells, MARGIN = 1, FUN = PercentAbove, threshold = 0)
  pct.exp.Normal <- apply(X = Normal.cells, MARGIN = 1, FUN = PercentAbove, threshold = 0)
  
  # Filter out genes expressed in <10% of cells in 1 group 
  alpha.min <- pmax(pct.exp.LOAD, pct.exp.Normal)
  genes.to.keep <- names(x = which(x = alpha.min >= 0.1)) 
  length(genes.to.keep) #count genes to keep
  
  genedata <- genedata[rownames(genedata)%in%genes.to.keep,]
  
  
  genedata <- genedata[,rownames(coldata)] #match order of data with metadata
  
  fData <- data.frame(primerid=rownames(genedata)) #feature data
  rm(LOAD.cells, Normal.cells)
  rm(data, data.subset, allgenes)
  
  #run Nebula
  
  count <- genedata
  rm(genedata)
  count [1:3, 1:3]
  sid <- coldata$sampID
  sid <- as.numeric(gsub ('rna', '', sid))
  cov <- coldata
  cov$Sex[cov$Sex=='Male'] <- 0
  cov$Sex[cov$Sex=='Female'] <- 1
  cov$Sex <- as.numeric (cov$Sex)
  cov[1:3, ]
  prep <- cov[, c(4,5,6,7,8,3)]
  prep [1:3, ]
  
  dim(count)
  dim(cov)
  count <- as.matrix(count)
  df = model.matrix (~Sex+age+PMI+seq.sat+cell.type.number.proportion+orig.ident, data=prep)
  head(df)
  re = nebula(count,sid,pred=df,method='HL')
  re
  re$summary$logFC_orig.identNormal
  re$summary [1:3, ]
  result <- re$summary
  result$log2FC <- log2(exp(result$logFC_orig.identNormal)) #convert natural log to log2
  result[1:3, ]
  nrow(result)
  result$fdr <- p.adjust (result$p_orig.identNormal, method= 'fdr', n=nrow(result))
  result <- result[, c('gene','log2FC', 'fdr')]
  result_sig <- result[which(result$fdr < 0.05), ]
  result_sig[1:3, ]
  result <- as.data.frame(result)
  result <- result[order(result$fdr), ]
  result_sig <- as.data.frame(result_sig)
  result_sig <- result_sig[order(result_sig$fdr), ]
  #save work
  setwd(out.path)
  write.csv(result, paste0('LOAD.', cluster.id, '_all_genes_NEBULA.csv'))
  write.csv(result_sig, paste0('LOAD.', cluster.id, '_sig_genes_NEBULA.csv'))
  
  rm (coldata, count, cov,df,fData, prep, re, result, result_sig)
  gc() }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  
}

print(
  sessionInfo()
)
