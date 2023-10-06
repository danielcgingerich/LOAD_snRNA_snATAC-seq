cond <- grepl('/hpc_path_identifier', getwd())
if (cond){
  source('path1/to/config.R')
} else {
  source('path2/to/config.R')
} ; rm(cond)
library(ggplot2) 
library(fitdistrplus) 
library(MASS) 
library(tidyr) 
library(gdata)
library(Seurat)
library(data.table)
library(EnvStats)
library(purrr)
library(dplyr)
library(sn)
library(matrixStats)
library(fmsb)
library(lme4)

#####################################
### Format and transform the data ###
#####################################

print('formatting and transforming the data...')
# split seurat object by cluster to use individually for DE analysis:
setwd(rna_deg_analysis_objects)
Oligo1 <- readRDS('Oligo1_for_DE.rds')
allgenes_JG <- t(as.data.frame(Oligo1@assays$RNA@counts))
allgenes_JG <- as.data.frame(allgenes_JG)
ngenes <- ncol(allgenes_JG)
allgenes_JG$orig.ident <- Oligo1@meta.data$orig.ident
allgenes_JG$sampID <- Oligo1@meta.data$sampID
allgenes_JG$Sex <- Oligo1@meta.data$Sex
allgenes_JG$age <- Oligo1@meta.data$age
allgenes_JG$PMI <- Oligo1@meta.data$PMI
allgenes_JG$seq.sat <- Oligo1@meta.data$seq.sat
allgenes_JG$wellKey <- paste(Oligo1@meta.data$sampID, rownames(Oligo1@meta.data), sep = "_")
rownames(allgenes_JG) <- allgenes_JG$wellKey
allgenes_JG$cell.type.number.proportion <- Oligo1@meta.data$cell.type.number.proportion
allgenes_JG <- allgenes_JG[,c(ngenes+7,ngenes+2,ngenes+1,ngenes+3,ngenes+4,ngenes+5,ngenes+6,ngenes+8,1:ngenes)]
genecounts_JG <- as.matrix(t(allgenes_JG[,c(-1,-2,-3,-4,-5,-6,-7,-8)]))
genecounts_JG <- log2(genecounts_JG + 1) #log2 transform
coldata_JG <- allgenes_JG[,1:8]
coldata_JG$orig.ident <- as.factor(coldata_JG$orig.ident)
coldata_JG$Sex <- as.factor(coldata_JG$Sex)
nrow(genecounts_JG) #36588 genes

#####################################
###        Filter the data        ###
#####################################

print('filtering the data...')

#Calculate %cells expressing each gene in each group
LOAD.cells <- as.matrix(t(allgenes_JG[allgenes_JG$orig.ident == 'LOAD', c(-1,-2,-3,-4,-5,-6,-7,-8)]))
ncol(LOAD.cells) #count LOAD cells
Normal.cells <- as.matrix(t(allgenes_JG[allgenes_JG$orig.ident == 'Normal', c(-1,-2,-3,-4,-5,-6,-7,-8)]))
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
genecounts_JG <- genecounts_JG[rownames(genecounts_JG)%in%genes.to.keep,]

#####################################
###               Run MAST        ###
#####################################

print('running MAST...')

genecounts_JG <- genecounts_JG[,rownames(coldata_JG)] #match order of data with metadata

fData_JG <- data.frame(primerid=rownames(genecounts_JG)) #feature data
sca_JG <- suppressMessages(MAST::FromMatrix(exprsArray=genecounts_JG, cData=coldata_JG, fData=fData_JG)) #makes a single-cell assay with expression data, metadata (e.g. Control_1_cell_1), and feature data
cdr2_JG <- colSums(SummarizedExperiment::assay(sca_JG)>0) #compute cellular detection rate (cdr)
SummarizedExperiment::colData(sca_JG)$ngeneson <- scale(cdr2_JG) #adds centered and scaled cdr to listData in colData of sca object
cond <-factor(SummarizedExperiment::colData(sca_JG)$orig.ident)
cond <-relevel(cond,"Normal") # set the reference level of orig.ident condition to be Normal
SummarizedExperiment::colData(sca_JG)$orig.ident <- cond
SummarizedExperiment::colData(sca_JG)$sampID <- factor(SummarizedExperiment::colData(sca_JG)$sampID) #set sampID to factor
SummarizedExperiment::colData(sca_JG)$Sex <- factor(SummarizedExperiment::colData(sca_JG)$Sex) #set Sex to factor

zlmCond_JG_prop <- suppressMessages(MAST::zlm(~ orig.ident + ngeneson + cell.type.number.proportion + seq.sat + age + Sex + PMI + (1 | sampID), sca_JG, method='glmer',ebayes = F, fitArgsD = list(nAGQ = 0), strictConvergence = FALSE)) #runs zero-inflated regression fitting a generalized linear mixed-effects model.

colnames(coef(zlmCond_JG_prop, 'D')) #check names of modeled coefficients

# Generate summary of results with likelihood ratio test
summaryCond_JG_prop <- suppressMessages(MAST::summary(zlmCond_JG_prop,
                                                      doLRT='orig.identLOAD'))

#format the results
summaryDt_JG_prop <- summaryCond_JG_prop$datatable

fcHurdle_JG_prop <- merge(summaryDt_JG_prop[contrast=='orig.identLOAD' & component=='H',.(primerid, `Pr(>Chisq)`)], 
                          summaryDt_JG_prop[contrast=='orig.identLOAD' & component=='logFC', .(primerid, coef, ci.hi, ci.lo)], 
                          by='primerid')

fcHurdle_JG_fdr_prop <- fcHurdle_JG_prop[,fdr:=p.adjust(`Pr(>Chisq)`, 'fdr')] #add fdr adjusted p-values

#Add % cells from each group expressing each gene
pct.cells.LOAD <- as.data.frame(pct.exp.LOAD)
pct.cells.LOAD <- subset(pct.cells.LOAD, rownames(pct.cells.LOAD)%in%fcHurdle_JG_fdr_prop$primerid)
pct.cells.LOAD <- tibble::rownames_to_column(pct.cells.LOAD, var = "primerid")
pct.cells.Normal <- as.data.frame(pct.exp.Normal)
pct.cells.Normal <- subset(pct.cells.Normal, rownames(pct.cells.Normal)%in%fcHurdle_JG_fdr_prop$primerid)
pct.cells.Normal <- tibble::rownames_to_column(pct.cells.Normal, var = "primerid")
fcHurdle_JG_fdr_prop_pct.cells <- merge(fcHurdle_JG_fdr_prop, pct.cells.LOAD, by.x = 'primerid', by.y = 'primerid')
fcHurdle_JG_fdr_prop_pct.cells <- merge(fcHurdle_JG_fdr_prop_pct.cells, pct.cells.Normal, by.x = 'primerid', by.y = 'primerid')
length(rownames(fcHurdle_JG_fdr_prop_pct.cells))

fcHurdle_JG_na_omit <- stats::na.omit(as.data.frame(fcHurdle_JG_fdr_prop_pct.cells)) #count genes omitted - NA values occur from convergence failures
length(rownames(fcHurdle_JG_na_omit))

fcHurdleSig_JG <- subset(fcHurdle_JG_na_omit, subset = fdr<0.05) #determine number of significant genes
length(rownames(fcHurdleSig_JG))

#Save work
setwd(rna_deg_analysis_objects)
write.csv(fcHurdle_JG_na_omit, 'Oligo1_all_genes.csv')
write.csv(fcHurdleSig_JG, 'Oligo1_sig_genes.csv')
saveRDS(summaryCond_JG_prop, 'Oligo1_LRT.rds')
saveRDS(zlmCond_JG_prop, 'Oligo1_zlm.rds')
