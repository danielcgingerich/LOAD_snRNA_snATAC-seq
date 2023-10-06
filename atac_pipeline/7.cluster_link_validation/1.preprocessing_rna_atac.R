cond <- grepl('/hpc_path_identifier', getwd())
if (cond){
  source('path1/to/config.R')
} else {
  source('path2/to/config.R')
} ; rm(cond)
library(EnsDb.Hsapiens.v75)

# load the RNA and ATAC data
setwd(atac_cluster_link_validation_objects)
counts <- Read10X_h5("pbmc_granulocyte_sorted_10k_filtered_feature_bc_matrix.h5")
fragpath <- paste0(atac_cluster_link_validation_objects, "pbmc_granulocyte_sorted_10k_atac_fragments.tsv.gz")
rna <- CreateSeuratObject(counts = counts$`Gene Expression`, assay = "rna")
atac <- CreateChromatinAssay(counts = counts$Peaks, sep = c(":", "-"), fragments = fragpath)
atac <- CreateSeuratObject(atac, assay = 'peaks')


# rna preprocess ----------------------------------------------------------
rna[["percent.mt"]] <- PercentageFeatureSet(rna, pattern = "^MT-")

VlnPlot(rna, features = c("nFeature_rna", "nCount_rna", "percent.mt"), ncol = 3)


plot1 <- FeatureScatter(rna, feature1 = "nCount_rna", feature2 = "percent.mt")
plot2 <- FeatureScatter(rna, feature1 = "nCount_rna", feature2 = "nFeature_rna")
plot1 + plot2

rna <- subset(rna, subset = nFeature_rna >= 1000 & nFeature_rna < 5500 & percent.mt < 20)

rna <- NormalizeData(rna, normalization.method = "LogNormalize", scale.factor = 10000)

rna <- FindVariableFeatures(rna, selection.method = "vst", nfeatures = 2000)

all.genes <- rownames(rna)
rna <- ScaleData(rna, features = all.genes)
rna <- RunPCA(rna, features = VariableFeatures(object = rna))

DimHeatmap(rna, dims = 1:15, cells = 500, balanced = TRUE)
gc()

ElbowPlot(rna, ndims = 50)
DepthCor(rna, reduction = 'pca')
rna <- RunUMAP(rna, dims = 1:30)

rna <- FindNeighbors(rna, dims = 1:2, reduction = 'umap')
rna <- FindClusters(rna, resolution = 0.15)

DimPlot(rna, reduction = "umap", group.by = 'seurat_clusters', label = T) # 0 and 4 same cluster
rna$seurat_clusters[rna$seurat_clusters == 4] <- 0


# atac preprocess ---------------------------------------------------------
# extract gene annotations from EnsDb
annotations <- rtracklayer::import(paste0(ref_genome_path, 'gencode.v32.primary_assembly.annotation.gtf'))
Annotation(atac) <- annotations 

atac <- NucleosomeSignal(object = atac)
Annotation(atac)$gene_biotype <- Annotation(atac)$gene_type

fragpath <- paste0(atac_cluster_link_validation_objects, "pbmc_granulocyte_sorted_10k_atac_fragments.tsv.gz")
atac@assays$peaks@fragments[[1]]@path <- fragpath
atac <- TSSEnrichment(object = atac, fast = FALSE)

metadata <- read.csv('pbmc_granulocyte_sorted_10k_per_barcode_metrics.csv', row.names = 1)
atac <- AddMetaData(atac, metadata = metadata)
atac$pct_reads_in_peaks <- atac$atac_peak_region_fragments / atac$atac_fragments * 100
atac$high.tss <- ifelse(atac$TSS.enrichment > 2, 'High', 'Low')
TSSPlot(atac, group.by = 'high.tss') + NoLegend()
atac$nucleosome_group <- ifelse(atac$nucleosome_signal > 4, 'NS > 4', 'NS < 4')

cond1 <- atac$atac_peak_region_fragments > 1000 & 
  atac$atac_peak_region_fragments < quantile(atac$atac_peak_region_fragments, 0.99)
cond2 <- atac$pct_reads_in_peaks > 15
cond3 <- atac$nucleosome_signal < 4 
cond4 <- atac$TSS.enrichment > 2
cells_use <- colnames(atac)[cond1 & cond2 & cond3 & cond4]
atac <- subset(atac, cells = cells_use)

setwd(atac_cluster_link_validation_objects)
saveRDS(rna, '1.preprocessing_rna_atac___rna_preprocessed.rds')
saveRDS(atac, '1.preprocessing_rna_atac___atac_preprocessed.rds')




