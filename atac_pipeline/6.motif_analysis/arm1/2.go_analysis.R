#GO analysis
cond <- grepl('/hpc_path_identifier', getwd())
if (cond){
  source('path1/to/config.R')
} else {
  source('path2/to/config.R')
} ; rm(cond)
library(ggplot2)
library(cowplot)
library(data.table)
library(topGO)
library(org.Hs.eg.db)
library(forcats)
library(ggplot2)

# topGO version 2.42.0 was used 

# import all genes tested per cluster
setwd(rna_deg_analysis_objects)

genes_tested <- list.files(pattern = '[[:digit:]]_all_genes.csv$')
genes_tested <- lapply(genes_tested, 
                       function(x){
                         csv <- read.csv(x)
                         csv$cluster <- gsub('_.*', '', x)
                         csv <- as.data.table(csv)
                         return(csv)
                       }) %>% bind_rows()

# import double dap data
setwd(atac_motif_enrichment_objects)
ccre <- list.files(pattern = '1.double_daps___.*.csv$')
ccre <- lapply(ccre, 
               function(x){
                 dr <- gsub('.*___|_ccans.csv', '', x)
                 ccre <- read.csv(x)
                 ccre$ccan_direction <- dr
                 ccre <- as.data.table(ccre)
                 return(ccre)
               }) %>% bind_rows()
ccre$degs <- gsub('.*\n', '', ccre$degs)



# double daps GO -------------------------------------------------------------
ccre <- ccre[, paste(degs, collapse = ', '), by = cluster]
gene_set <- strsplit(ccre$V1, split = ', ')
gene_set <- lapply(1:length(gene_set),
                   function(x){
                     name <- ccre$cluster[x]
                     genes <- gene_set[[x]]
                     data.table(cluster = name, gene = unique(genes))
                   }) %>% bind_rows()

# limit GO to gene sets of at least size 3
filt <- gene_set[, .N, by = cluster][N >= 3, cluster]
gene_set <- gene_set[cluster %in% filt, ]

cluster_ix <- unique(gene_set$cluster)
results <- results_cc <- results_bp <- results_mf <- list()
enr_cc <- enr_bp <- enr_mf <- list()
topgo_bp <- topgo_cc <- topgo_mf <- list()
for (i in cluster_ix){
  message(i)
  genes_ix <- gene_set[cluster == i, unique(gene)]
  all_genes_ix <- genes_tested[cluster == i, primerid]
  go_input <- ifelse(all_genes_ix %in% genes_ix, 1, 0)
  names(go_input) <- all_genes_ix
  
  topgo_bp[[i]] <- new('topGOdata', ontology = 'BP', allGenes = go_input, 
                    geneSelectionFun = function(x){x==1}, nodeSize = 3,
                    annot = annFUN.org, mapping = 'org.Hs.eg.db', 
                    ID = 'symbol')
  topgo_cc[[i]] <- new('topGOdata', ontology = 'CC', allGenes = go_input, 
                       geneSelectionFun = function(x){x==1}, nodeSize = 3,
                       annot = annFUN.org, mapping = 'org.Hs.eg.db', 
                       ID = 'symbol')
  topgo_mf[[i]] <- new('topGOdata', ontology = 'MF', allGenes = go_input, 
                       geneSelectionFun = function(x){x==1}, nodeSize = 3,
                       annot = annFUN.org, mapping = 'org.Hs.eg.db', 
                       ID = 'symbol')
  
  enr_bp[[i]] <- runTest(topgo_bp[[i]], algorithm = 'parentchild', statistic = 'fisher')
  enr_cc[[i]] <- runTest(topgo_cc[[i]], algorithm = 'parentchild', statistic = 'fisher')
  enr_mf[[i]] <- runTest(topgo_mf[[i]], algorithm = 'parentchild', statistic = 'fisher')
  
  results_bp[[i]] <- GenTable(topgo_bp[[i]], fisher =  enr_bp[[i]], 
                      topNodes = length(topgo_bp[[i]]@graph@nodes), numChar = 100)
  results_mf[[i]] <- GenTable(topgo_mf[[i]], fisher =  enr_mf[[i]], 
                           topNodes = length(topgo_mf[[i]]@graph@nodes), numChar = 100)
  results_cc[[i]] <- GenTable(topgo_cc[[i]], fisher =  enr_cc[[i]], 
                           topNodes = length(topgo_cc[[i]]@graph@nodes), numChar = 100)
  
  results_cc[[i]] <- as.data.table(results_cc[[i]])
  results_cc[[i]]$cluster <- i
  results_cc[[i]]$ontology_type <- 'cc'
  results_cc[[i]]$fisher <- as.numeric(results_cc[[i]]$fisher)
  results_cc[[i]]$fdr <- p.adjust(results_cc[[i]]$fisher, method = 'fdr')
  
  results_bp[[i]] <- as.data.table(results_bp[[i]])
  results_bp[[i]]$cluster <- i
  results_bp[[i]]$ontology_type <- 'bp'
  results_bp[[i]]$fisher <- as.numeric(results_bp[[i]]$fisher)
  results_bp[[i]]$fdr <- p.adjust(results_bp[[i]]$fisher, method = 'fdr')
  
  results_mf[[i]] <- as.data.table(results_mf[[i]])
  results_mf[[i]]$cluster <- i
  results_mf[[i]]$ontology_type <- 'mf'
  results_mf[[i]]$fisher <- as.numeric(results_mf[[i]]$fisher)
  results_mf[[i]]$fdr <- p.adjust(results_mf[[i]]$fisher, method = 'fdr')
  
  results[[i]] <- rbind(results_cc[[i]], results_bp[[i]], results_mf[[i]])
  results[[i]]$fdr_all <- p.adjust(results[[i]]$fisher, method = 'fdr')
}

results <- bind_rows(results)
results$log2FE <- log2(results[, Significant/Expected])

setwd(atac_motif_enrichment_objects)
saveRDS(results, '2.go_analysis___results.rds')


# plot results ------------------------------------------------------------
setwd(atac_motif_enrichment_objects)
results <- readRDS('2.go_analysis___results.rds')
results <- results[Significant >= 3, ]
cluster_ix <- unique(results$cluster)

for (i in cluster_ix){
  message('plotting ', i)
  gg_go <- results[cluster == i, ]
  gg_go$fisher <- as.numeric(gg_go$fisher)
  gg_go <- gg_go[order(gg_go$fisher), .SD[1:10] , by = ontology_type]
  gg_go <- gg_go[fisher != 1, ]
  
  gg_mf <- results[cluster == i & ontology_type == 'mf', ]
  gg_mf$fisher <- as.numeric(gg_mf$fisher)
  gg_mf <- gg_mf[order(gg_mf$fisher), ][10:1, ]
  gg_mf$Term <- factor(gg_mf$Term, levels = gg_mf$Term)
  gg_mf <- gg_mf[fisher != 1, ]
  
  
  gg_bp <- results[cluster == i & ontology_type == 'bp', ]
  gg_bp$fisher <- as.numeric(gg_bp$fisher)
  gg_bp <- gg_bp[order(gg_bp$fisher), ][10:1, ]
  gg_bp$Term <- factor(gg_bp$Term, levels = gg_bp$Term)
  gg_bp <- gg_bp[fisher != 1, ]
  
  
  gg_cc <- results[cluster == i & ontology_type == 'cc', ]
  gg_cc$fisher <- as.numeric(gg_cc$fisher)
  gg_cc <- gg_cc[order(gg_cc$fisher), ][10:1, ]
  gg_cc$Term <- factor(gg_cc$Term, levels = gg_cc$Term)
  gg_cc <- gg_cc[fisher != 1, ]
  
  # x limit for plot
  vals <- c(gg_mf$fisher, gg_bp$fisher, gg_cc$fisher)
  vals <- sort(vals)
  vals[vals == 0] <- vals[vals != 0][1] # replace 0s with smallest non zero p val
  vals <- -log10(vals)
  x_limit <- xlim(0, max(vals)+0.1)
  
  # alt 2 
  salmon <- '#F8766D'
  green <- '#7CAE00'
  cyan <- '#00BFC4'
  
  custom_theme <- theme(axis.line.x = element_blank(), 
                        axis.ticks = element_blank(), 
                        axis.text.x = element_blank(), 
                        axis.title.x = element_blank(),
                        axis.title = element_text(face = 'bold'), 
                        axis.text.y = element_text(size = 10, angle = 0), 
                        axis.title.y = element_text(size = 10))
  ggplot() + 
    geom_bar(data = gg_mf, mapping = aes(x = -log10(fisher), y = Term), stat = 'identity', fill = salmon) + 
    ylab('Molecular Function') + theme_cowplot() + custom_theme + geom_vline(xintercept = -log10(0.05)) + x_limit + 
    ggtitle(i) -> mf
  
  ggplot() + 
    geom_bar(data = gg_cc, mapping = aes(x = -log10(fisher), y = Term), stat = 'identity', fill = green) + 
    ylab('Cellular Component') + theme_cowplot() + custom_theme + geom_vline(xintercept = -log10(0.05)) + x_limit -> cc
  
  ggplot() + 
    geom_bar(data = gg_bp, mapping = aes(x = -log10(fisher), y = Term), stat = 'identity', fill = cyan) + 
    ylab('Biological Process') + xlab( expression(bold(-log[10](p))) ) + 
    theme_cowplot() + theme(axis.title = element_text(face = 'bold', size = 12)) + 
    theme(axis.text.y = element_text(size = 10, angle = 0), axis.title.y = element_text(size = 10)) + 
    geom_vline(xintercept = -log10(0.05))  + x_limit-> bp
  
  img <- plot_grid(plotlist = list(mf, cc, bp), nrow = 3, align = 'v') ; img
  
  img2 <- bp
  setwd(atac_motif_enrichment_objects)
  plot_name <- paste0('2.go_analysis___', i, '_double_dap_go_plot___bp_cc_mf.png')
  png(filename = plot_name, width = 700, height = 1500, res = 100)
  print(img)
  dev.off()
  message('done! plot saved in ', getwd(), '/', plot_name)
}


