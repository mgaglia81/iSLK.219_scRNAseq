library(dplyr)
library(Seurat)
library(cowplot)
library(ggplot2)

# If data needs to be loaded from RDS file
#islk <- readRDS('~/Documents/Genomics/islk-merged.rds')
#islk.list <- SplitObject(islk, split.by = 'orig.ident')

### Run this to load data from RDS file
#cell.combined <- readRDS('~/Documents/Genomics/reanalyze_with_insertion/cell-combined-final.rds')
#cell.islkh <- subset(cell.combined, orig.ident == 'islkh')
#cell.islkdox <- subset(cell.combined, orig.ident == 'islkdox')
#cell.islkIDN <- subset(cell.combined, orig.ident == 'islkIDN')
#cell.islkIFN <- subset(cell.combined, orig.ident == 'islkIFN')

#Colors for clusters according to Figure 5E
colors <- c('1' = '#D6D6D5', '10' = '#2661FD', '13' = '#0F80FF', '2' = '#76D5FE', '6' = '#33E866', '4' = '#2EAD1F',
            '3' = '#87E12C', '0' = '#FCFA43', '5' = '#EEC85F', '8' = '#FE9705', '7' = '#CA5915', '12' = '#FD2513',
            '9' = '#B04CF2', '11' = '#8336FF')

# Read 10X data
print('Reading 10X data...')
# Old analysis before adding GFP/RFP cassette
# islkh.data <- Read10X(data.dir = '~/Documents/Genomics/iSLKH_iSLK219_count/outs/filtered_feature_bc_matrix/')
# islkdox.data <- Read10X(data.dir = '~/Documents/Genomics/iSLK219_dox_count/outs/filtered_feature_bc_matrix/')
# islk_IDN.data <- Read10X(data.dir = '~/Documents/Genomics/iSLK219_IDN_count/outs/filtered_feature_bc_matrix/')
# islk_IFN.data <- Read10X(data.dir = '~/Documents/Genomics/iSLK219_IDN_IFN_count/outs/filtered_feature_bc_matrix/')

# New analysis with GFP/RFP cassette
islkh.data <- Read10X(data.dir = '~/Documents/Genomics/reanalyze_with_insertion/iSLKH_iSLK219_insert_count/outs/filtered_feature_bc_matrix/')
islkdox.data <- Read10X(data.dir = '~/Documents/Genomics/reanalyze_with_insertion/iSLK219_dox_insert_count/outs/filtered_feature_bc_matrix/')
islk_IDN.data <- Read10X(data.dir = '~/Documents/Genomics/reanalyze_with_insertion/iSLK219_IDN_insert_count/outs/filtered_feature_bc_matrix/')
islk_IFN.data <- Read10X(data.dir = '~/Documents/Genomics/reanalyze_with_insertion/iSLK219_IDN_IFN_insert_count/outs/filtered_feature_bc_matrix/')

print('Creating Seurat objects...')
islkh <- CreateSeuratObject(counts = islkh.data, project = 'islkh')
islkdox <- CreateSeuratObject(counts = islkdox.data, project = 'islkdox')
islkIDN <- CreateSeuratObject(counts = islk_IDN.data, project = 'islkIDN')
islkIFN <- CreateSeuratObject(counts = islk_IFN.data, project = 'islkIFN')

print('Normalizing data...')
islk.list <- c(islkh, islkdox, islkIDN, islkIFN)
islk.list <- lapply(X = islk.list, FUN = function(x) {
  x[['percent.mt']] <- PercentageFeatureSet(x, pattern = '^MT-')
  x <- subset(x, subset = percent.mt < 20)
  x <- subset(x, subset = nFeature_RNA < 7500)
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

# Unfiltered QC data
islkunf <- merge(islkh, y = c(islkdox, islkIDN, islkIFN),
              add.cell.ids = c('iSLKH_iSLK219','iSLK219+dox','iSLK219+dox+IDN','iSLK219+dox+IDN+IFN'),
              project = 'expt_283')

islkunf[['percent.mt']] <- PercentageFeatureSet(islkunf, pattern = '^MT-')

print('Integrating data...')
all.genes <- rownames(islkunf)
cell.anchors <- FindIntegrationAnchors(object.list = islk.list, dims = 1:30)
cell.combined <- IntegrateData(anchorset = cell.anchors, dims = 1:30)
DefaultAssay(cell.combined) <- "integrated"

# Add sample ID to metadata
cell.combined$sample <- 'temp'
cell.combined$sample[cell.combined$orig.ident=='islkh'] <- 'A'
cell.combined$sample[cell.combined$orig.ident=='islkdox'] <- 'B'
cell.combined$sample[cell.combined$orig.ident=='islkIDN'] <- 'C'
cell.combined$sample[cell.combined$orig.ident=='islkIFN'] <- 'D'

islkunf$sample <- 'temp'
islkunf$sample[islkunf$orig.ident=='islkh'] <- 'A'
islkunf$sample[islkunf$orig.ident=='islkdox'] <- 'B'
islkunf$sample[islkunf$orig.ident=='islkIDN'] <- 'C'
islkunf$sample[islkunf$orig.ident=='islkIFN'] <- 'D'

# Visualize QC metrics as a violin plot
print('Performing QC analysis...')
p1 <- VlnPlot(cell.combined, features = c('nFeature_RNA'), pt.size = 0.001, group.by = 'sample') + NoLegend() +
  stat_summary(fun.y = median, geom = 'point', size = 15, color = 'white', shape = 95)
p2 <- VlnPlot(cell.combined, features = c('nCount_RNA'), pt.size = 0.001, group.by = 'sample') + NoLegend() +
  stat_summary(fun.y = median, geom = 'point', size = 15, color = 'white', shape = 95)
p3 <- VlnPlot(cell.combined, features = c('percent.mt'), pt.size = 0.001, group.by = 'sample') + NoLegend() +
  stat_summary(fun.y = median, geom = 'point', size = 15, color = 'white', shape = 95)
qc.plot <- plot_grid(p1, p2, p3, ncol = 3)
postscript(file = paste('~/Documents/Genomics/reanalyze_with_insertion/qcplotv2.eps', sep = ''), width = 10, height = 8)
print(qc.plot)
dev.off()
qc.data <- matrix(c(length(cell.combined$sample), length(subset(cell.combined$sample, cell.combined$sample=='A')),
                 length(subset(cell.combined$sample, cell.combined$sample=='B')),
                 length(subset(cell.combined$sample, cell.combined$sample=='C')),
                 length(subset(cell.combined$sample, cell.combined$sample=='D')),
                 median(cell.combined$nFeature_RNA), median(cell.combined$nFeature_RNA[cell.combined$sample=='A']),
                 median(cell.combined$nFeature_RNA[cell.combined$sample=='B']),
                 median(cell.combined$nFeature_RNA[cell.combined$sample=='C']),
                 median(cell.combined$nFeature_RNA[cell.combined$sample=='D']),
                 median(cell.combined$nCount_RNA), median(cell.combined$nCount_RNA[cell.combined$sample=='A']),
                 median(cell.combined$nCount_RNA[cell.combined$sample=='B']),
                 median(cell.combined$nCount_RNA[cell.combined$sample=='C']),
                 median(cell.combined$nCount_RNA[cell.combined$sample=='D']),
                 median(cell.combined$percent.mt), median(cell.combined$percent.mt[cell.combined$sample=='A']),
                 median(cell.combined$percent.mt[cell.combined$sample=='B']),
                 median(cell.combined$percent.mt[cell.combined$sample=='C']),
                 median(cell.combined$percent.mt[cell.combined$sample=='D'])), ncol = 5, byrow = T)
colnames(qc.data) <- c('total', 'A', 'B', 'C', 'D')
rownames(qc.data) <- c('cell count', 'median unique RNAs', 'median reads', 'percentage reads mapped to mtDNA')
write.table(qc.data, file = '~/Documents/Genomics/reanalyze_with_insertion/qcdata.txt')

# Unfiltered qc data
p1 <- VlnPlot(islkunf, features = c('nFeature_RNA'), pt.size = 0.001, group.by = 'sample') + NoLegend() +
  stat_summary(fun.y = median, geom = 'point', size = 15, color = 'white', shape = 95) +
  geom_hline(yintercept = 7500, linetype = 'dashed', color = 'red')
p2 <- VlnPlot(islkunf, features = c('nCount_RNA'), pt.size = 0.001, group.by = 'sample') + NoLegend() +
  stat_summary(fun.y = median, geom = 'point', size = 15, color = 'white', shape = 95)
p3 <- VlnPlot(islkunf, features = c('percent.mt'), pt.size = 0.001, group.by = 'sample') + NoLegend() +
  stat_summary(fun.y = median, geom = 'point', size = 15, color = 'white', shape = 95) +
  geom_hline(yintercept = 20, linetype = 'dashed', color = 'red')
qc.unfiltered.plot <- plot_grid(p1, p2, p3, ncol = 3)
postscript(file = paste('~/Documents/Genomics/reanalyze_with_insertion/qcunfiltered.eps', sep = ''), width = 10, height = 8)
print(qc.unfiltered.plot)
dev.off()
qc.data.unfiltered <- matrix(c(length(islkunf$sample), length(subset(islkunf$sample, islkunf$sample=='A')),
                               length(subset(islkunf$sample, islkunf$sample=='B')),
                               length(subset(islkunf$sample, islkunf$sample=='C')),
                               length(subset(islkunf$sample, islkunf$sample=='D')),
                               median(islkunf$nFeature_RNA), median(islkunf$nFeature_RNA[islkunf$sample=='A']),
                               median(islkunf$nFeature_RNA[islkunf$sample=='B']),
                               median(islkunf$nFeature_RNA[islkunf$sample=='C']),
                               median(islkunf$nFeature_RNA[islkunf$sample=='D']),
                               median(islkunf$nCount_RNA), median(islkunf$nCount_RNA[islkunf$sample=='A']),
                               median(islkunf$nCount_RNA[islkunf$sample=='B']),
                               median(islkunf$nCount_RNA[islkunf$sample=='C']),
                               median(islkunf$nCount_RNA[islkunf$sample=='D']),
                               median(islkunf$percent.mt), median(islkunf$percent.mt[islkunf$sample=='A']),
                               median(islkunf$percent.mt[islkunf$sample=='B']),
                               median(islkunf$percent.mt[islkunf$sample=='C']),
                               median(islkunf$percent.mt[islkunf$sample=='D'])), ncol = 5, byrow = T)
colnames(qc.data.unfiltered) <- c('total', 'A', 'B', 'C', 'D')
rownames(qc.data.unfiltered) <- c('cell count', 'median unique RNAs', 'median reads', 'percentage reads mapped to mtDNA')
write.table(qc.data.unfiltered, file = '~/Documents/Genomics/reanalyze_with_insertion/qcdataunfiltered.txt')

# Run the standard workflow for visualization and clustering
print('Clearing objects from environment...')
rm(islkh.data, islkdox.data, islk_IDN.data, islk_IFN.data, islkh, islkdox, islkIDN, islkIFN, islk.list, islkunf, cell.anchors)
print('Scaling data...')
cell.combined <- ScaleData(cell.combined, verbose = TRUE, features = all.genes)
cell.combined <- RunPCA(cell.combined, npcs = 30, verbose = FALSE)
# t-SNE and Clustering
print('Performing dimensional reduction and clustering...')
cell.combined <- RunUMAP(cell.combined, reduction = "pca", dims = 1:30)
cell.combined <- FindNeighbors(cell.combined, reduction = "pca", dims = 1:30)
cell.combined <- FindClusters(cell.combined, resolution = 0.5)
# Visualization
p1 <- DimPlot(cell.combined, reduction = "umap", label = TRUE) + NoLegend()
p2 <- DimPlot(cell.combined, reduction = "umap", group.by = "sample") + NoLegend()
p3 <- plot_grid(p1, p2)

# Differential gene expression
print('Scaling RNA data...')
DefaultAssay(cell.combined) <- "RNA"
cell.combined <- ScaleData(cell.combined, verbose = TRUE, features = all.genes)
print('Performing differential gene expression analysis...')
# cluster14.markers <- FindConservedMarkers(cell.combined, ident.1 = 14, grouping.var = "orig.ident", verbose = FALSE)
# logfc is natural log
integration.markers <- FindAllMarkers(cell.combined, only.pos = T, min.pct = 0.25, logfc.threshold = 0.25)
cell.combined.top10 <- integration.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
cell.combined.top5 <- integration.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)
cell.combined.heatmap <- DoHeatmap(cell.combined, features = cell.combined.top5$gene, assay = 'integrated')

print('Saving Seurat object to RDS file...')
saveRDS(cell.combined, file = '~/Documents/Genomics/reanalyze_with_insertion/cell-combined.rds')

# Differential gene expression in IFNB expressing cells

feature.list <- read.table('~/Documents/Genomics/reanalyze_with_insertion/iSLK219_dox_insert_count/outs/filtered_feature_bc_matrix/features.tsv.gz', stringsAsFactors = F)
feature.list <- feature.list[,1:2]
colnames(feature.list) <- c('ensembl', 'gene')
feature.list$gene <- make.unique(feature.list$gene)

IFNBcell.barcodes <- rownames(subset(cell.combined, subset = IFNB1 > 0)@meta.data)
cell.IFNB <- cell.combined
Idents(cell.IFNB, cells = IFNBcell.barcodes) <- 'IFNB1cells'
IFNBcell.markers <- FindMarkers(cell.IFNB, ident.1 = 'IFNB1cells', min.pct = 0.25)
IFNBcell.posmarkers <- subset(IFNBcell.markers, subset = avg_logFC >= 0)
IFNBcell.negmarkers <- subset(IFNBcell.markers, subset = avg_logFC < 0)

IFNBcell.markers$gene <- rownames(IFNBcell.markers)
IFNBcell.posmarkers$gene <- rownames(IFNBcell.posmarkers)
IFNBcell.negmarkers$gene <- rownames(IFNBcell.negmarkers)

IFNBcell.markers <- merge(IFNBcell.markers, feature.list, by.x = 'gene', by.y = 'gene')
IFNBcell.posmarkers <- merge(IFNBcell.posmarkers, feature.list, by.x = 'gene', by.y = 'gene')
IFNBcell.negmarkers <- merge(IFNBcell.negmarkers, feature.list, by.x = 'gene', by.y = 'gene')

IFNBcell.markers <- IFNBcell.markers[order(IFNBcell.markers$p_val_adj),]
IFNBcell.posmarkers <- IFNBcell.posmarkers[order(IFNBcell.posmarkers$p_val_adj),]
IFNBcell.negmarkers <- IFNBcell.negmarkers[order(IFNBcell.negmarkers$p_val_adj),]

write.table(matrix(c(IFNBcell.markers$gene, IFNBcell.markers$ensembl, IFNBcell.markers$p_val_adj), ncol = 3),
            file = '~/Documents/Genomics/reanalyze_with_insertion/IFN_differential_expression/IFNBmarkers.csv', sep = ',',
            row.names = F, col.names = c('Gene', 'Ensemble', 'P-value adj'), quote = F)
write.table(matrix(c(IFNBcell.posmarkers$gene, IFNBcell.posmarkers$ensembl, IFNBcell.posmarkers$p_val_adj), ncol = 3),
            file = '~/Documents/Genomics/reanalyze_with_insertion/IFN_differential_expression/IFNBposmarkers.csv', sep = ',',
            row.names = F, col.names = c('Gene', 'Ensemble', 'P-value adj'), quote = F)
write.table(matrix(c(IFNBcell.negmarkers$gene, IFNBcell.negmarkers$ensembl, IFNBcell.negmarkers$p_val_adj), ncol = 3),
            file = '~/Documents/Genomics/reanalyze_with_insertion/IFN_differential_expression/IFNBnegmarkers.csv', sep = ',',
            row.names = F, col.names = c('Gene', 'Ensemble', 'P-value adj'), quote = F)

# Differential gene expression in IFNB expressing cells in clusters 9 and 11

IFNB3.barcodes <- rownames(subset(cell.combined, subset = ((seurat_clusters == 9 | seurat_clusters == 11) & IFNB1 > 0))@meta.data)
cell.IFNB <- cell.combined
cell.IFNB$IFNB1 <- ifelse(colnames(cell.IFNB) %in% IFNB3.barcodes, 'IFNB1+', 'IFNB1-')
Idents(cell.IFNB, cells = IFNB3.barcodes) <- 'IFNB1cells'
IFNB3.markers <- FindMarkers(cell.IFNB, ident.1 = 'IFNB1cells', min.pct = 0.25)
IFNB3.posmarkers <- subset(IFNB3.markers, subset = avg_logFC >= 0)
IFNB3.negmarkers <- subset(IFNB3.markers, subset = avg_logFC < 0)

IFNB3.markers$gene <- rownames(IFNB3.markers)
IFNB3.posmarkers$gene <- rownames(IFNB3.posmarkers)
IFNB3.negmarkers$gene <- rownames(IFNB3.negmarkers)

IFNB3.markers <- merge(IFNB3.markers, feature.list, by.x = 'gene', by.y = 'gene')
IFNB3.posmarkers <- merge(IFNB3.posmarkers, feature.list, by.x = 'gene', by.y = 'gene')
IFNB3.negmarkers <- merge(IFNB3.negmarkers, feature.list, by.x = 'gene', by.y = 'gene')

IFNB3.markers <- IFNB3.markers[order(IFNB3.markers$p_val_adj),]
IFNB3.posmarkers <- IFNB3.posmarkers[order(IFNB3.posmarkers$p_val_adj),]
IFNB3.negmarkers <- IFNB3.negmarkers[order(IFNB3.negmarkers$p_val_adj),]

write.table(matrix(c(IFNB3.markers$gene, IFNB3.markers$ensembl, IFNB3.markers$p_val_adj), ncol = 3),
            file = '~/Documents/Genomics/reanalyze_with_insertion/IFN_differential_expression/IFNB3markers.csv', sep = ',',
            row.names = F, col.names = c('Gene', 'Ensemble', 'P-value adj'), quote = F)
write.table(matrix(c(IFNB3.posmarkers$gene, IFNB3.posmarkers$ensembl, IFNB3.posmarkers$p_val_adj), ncol = 3),
            file = '~/Documents/Genomics/reanalyze_with_insertion/IFN_differential_expression/IFNB3posmarkers.csv', sep = ',',
            row.names = F, col.names = c('Gene', 'Ensemble', 'P-value adj'), quote = F)
write.table(matrix(c(IFNB3.negmarkers$gene, IFNB3.negmarkers$ensembl, IFNB3.negmarkers$p_val_adj), ncol = 3),
            file = '~/Documents/Genomics/reanalyze_with_insertion/IFN_differential_expression/IFNB3negmarkers.csv', sep = ',',
            row.names = F, col.names = c('Gene', 'Ensemble', 'P-value adj'), quote = F)

# Differential gene expression in IFNB expressing cells outside of cluster 9 and 11

IFNBother.barcodes <- rownames(subset(cell.combined, subset = ((seurat_clusters != 9 & seurat_clusters != 11) & IFNB1 > 0))@meta.data)
cell.IFNB <- cell.combined
cell.IFNB$IFNB1 <- ifelse(colnames(cell.IFNB) %in% IFNBother.barcodes, 'IFNB1+', 'IFNB1-')
Idents(cell.IFNB, cells = IFNBother.barcodes) <- 'IFNB1cells'
IFNBother.markers <- FindMarkers(cell.IFNB, ident.1 = 'IFNB1cells', min.pct = 0.25)
IFNBother.posmarkers <- subset(IFNBother.markers, subset = avg_logFC >= 0)
IFNBother.negmarkers <- subset(IFNBother.markers, subset = avg_logFC < 0)

IFNBother.markers$gene <- rownames(IFNBother.markers)
IFNBother.posmarkers$gene <- rownames(IFNBother.posmarkers)
IFNBother.negmarkers$gene <- rownames(IFNBother.negmarkers)

IFNBother.markers <- merge(IFNBother.markers, feature.list, by.x = 'gene', by.y = 'gene')
IFNBother.posmarkers <- merge(IFNBother.posmarkers, feature.list, by.x = 'gene', by.y = 'gene')
IFNBother.negmarkers <- merge(IFNBother.negmarkers, feature.list, by.x = 'gene', by.y = 'gene')

IFNBother.markers <- IFNBother.markers[order(IFNBother.markers$p_val_adj),]
IFNBother.posmarkers <- IFNBother.posmarkers[order(IFNBother.posmarkers$p_val_adj),]
IFNBother.negmarkers <- IFNBother.negmarkers[order(IFNBother.negmarkers$p_val_adj),]

write.table(matrix(c(IFNBother.markers$gene, IFNBother.markers$ensembl, IFNBother.markers$p_val_adj), ncol = 3),
            file = '~/Documents/Genomics/reanalyze_with_insertion/IFN_differential_expression/IFNBothermarkers.csv', sep = ',',
            row.names = F, col.names = c('Gene', 'Ensemble', 'P-value adj'), quote = F)
write.table(matrix(c(IFNBother.posmarkers$gene, IFNBother.posmarkers$ensembl, IFNBother.posmarkers$p_val_adj), ncol = 3),
            file = '~/Documents/Genomics/reanalyze_with_insertion/IFN_differential_expression/IFNBotherposmarkers.csv', sep = ',',
            row.names = F, col.names = c('Gene', 'Ensemble', 'P-value adj'), quote = F)
write.table(matrix(c(IFNBother.negmarkers$gene, IFNBother.negmarkers$ensembl, IFNBother.negmarkers$p_val_adj), ncol = 3),
            file = '~/Documents/Genomics/reanalyze_with_insertion/IFN_differential_expression/IFNBothernegmarkers.csv', sep = ',',
            row.names = F, col.names = c('Gene', 'Ensemble', 'P-value adj'), quote = F)

#Differential gene expression in IFNB expressing cells in the antiIFN sample

IFNBIFNcell.barcodes <- rownames(subset(cell.islkIFN, subset = IFNB1 > 0)@meta.data)
cell.IFNBIFN <- cell.islkIFN
cell.IFNBIFN$IFNB1 <- ifelse(colnames(cell.IFNBIFN) %in% IFNBIFNcell.barcodes, 'IFNB1+', 'IFNB1-')
Idents(cell.IFNBIFN, cells = IFNBIFNcell.barcodes) <- 'IFNB1cells'
IFNBIFNcell.markers <- FindMarkers(cell.IFNBIFN, ident.1 = 'IFNB1cells', min.pct = 0.25)
IFNBIFNcell.posmarkers <- subset(IFNBIFNcell.markers, subset = avg_logFC >= 0)
IFNBIFNcell.negmarkers <- subset(IFNBIFNcell.markers, subset = avg_logFC < 0)

IFNBIFNcell.markers$gene <- rownames(IFNBIFNcell.markers)
IFNBIFNcell.posmarkers$gene <- rownames(IFNBIFNcell.posmarkers)
IFNBIFNcell.negmarkers$gene <- rownames(IFNBIFNcell.negmarkers)

IFNBIFNcell.markers <- merge(IFNBIFNcell.markers, feature.list, by.x = 'gene', by.y = 'gene')
IFNBIFNcell.posmarkers <- merge(IFNBIFNcell.posmarkers, feature.list, by.x = 'gene', by.y = 'gene')
IFNBIFNcell.negmarkers <- merge(IFNBIFNcell.negmarkers, feature.list, by.x = 'gene', by.y = 'gene')

IFNBIFNcell.markers <- IFNBIFNcell.markers[order(IFNBIFNcell.markers$p_val_adj),]
IFNBIFNcell.posmarkers <- IFNBIFNcell.posmarkers[order(IFNBIFNcell.posmarkers$p_val_adj),]
IFNBIFNcell.negmarkers <- IFNBIFNcell.negmarkers[order(IFNBIFNcell.negmarkers$p_val_adj),]

write.table(matrix(c(IFNBIFNcell.markers$gene, IFNBIFNcell.markers$ensembl, IFNBIFNcell.markers$p_val_adj), ncol = 3),
            file = '~/Documents/Genomics/reanalyze_with_insertion/IFN_differential_expression/IFNBIFNmarkers.csv', sep = ',',
            row.names = F, col.names = c('Gene', 'Ensemble', 'P-value adj'), quote = F)
write.table(matrix(c(IFNBIFNcell.posmarkers$gene, IFNBIFNcell.posmarkers$ensembl, IFNBIFNcell.posmarkers$p_val_adj), ncol = 3),
            file = '~/Documents/Genomics/reanalyze_with_insertion/IFN_differential_expression/IFNBIFNposmarkers.csv', sep = ',',
            row.names = F, col.names = c('Gene', 'Ensemble', 'P-value adj'), quote = F)
write.table(matrix(c(IFNBIFNcell.negmarkers$gene, IFNBIFNcell.negmarkers$ensembl, IFNBIFNcell.negmarkers$p_val_adj), ncol = 3),
            file = '~/Documents/Genomics/reanalyze_with_insertion/IFN_differential_expression/IFNBIFNnegmarkers.csv', sep = ',',
            row.names = F, col.names = c('Gene', 'Ensemble', 'P-value adj'), quote = F)

# Differential gene expression in IFNB expressing cells in cluster 9 and 11 of the antiIFN sample

IFNBIFN3.barcodes <- rownames(subset(cell.islkIFN, subset = ((seurat_clusters == 9 | seurat_clusters == 11) & IFNB1 >0))@meta.data)
cell.IFNBIFN <- cell.islkIFN
cell.IFNBIFN$IFNB1 <- ifelse(colnames(cell.IFNBIFN) %in% IFNBIFN3.barcodes, 'IFNB1+', 'IFNB1-')
Idents(cell.IFNBIFN, cells = IFNBIFN3.barcodes) <- 'IFNB1cells'
IFNBIFN3.markers <- FindMarkers(cell.IFNBIFN, ident.1 = 'IFNB1cells', min.pct = 0.25)
IFNBIFN3.posmarkers <- subset(IFNBIFN3.markers, subset = avg_logFC >= 0)
IFNBIFN3.negmarkers <- subset(IFNBIFN3.markers, subset = avg_logFC < 0)

IFNBIFN3.markers$gene <- rownames(IFNBIFN3.markers)
IFNBIFN3.posmarkers$gene <- rownames(IFNBIFN3.posmarkers)
IFNBIFN3.negmarkers$gene <- rownames(IFNBIFN3.negmarkers)

IFNBIFN3.markers <- merge(IFNBIFN3.markers, feature.list, by.x = 'gene', by.y = 'gene')
IFNBIFN3.posmarkers <- merge(IFNBIFN3.posmarkers, feature.list, by.x = 'gene', by.y = 'gene')
IFNBIFN3.negmarkers <- merge(IFNBIFN3.negmarkers, feature.list, by.x = 'gene', by.y = 'gene')

IFNBIFN3.markers <- IFNBIFN3.markers[order(IFNBIFN3.markers$p_val_adj),]
IFNBIFN3.posmarkers <- IFNBIFN3.posmarkers[order(IFNBIFN3.posmarkers$p_val_adj),]
IFNBIFN3.negmarkers <- IFNBIFN3.negmarkers[order(IFNBIFN3.negmarkers$p_val_adj),]

write.table(matrix(c(IFNBIFN3.markers$gene, IFNBIFN3.markers$ensembl, IFNBIFN3.markers$p_val_adj), ncol = 3),
            file = '~/Documents/Genomics/reanalyze_with_insertion/IFN_differential_expression/IFNBIFN3markers.csv', sep = ',',
            row.names = F, col.names = c('Gene', 'Ensemble', 'P-value adj'), quote = F)
write.table(matrix(c(IFNBIFN3.posmarkers$gene, IFNBIFN3.posmarkers$ensembl, IFNBIFN3.posmarkers$p_val_adj), ncol = 3),
            file = '~/Documents/Genomics/reanalyze_with_insertion/IFN_differential_expression/IFNBIFN3posmarkers.csv', sep = ',',
            row.names = F, col.names = c('Gene', 'Ensemble', 'P-value adj'), quote = F)
write.table(matrix(c(IFNBIFN3.negmarkers$gene, IFNBIFN3.negmarkers$ensembl, IFNBIFN3.negmarkers$p_val_adj), ncol = 3),
            file = '~/Documents/Genomics/reanalyze_with_insertion/IFN_differential_expression/IFNBIFN3negmarkers.csv', sep = ',',
            row.names = F, col.names = c('Gene', 'Ensemble', 'P-value adj'), quote = F)

# Differential gene expression in IFNB expressing cells in cluster 9 and 11 of the IDN sample

IFNBIDN3.barcodes <- rownames(subset(cell.islkIDN, subset = ((seurat_clusters == 9 | seurat_clusters == 11) & IFNB1 >0))@meta.data)
cell.IFNBIDN <- subset(cell.islkIDN, subset = (seurat_clusters == 9 | seurat_clusters == 11))
cell.IFNBIDN$IFNB1exp <- ifelse(colnames(cell.IFNBIDN) %in% IFNBIDN3.barcodes, 'IFNB1+', 'IFNB1-')
#Idents(cell.IFNBIDN, cells = IFNBIDN3.barcodes) <- 'IFNB1cells'
IFNBIDN3.markers <- FindMarkers(cell.IFNBIDN, ident.1 = 'IFNB1+', group.by = 'IFNB1exp', min.pct = 0.25)
IFNBIDN3.posmarkers <- subset(IFNBIDN3.markers, subset = avg_logFC >= 0)
IFNBIDN3.negmarkers <- subset(IFNBIDN3.markers, subset = avg_logFC < 0)

IFNBIDN3.markers$gene <- rownames(IFNBIDN3.markers)
IFNBIDN3.posmarkers$gene <- rownames(IFNBIDN3.posmarkers)
IFNBIDN3.negmarkers$gene <- rownames(IFNBIDN3.negmarkers)

IFNBIDN3.markers <- merge(IFNBIDN3.markers, feature.list, by.x = 'gene', by.y = 'gene')
IFNBIDN3.posmarkers <- merge(IFNBIDN3.posmarkers, feature.list, by.x = 'gene', by.y = 'gene')
IFNBIDN3.negmarkers <- merge(IFNBIDN3.negmarkers, feature.list, by.x = 'gene', by.y = 'gene')

IFNBIDN3.markers <- IFNBIDN3.markers[order(IFNBIDN3.markers$p_val_adj),]
IFNBIDN3.posmarkers <- IFNBIDN3.posmarkers[order(IFNBIDN3.posmarkers$p_val_adj),]
IFNBIDN3.negmarkers <- IFNBIDN3.negmarkers[order(IFNBIDN3.negmarkers$p_val_adj),]

write.table(matrix(c(IFNBIDN3.markers$gene, IFNBIDN3.markers$ensembl, IFNBIDN3.markers$p_val_adj), ncol = 3),
            file = '~/Documents/Genomics/reanalyze_with_insertion/IFN_differential_expression/IFNBIDN3markersnew.csv', sep = ',',
            row.names = F, col.names = c('Gene', 'Ensemble', 'P-value adj'), quote = F)
write.table(matrix(c(IFNBIDN3.posmarkers$gene, IFNBIDN3.posmarkers$ensembl, IFNBIDN3.posmarkers$p_val_adj), ncol = 3),
            file = '~/Documents/Genomics/reanalyze_with_insertion/IFN_differential_expression/IFNBIDN3posmarkersnew.csv', sep = ',',
            row.names = F, col.names = c('Gene', 'Ensemble', 'P-value adj'), quote = F)
write.table(matrix(c(IFNBIDN3.negmarkers$gene, IFNBIDN3.negmarkers$ensembl, IFNBIDN3.negmarkers$p_val_adj), ncol = 3),
            file = '~/Documents/Genomics/reanalyze_with_insertion/IFN_differential_expression/IFNBIDN3negmarkersnew.csv', sep = ',',
            row.names = F, col.names = c('Gene', 'Ensemble', 'P-value adj'), quote = F)

marker.features <- c('IFNB1', 'IFNL1', 'RHEBL1', 'PMAIP1', 'IFIT2', 'CCL5', 'HIP1R', 'SPRY2', 'FZD4', 'IFIT3', 'ZC3HAV1', 'OASL', 'IL6',
                     'ZFP36L2', 'BBC3', 'CDKN2C', 'CITED2', 'IL27', 'NEURL3', 'SIX1', 'PPM1K', 'SPAG9', 'KCNV1', 'PELI1')

DoHeatmap(cell.IFNBIDN, features = marker.features, size = 3.2, draw.lines = T) + scale_fill_gradient2(low = 'blue', mid = 'white', high = 'red', midpoint = 0)

# Differential gene expression in IFNB expressing cells outside of cluster 9 and 11 of the antiIFN sample

IFNBIFNother.barcodes <- rownames(subset(cell.islkIFN, subset = ((seurat_clusters != 9 & seurat_clusters != 11) & IFNB1 >0))@meta.data)
cell.IFNBIFN <- cell.islkIFN
cell.IFNBIFN$IFNB1 <- ifelse(colnames(cell.IFNBIFN) %in% IFNBIFNother.barcodes, 'IFNB1+', 'IFNB1-')
Idents(cell.IFNBIFN, cells = IFNBIFNother.barcodes) <- 'IFNB1cells'
IFNBIFNother.markers <- FindMarkers(cell.IFNBIFN, ident.1 = 'IFNB1cells', min.pct = 0.25)
IFNBIFNother.posmarkers <- subset(IFNBIFNother.markers, subset = avg_logFC >= 0)
IFNBIFNother.negmarkers <- subset(IFNBIFNother.markers, subset = avg_logFC < 0)

IFNBIFNother.markers$gene <- rownames(IFNBIFNother.markers)
IFNBIFNother.posmarkers$gene <- rownames(IFNBIFNother.posmarkers)
IFNBIFNother.negmarkers$gene <- rownames(IFNBIFNother.negmarkers)

IFNBIFNother.markers <- merge(IFNBIFNother.markers, feature.list, by.x = 'gene', by.y = 'gene')
IFNBIFNother.posmarkers <- merge(IFNBIFNother.posmarkers, feature.list, by.x = 'gene', by.y = 'gene')
IFNBIFNother.negmarkers <- merge(IFNBIFNother.negmarkers, feature.list, by.x = 'gene', by.y = 'gene')

IFNBIFNother.markers <- IFNBIFNother.markers[order(IFNBIFNother.markers$p_val_adj),]
IFNBIFNother.posmarkers <- IFNBIFNother.posmarkers[order(IFNBIFNother.posmarkers$p_val_adj),]
IFNBIFNother.negmarkers <- IFNBIFNother.negmarkers[order(IFNBIFNother.negmarkers$p_val_adj),]

write.table(matrix(c(IFNBIFNother.markers$gene, IFNBIFNother.markers$ensembl, IFNBIFNother.markers$p_val_adj), ncol = 3),
            file = '~/Documents/Genomics/reanalyze_with_insertion/IFN_differential_expression/IFNBIFNothermarkers.csv', sep = ',',
            row.names = F, col.names = c('Gene', 'Ensemble', 'P-value adj'), quote = F)
write.table(matrix(c(IFNBIFNother.posmarkers$gene, IFNBIFNother.posmarkers$ensembl, IFNBIFNother.posmarkers$p_val_adj), ncol = 3),
            file = '~/Documents/Genomics/reanalyze_with_insertion/IFN_differential_expression/IFNBIFNotherposmarkers.csv', sep = ',',
            row.names = F, col.names = c('Gene', 'Ensemble', 'P-value adj'), quote = F)
write.table(matrix(c(IFNBIFNother.negmarkers$gene, IFNBIFNother.negmarkers$ensembl, IFNBIFNother.negmarkers$p_val_adj), ncol = 3),
            file = '~/Documents/Genomics/reanalyze_with_insertion/IFN_differential_expression/IFNBIFNothernegmarkers.csv', sep = ',',
            row.names = F, col.names = c('Gene', 'Ensemble', 'P-value adj'), quote = F)


# Write to table the number of cells per sample by cluster

for (i in 0:(max(as.numeric(cell.combined$seurat_clusters))-1)) {
  print(paste0('Writing txt for cluster ', as.character(i)))
  subset <- c(length(subset(cell.combined$seurat_clusters, cell.combined$seurat_clusters==i&cell.combined$sample=='A')),
             length(subset(cell.combined$seurat_clusters, cell.combined$seurat_clusters==i&cell.combined$sample=='B')),
             length(subset(cell.combined$seurat_clusters, cell.combined$seurat_clusters==i&cell.combined$sample=='C')),
             length(subset(cell.combined$seurat_clusters, cell.combined$seurat_clusters==i&cell.combined$sample=='D')))
  subset <- matrix(c('islkH,', 'islkdox,', 'islkIDN,', 'islkIFN,', subset), nrow = 4, ncol = 2)
  write.table(subset, file = paste('~/Documents/Genomics/reanalyze_with_insertion/Cluster_Counts/', as.character(i),'.csv', sep = ''),
              row.names = F, col.names = F, quote = F)
}

ISGs <- c('IRF7','CXCL10','IFIT2','IDO1','CXCL11','GBP1','IFIT1','RSAD2','OAS1','IFIT3','USP18','IFIT5','DDX58','IFIH1','NMI','TNFSF10',
          'AIM2','CCL5','IRF9','IFI16','DDX60','RTP4','OASL','MX1','SP110','ZC3HAV1','TLR3','FFAR2','DEC1','OAS2','GBP2','IL15',
          'PMAIP1','IFI44', 'IFNB1','ORF34','ORFK12','ORF42','GAPDH')
KSHVORFs <- c('ORFK1', 'ORF4', 'ORF6', 'ORF11', 'vIL6', 'ORFK3', 'ORF70', 'ORFK4', 'ORFK4.1', '1.4kb', 'ORFK5', 'ORFK6', 'PAN',
              'ORF16', 'ORF17.5', 'ORF18', 'ORF34', 'ORF35', 'ORF37', 'ORF38', 'ORF39', 'ORF45', 'ORF46', 'ORF47', 'ORF50', 'ORFK8',
              'ORF57', 'ORF58', 'ORF59', 'ORF60', 'ORF61', 'ORFK12', 'ORF71', 'ORF72', 'ORF73', 'ORF8', 'ORF9', 'ORF10', 'ORFK3A',
              'ORF17', 'ORF21', 'ORF22', 'ORF25', 'ORF26', 'ORF27', 'ORF28', 'ORF30', 'ORF31', 'ORF33', 'ORF42', 'ORF43', 'ORF44',
              'ORF45.1', 'ORFK8.1', 'ORF52', 'ORF53', 'ORF54', 'ORF55', 'ORF62', 'ORF65', 'ORF67', 'ORF67.5', 'ORF69', 'ORFK14',
              'ORF74', 'ORF75')

KSHVORFs.update <- c('ORFK1', 'ORF4', 'ORF6', 'ORF11', 'vIL6', 'ORFK3', 'ORF70', 'ORFK4', 'ORFK4.1', '1.4kb', 'ORFK5', 'ORFK6', 'PAN',
                     'ORF16', 'ORF17.5', 'ORF34', 'ORF35', 'ORF38', 'ORF39', 'ORF45', 'ORF46', 'ORF47', 'ORF50', 'ORFK8',
                     'ORF57', 'ORF58', 'ORF59', 'ORF60', 'ORF61', 'ORFK12', 'ORF71', 'ORF72', 'ORF73', 'ORF8', 'ORFK3A',
                     'ORF17', 'ORF21', 'ORF25', 'ORF26', 'ORF28', 'ORF33', 'ORF42', 'ORF45.1', 'ORFK8.1', 'ORF52', 'ORF53', 'ORF55',
                     'ORF62', 'ORF65', 'ORF69', 'ORFK14', 'ORF74', 'ORF75')

KSHVORFs.detected <- c('GFP', 'DSRED', 'PAC', 'vIL6', 'PAN', 'ORFK4.1', 'ORFK12', 'ORF73', 'ORF50', 'ORF57', 'ORFK1', 'ORF6', '1.4kb', 'ORFK5', 'ORFK6', 'ORF16',
                       'ORF34', 'ORF39', 'ORF8', 'ORF17', 'ORF21', 'ORF25', 'ORF28', 'ORF42', 'ORF53', 'ORF55', 'ORF62', 'ORF65', 'ORF69', 'ORF75')

KSHV.latent <- c('GFP', 'PAC', 'ORF73')
KSHV.early <- c('vIL6', 'ORF50','ORF57','ORFK1','ORF6','ORFK5', 'ORFK6', 'ORF16', 'ORF34', 'ORF39', 'ORF8', 'ORF53', 'ORF55', 'ORF69')
KSHV.intermediate <- c('ORFK4.1', 'ORF42', 'ORF75')
KSHV.late <- c('DSRED','ORFK12', 'ORF17', 'ORF21', 'ORF25', 'ORF28', 'ORF62', 'ORF65')

#reorder by viral stage --> Reporters, Latent, early, late
KSHVORFs.detected <- c('GFP', 'DSRED', 'PAC', 'ORF73', 'vIL6', 'PAN', 'ORFK4.1', 'ORFK12', 'ORF50', 'ORF57', 'ORFK1', 'ORF6', '1.4kb', 'ORFK5', 'ORFK6', 'ORF16',
                       'ORF34', 'ORF39', 'ORF8', 'ORF17', 'ORF21', 'ORF25', 'ORF28', 'ORF42', 'ORF53', 'ORF55', 'ORF62', 'ORF65', 'ORF69', 'ORF75')

IFNgenes.list <- c('DSRED','IFNB1', 'IFNL1', 'IFNL2', 'IFNG', 'IFNAR1', 'IFNAR2', 'TYK2', 'JAK1', 'JAK2', 'JAK3', 'STAT1', 'STAT2', 'IRF9',
                   'CGAS', 'IFI16', 'TMEM173', 'DDX58', 'MAVS', 'TLR3', 'TICAM1', 'TBK1', 'IKBKE', 'IRF3', 'IRF7',
                   'CYLD', 'RNF125', 'IFIH1', 'DDX25', 'TRIM25', 'TRADD', 'SIKE1', 'TRAF3', 'TANK', 'CASP8', 'CASP10', 'NFKBIB', 'NFKB1', 'DDX3X',
                   'ISG15', 'IFIT1', 'CCL5', 'IFIT2', 'ISG20', 'MYD88', 'OAS1', 'MX1', 'OASL', 'IFI6')

KSHV.reporters <- c('GFP', 'DSRED', 'PAC')
KSHV.latent2 <- c('ORFK12', 'ORF73')
KSHV.early2 <- c('1.4kb', 'PAN', 'vIL6', 'ORFK1', 'ORFK4.1', 'ORFK5', 'ORFK6', 'ORF6', 'ORF8', 'ORF16', 'ORF34', 'ORF39', 'ORF50', 'ORF53', 'ORF55', 'ORF57', 'ORF69')
KSHV.late2 <- c('ORF17', 'ORF21', 'ORF25', 'ORF28', 'ORF42', 'ORF62', 'ORF65', 'ORF75')
KSHVORFs.detected <- c(KSHV.reporters, KSHV.latent2, KSHV.early2, KSHV.late2)

# Average based on the expression within just the individual samples
cell.islkh <- subset(cell.combined, orig.ident == 'islkh')
cell.islkdox <- subset(cell.combined, orig.ident == 'islkdox')
cell.islkIDN <- subset(cell.combined, orig.ident == 'islkIDN')
cell.islkIFN <- subset(cell.combined, orig.ident == 'islkIFN')

islkh.averages <- AverageExpression(cell.islkh, return.seurat = T)
islkdox.averages <- AverageExpression(cell.islkdox, return.seurat = T)
islkIDN.averages <- AverageExpression(cell.islkIDN, return.seurat = T)
islkIFN.averages <- AverageExpression(cell.islkIFN, return.seurat = T)

# Combined averages of the clusters for the samples
cell.averages <- AverageExpression(cell.combined, return.seurat = T)

# Averages of the clusters separated by each sample
cell.averages2 <- AverageExpression(cell.combined, return.seurat = T, add.ident = c('orig.ident', 'seurat_clusters'))

cell.averages2$sample <- c('islkh', 'islkh', 'islkh', 'islkh', 'islkh', 'islkh', 'islkh', 'islkh', 'islkh', 'islkh', 'islkh', 'islkh', 'islkh',
                           'islkdox', 'islkdox', 'islkdox', 'islkdox', 'islkdox', 'islkdox', 'islkdox', 'islkdox', 'islkdox', 'islkdox',
                           'islkdox', 'islkdox', 'islkdox', 'islkdox', 'islkIDN', 'islkIDN', 'islkIDN', 'islkIDN', 'islkIDN', 'islkIDN',
                           'islkIDN', 'islkIDN', 'islkIDN', 'islkIDN', 'islkIDN', 'islkIDN', 'islkIDN', 'islkIDN', 'islkIFN', 'islkIFN',
                           'islkIFN', 'islkIFN', 'islkIFN', 'islkIFN', 'islkIFN', 'islkIFN', 'islkIFN', 'islkIFN', 'islkIFN', 'islkIFN',
                           'islkIFN', 'islkIFN')
cell.averages2$sample[cell.averages2$sample =='islkh'] <- 'A_islkH' # Plot islkH sample first

cell.averages2$names <- rownames(cell.averages2@meta.data)
cell.averages2$cluster <- as.character(cell.averages2$orig.ident)
for (x in seq(from = 0, to = 9)) {
  cell.averages2$cluster[cell.averages2$cluster == as.character(x)] <- paste('0',as.character(x), sep = '')
}
cell.averages2@meta.data <- arrange(cell.averages2@meta.data, cell.averages2$sample, cell.averages2$cluster)
cell.averages2.names <- cell.averages2$names
rownames(cell.averages2@meta.data) <- cell.averages2.names

# Average based on overall expression across all the different samples
islkh.averages2 <- subset(cell.averages2, subset = sample == 'A_islkH')
islkdox.averages2 <- subset(cell.averages2, subset = sample == 'islkdox')
islkIDN.averages2 <- subset(cell.averages2, subset = sample == 'islkIDN')
islkIFN.averages2 <- subset(cell.averages2, subset = sample == 'islkIFN')

# Get dataframe of scale.data average expression
islkh.average.data <- islkh.averages2 %>% GetAssayData(slot = 'scale.data') %>% as.data.frame()
islkdox.average.data <- islkdox.averages2 %>% GetAssayData(slot = 'scale.data') %>% as.data.frame()
islkIDN.average.data <- islkIDN.averages2 %>% GetAssayData(slot = 'scale.data') %>% as.data.frame()
islkIFN.average.data <- islkIFN.averages2 %>% GetAssayData(slot = 'scale.data') %>% as.data.frame()

colnames(islkh.average.data) <- islkh.averages2$orig.ident
islkh.average.data['12'] <- NA
colnames(islkdox.average.data) <- islkdox.averages2$orig.ident
colnames(islkIDN.average.data) <- islkIDN.averages2$orig.ident
colnames(islkIFN.average.data) <- islkIFN.averages2$orig.ident
islkh.average.data <- islkh.average.data[, clusters]
islkdox.average.data <- islkdox.average.data[, clusters]
islkIDN.average.data <- islkIDN.average.data[, clusters]
islkIFN.average.data <- islkIFN.average.data[, clusters]
islkh.average.data$gene <- rownames(islkh.average.data)
islkdox.average.data$gene <- rownames(islkdox.average.data)
islkIDN.average.data$gene <- rownames(islkIDN.average.data)
islkIFN.average.data$gene <- rownames(islkIFN.average.data)

# Standard deviation on overall expression across all the different samples
islkh.sd <- as.data.frame(matrix(nrow = nrow(cell.islkh), ncol = length(unique(cell.combined$seurat_clusters))))
rownames(islkh.sd) <- rownames(cell.combined)
islkdox.sd <- as.data.frame(matrix(nrow = nrow(cell.islkdox), ncol = length(unique(cell.combined$seurat_clusters))))
rownames(islkdox.sd) <- rownames(cell.combined)
islkIDN.sd <- as.data.frame(matrix(nrow = nrow(cell.islkIDN), ncol = length(unique(cell.combined$seurat_clusters))))
rownames(islkIDN.sd) <- rownames(cell.combined)
islkIFN.sd <- as.data.frame(matrix(nrow = nrow(cell.islkIFN), ncol = length(unique(cell.combined$seurat_clusters))))
rownames(islkIFN.sd) <- rownames(cell.combined)
clusters <- c('1', '10', '13', '2', '6', '4', '3', '0', '5', '8', '7', '12', '9', '11')

pb <- txtProgressBar(min = 0, max = 14, style = 3)
for (i in 0:13) {
  cat(paste0('Getting standard deviations for cluster ', as.character(i), '...\n'))
  if (i != 12) {
    islkh.temp <- subset(cell.islkh, seurat_clusters == as.character(i)) %>% GetAssayData(slot = 'scale.data') %>% as.data.frame()
    islkh.temp <- transform(islkh.temp, sd = apply(islkh.temp, 1, sd))
    islkh.sd[,i+1] <- islkh.temp$sd
  }
  colnames(islkh.sd)[i+1] <- as.character(i)
  islkdox.temp <- subset(cell.islkdox, seurat_clusters == as.character(i)) %>% GetAssayData(slot = 'scale.data') %>% as.data.frame()
  islkdox.temp <- transform(islkdox.temp, sd = apply(islkdox.temp, 1, sd))
  islkdox.sd[,i+1] <- islkdox.temp$sd
  colnames(islkdox.sd)[i+1] <- as.character(i)
  islkIDN.temp <- subset(cell.islkIDN, seurat_clusters == as.character(i)) %>% GetAssayData(slot = 'scale.data') %>% as.data.frame()
  islkIDN.temp <- transform(islkIDN.temp, sd = apply(islkIDN.temp, 1, sd))
  islkIDN.sd[,i+1] <- islkIDN.temp$sd
  colnames(islkIDN.sd)[i+1] <- as.character(i)
  islkIFN.temp <- subset(cell.islkIFN, seurat_clusters == as.character(i)) %>% GetAssayData(slot = 'scale.data') %>% as.data.frame()
  islkIFN.temp <- transform(islkIFN.temp, sd = apply(islkIFN.temp, 1, sd))
  islkIFN.sd[,i+1] <- islkIFN.temp$sd
  colnames(islkIFN.sd)[i+1] <- as.character(i)
  setTxtProgressBar(pb, pb$getVal()+1)
  cat('\n')
}
# Reorder columns
islkh.sd <- islkh.sd[,clusters]
islkdox.sd <- islkdox.sd[,clusters]
islkIDN.sd <- islkIDN.sd[,clusters]
islkIFN.sd <- islkIFN.sd[,clusters]
islkh.sd$gene <- rownames(islkh.sd)
islkdox.sd$gene <- rownames(islkdox.sd)
islkIDN.sd$gene <- rownames(islkIDN.sd)
islkIFN.sd$gene <- rownames(islkIFN.sd)
close(pb)

nfkb.genes <- c('IFNB1', 'NFKB1', 'NFKB2', 'RELA', 'RELB', 'CGAS', 'TMEM173', 'TBK1', 'IKBKE', 'IRF3', 'IRF7')

cell.combined.counts <- table(Idents(cell.combined), cell.combined$sample)
colnames(cell.combined.counts) <- c('iSLKH', 'iSLKdox', 'iSLKIDN', 'iSLKIFN')

# Write tables for average scale.data expression and standard deviation
write.csv(subset(islkh.average.data, gene %in% KSHVORFs.detected), '~/Documents/Genomics/reanalyze_with_insertion/cluster_avg_stdev/islkh_avg_ORFs.csv')
write.csv(subset(islkdox.average.data, gene %in% KSHVORFs.detected), '~/Documents/Genomics/reanalyze_with_insertion/cluster_avg_stdev/islkdox_avg_ORFs.csv')
write.csv(subset(islkIDN.average.data, gene %in% KSHVORFs.detected), '~/Documents/Genomics/reanalyze_with_insertion/cluster_avg_stdev/islkIDN_avg_ORFs.csv')
write.csv(subset(islkIFN.average.data, gene %in% KSHVORFs.detected), '~/Documents/Genomics/reanalyze_with_insertion/cluster_avg_stdev/islkIFN_avg_ORFs.csv')
write.csv(subset(islkh.average.data, gene %in% nfkb.genes), '~/Documents/Genomics/reanalyze_with_insertion/cluster_avg_stdev/islkh_avg_nfkb.csv')
write.csv(subset(islkdox.average.data, gene %in% nfkb.genes), '~/Documents/Genomics/reanalyze_with_insertion/cluster_avg_stdev/islkdox_avg_nfkb.csv')
write.csv(subset(islkIDN.average.data, gene %in% nfkb.genes), '~/Documents/Genomics/reanalyze_with_insertion/cluster_avg_stdev/islkIDN_avg_nfkb.csv')
write.csv(subset(islkIFN.average.data, gene %in% nfkb.genes), '~/Documents/Genomics/reanalyze_with_insertion/cluster_avg_stdev/islkIFN_avg_nfkb.csv')
write.csv(subset(islkh.sd, gene %in% KSHVORFs.detected), '~/Documents/Genomics/reanalyze_with_insertion/cluster_avg_stdev/islkh_sd_ORFs.csv')
write.csv(subset(islkdox.sd, gene %in% KSHVORFs.detected), '~/Documents/Genomics/reanalyze_with_insertion/cluster_avg_stdev/islkdox_sd_ORFs.csv')
write.csv(subset(islkIDN.sd, gene %in% KSHVORFs.detected), '~/Documents/Genomics/reanalyze_with_insertion/cluster_avg_stdev/islkIDN_sd_ORFs.csv')
write.csv(subset(islkIFN.sd, gene %in% KSHVORFs.detected), '~/Documents/Genomics/reanalyze_with_insertion/cluster_avg_stdev/islkIFN_sd_ORFs.csv')
write.csv(subset(islkh.sd, gene %in% nfkb.genes), '~/Documents/Genomics/reanalyze_with_insertion/cluster_avg_stdev/islkh_sd_nfkb.csv')
write.csv(subset(islkdox.sd, gene %in% nfkb.genes), '~/Documents/Genomics/reanalyze_with_insertion/cluster_avg_stdev/islkdox_sd_nfkb.csv')
write.csv(subset(islkIDN.sd, gene %in% nfkb.genes), '~/Documents/Genomics/reanalyze_with_insertion/cluster_avg_stdev/islkIDN_sd_nfkb.csv')
write.csv(subset(islkIFN.sd, gene %in% nfkb.genes), '~/Documents/Genomics/reanalyze_with_insertion/cluster_avg_stdev/islkIFN_sd_nfkb.csv')
write.csv(cell.combined.counts, '~/Documents/Genomics/reanalyze_with_insertion/cluster_avg_stdev/cluster_counts.csv')


umap.names <- c('','-by-sample','-split')
u0 <- DimPlot(cell.combined, reduction = 'umap', pt.size = 0.4, label = T, label.size = 7, cols = colors) + NoAxes()
u1 <- DimPlot(cell.combined, reduction = 'umap', pt.size = 0.4, group.by = 'sample') + NoAxes()
u2 <- DimPlot(cell.combined, reduction = 'umap', pt.size = 0.4, split.by = 'sample', label = T, label.size = 7, cols = colors) + NoAxes()
postscript(file = paste('~/Documents/Genomics/reanalyze_with_insertion/UMAP_Plots/cell-combined.eps', sep = ''), width = 8, height = 6)
print(u0)
dev.off()
postscript(file = paste('~/Documents/Genomics/reanalyze_with_insertion/UMAP_Plots/cell-combined-by-sample.eps', sep = ''), width = 8, height = 6)
print(u1)
dev.off()
postscript(file = paste('~/Documents/Genomics/reanalyze_with_insertion/UMAP_Plots/cell-combined-split.eps', sep = ''), width = 10, height = 4)
print(u2)
dev.off()

# Manually order clusters by percentage of genes mapping to KSHV genome
cell.combined$seurat_clusters <- factor(cell.combined$seurat_clusters, levels = c('10', '1', '13', '6', '4', '2', '3', '0', '9', '5', '11', '8', '7', '12'))
cell.islkh$seurat_clusters <- factor(cell.islkh$seurat_clusters, levels = c('10', '1', '13', '6', '4', '2', '3', '0', '9', '5', '11', '8', '7', '12'))
cell.islkdox$seurat_clusters <- factor(cell.islkdox$seurat_clusters, levels = c('10', '1', '13', '6', '4', '2', '3', '0', '9', '5', '11', '8', '7', '12'))
cell.islkIDN$seurat_clusters <- factor(cell.islkIDN$seurat_clusters, levels = c('10', '1', '13', '6', '4', '2', '3', '0', '9', '5', '11', '8', '7', '12'))
cell.islkIFN$seurat_clusters <- factor(cell.islkIFN$seurat_clusters, levels = c('10', '1', '13', '6', '4', '2', '3', '0', '9', '5', '11', '8', '7', '12'))

# Manually order clusters to follow pie graphs from figure 3F
cell.combined$seurat_clusters <- factor(cell.combined$seurat_clusters, levels = c('1', '10', '13', '2', '6', '4', '3', '0', '5', '8', '7', '12', '9', '11'))
cell.islkh$seurat_clusters <- factor(cell.islkh$seurat_clusters, levels = c('1', '10', '13', '2', '6', '4', '3', '0', '5', '8', '7', '12', '9', '11'))
cell.islkdox$seurat_clusters <- factor(cell.islkdox$seurat_clusters, levels = c('1', '10', '13', '2', '6', '4', '3', '0', '5', '8', '7', '12', '9', '11'))
cell.islkIDN$seurat_clusters <- factor(cell.islkIDN$seurat_clusters, levels = c('1', '10', '13', '2', '6', '4', '3', '0', '5', '8', '7', '12', '9', '11'))
cell.islkIFN$seurat_clusters <- factor(cell.islkIFN$seurat_clusters, levels = c('1', '10', '13', '2', '6', '4', '3', '0', '5', '8', '7', '12', '9', '11'))

levels(islkh.averages) <- c('1', '10', '13', '2', '6', '4', '3', '0', '5', '8', '7', '12', '9', '11')
levels(islkdox.averages) <- c('1', '10', '13', '2', '6', '4', '3', '0', '5', '8', '7', '12', '9', '11')
levels(islkIDN.averages) <- c('1', '10', '13', '2', '6', '4', '3', '0', '5', '8', '7', '12', '9', '11')
levels(islkIFN.averages) <- c('1', '10', '13', '2', '6', '4', '3', '0', '5', '8', '7', '12', '9', '11')
levels(islkh.averages2) <- c('1', '10', '13', '2', '6', '4', '3', '0', '5', '8', '7', '12', '9', '11')
levels(islkdox.averages2) <- c('1', '10', '13', '2', '6', '4', '3', '0', '5', '8', '7', '12', '9', '11')
levels(islkIDN.averages2) <- c('1', '10', '13', '2', '6', '4', '3', '0', '5', '8', '7', '12', '9', '11')
levels(islkIFN.averages2) <- c('1', '10', '13', '2', '6', '4', '3', '0', '5', '8', '7', '12', '9', '11')

p0 <- DoHeatmap(cell.combined, features = KSHVORFs.detected, size = 3.2, draw.lines = T, group.by = 'seurat_clusters', raster = F, group.colors = colors) + scale_fill_gradient2(low = 'blue', mid = 'white', high = 'red', midpoint = 0)
p1 <- DoHeatmap(cell.islkh, features = KSHVORFs.detected, size = 3.2, draw.lines = T, group.by = 'seurat_clusters', raster = F, group.colors = colors) + scale_fill_gradient2(low = 'blue', mid = 'white', high = 'red', midpoint = 0)
p2 <- DoHeatmap(cell.islkdox, features = KSHVORFs.detected, size = 3.2, draw.lines = T, group.by = 'seurat_clusters', raster = F, group.colors = colors) + scale_fill_gradient2(low = 'blue', mid = 'white', high = 'red', midpoint = 0)
p3 <- DoHeatmap(cell.islkIDN, features = KSHVORFs.detected, size = 3.2, draw.lines = T, group.by = 'seurat_clusters', raster = F, group.colors = colors) + scale_fill_gradient2(low = 'blue', mid = 'white', high = 'red', midpoint = 0)
p4 <- DoHeatmap(cell.islkIFN, features = KSHVORFs.detected, size = 3.2, draw.lines = T, group.by = 'seurat_clusters', raster = F, group.colors = colors) + scale_fill_gradient2(low = 'blue', mid = 'white', high = 'red', midpoint = 0)
p5 <- DoHeatmap(cell.combined, features = KSHVORFs.detected, size = 3.2, draw.lines = T, group.by = 'sample', raster = F, group.colors = colors) + scale_fill_gradient2(low = 'blue', mid = 'white', high = 'red', midpoint = 0)
p6 <- DoHeatmap(cell.islkIDN, features = c('IFNB1', 'RHEBL1', 'NFKB1', 'NFKB2', 'RELA', 'RELB', 'CGAS', 'TMEM173', 'TBK1', 'IKBKE', 'IRF3', 'IRF7'), size = 3.2, draw.lines = T, group.by = 'seurat_clusters', raster = F, group.colors = colors) + scale_fill_gradient2(low = 'blue', mid = 'white', high = 'red', midpoint = 0)
p6 <- DoHeatmap(cell.islkIFN, features = c('IFNB1', 'RHEBL1', 'NFKB1', 'NFKB2', 'RELA', 'RELB', 'CGAS', 'TMEM173', 'TBK1', 'IKBKE', 'IRF3', 'IRF7'), size = 3.2, draw.lines = T, group.by = 'seurat_clusters', raster = F, group.colors = colors) + scale_fill_gradient2(low = 'blue', mid = 'white', high = 'red', midpoint = 0)

# clusters of interest for Fig 4E
early.lytic.clusters <- c(0, 5, 7, 8, 9, 11)
IDN.cells <- colnames(subset(cell.islkIDN, subset = seurat_clusters %in% early.lytic.clusters))
IFN.cells <- colnames(subset(cell.islkIFN, subset = seurat_clusters %in% early.lytic.clusters))

p7 <- DoHeatmap(cell.islkIDN, cells = IDN.cells, features = c('IFNB1', KSHVORFs.detected), size = 3.2, draw.lines = F, group.by = 'seurat_clusters', raster = F, group.colors = colors) + scale_fill_gradient2(low = 'blue', mid = 'white', high = 'red', midpoint = 0)
p8 <- DoHeatmap(cell.islkIFN, cells = IFN.cells, features = c('IFNB1', KSHVORFs.detected), size = 3.2, draw.lines = F, group.by = 'seurat_clusters', raster = F, group.colors = colors) + scale_fill_gradient2(low = 'blue', mid = 'white', high = 'red', midpoint = 0)

postscript(file = paste('~/Documents/Genomics/reanalyze_with_insertion/KSHVORFS_Heatmaps/cell-combined.eps', sep = ''), width = 8, height = 6)
print(p0)
dev.off()
postscript(file = paste('~/Documents/Genomics/reanalyze_with_insertion/KSHVORFS_Heatmaps/iSLKH.eps', sep = ''), width = 8, height = 6)
print(p1)
dev.off()
postscript(file = paste('~/Documents/Genomics/reanalyze_with_insertion/KSHVORFS_Heatmaps/iSLKdox.eps', sep = ''), width = 8, height = 6)
print(p2)
dev.off()
postscript(file = paste('~/Documents/Genomics/reanalyze_with_insertion/KSHVORFS_Heatmaps/iSLKIDN.eps', sep = ''), width = 8, height = 6)
print(p3)
dev.off()
postscript(file = paste('~/Documents/Genomics/reanalyze_with_insertion/KSHVORFS_Heatmaps/iSLKIFN.eps', sep = ''), width = 8, height = 6)
print(p4)
dev.off()
postscript(file = paste('~/Documents/Genomics/reanalyze_with_insertion/KSHVORFS_Heatmaps/cell-combined-group-sample.eps', sep = ''), width = 8, height = 6)
print(p5)
dev.off()
postscript(file = paste('~/Documents/Genomics/reanalyze_with_insertion/IFN_differential_expression/NFKB_IFN_heatmap.eps', sep = ''), width = 8, height = 6)
print(p6)
dev.off()
postscript(file = paste('~/Documents/Genomics/reanalyze_with_insertion/IFN_differential_expression/islkIDN_early_lytic_ORFS.eps', sep = ''), width = 8, height = 6)
print(p7)
dev.off()
postscript(file = paste('~/Documents/Genomics/reanalyze_with_insertion/IFN_differential_expression/islkIFN_early_lytic_ORFS.eps', sep = ''), width = 8, height = 6)
print(p8)
dev.off()

p1.average <- DoHeatmap(islkh.averages, features = KSHVORFs.detected, size = 3.2, draw.lines = F) + scale_fill_gradient2(low = 'blue', mid = 'white', high = 'red', midpoint = 0)
p2.average <- DoHeatmap(islkdox.averages, features = KSHVORFs.detected, size = 3.2, draw.lines = F) + scale_fill_gradient2(low = 'blue', mid = 'white', high = 'red', midpoint = 0)
p3.average <- DoHeatmap(islkIDN.averages, features = KSHVORFs.detected, size = 3.2, draw.lines = F) + scale_fill_gradient2(low = 'blue', mid = 'white', high = 'red', midpoint = 0)
p4.average <- DoHeatmap(islkIFN.averages, features = KSHVORFs.detected, size = 3.2, draw.lines = F) + scale_fill_gradient2(low = 'blue', mid = 'white', high = 'red', midpoint = 0)
p5.average <- DoHeatmap(islkIDN.averages, features = c('IFNB1', 'RHEBL1', 'NFKB1', 'NFKB2', 'RELA', 'RELB', 'CGAS', 'TMEM173', 'TBK1', 'IKBKE', 'IRF3', 'IRF7'), size = 3.2, draw.lines = F, raster = F, group.colors = colors) + scale_fill_gradient2(low = 'blue', mid = 'white', high = 'red', midpoint = 0)
p6.average <- DoHeatmap(islkIFN.averages, features = c('IFNB1', 'RHEBL1', 'NFKB1', 'NFKB2', 'RELA', 'RELB', 'CGAS', 'TMEM173', 'TBK1', 'IKBKE', 'IRF3', 'IRF7'), size = 3.2, draw.lines = F, raster = F, group.colors = colors) + scale_fill_gradient2(low = 'blue', mid = 'white', high = 'red', midpoint = 0)
p7.average <- DoHeatmap(islkh.averages2, features = KSHVORFs.detected, size = 3.2, draw.lines = F) + scale_fill_gradient2(low = 'blue', mid = 'white', high = 'red', midpoint = 0)
p8.average <- DoHeatmap(islkdox.averages2, features = KSHVORFs.detected, size = 3.2, draw.lines = F) + scale_fill_gradient2(low = 'blue', mid = 'white', high = 'red', midpoint = 0)
p9.average <- DoHeatmap(islkIDN.averages2, features = KSHVORFs.detected, size = 3.2, draw.lines = F) + scale_fill_gradient2(low = 'blue', mid = 'white', high = 'red', midpoint = 0)
p10.average <- DoHeatmap(islkIFN.averages2, features = KSHVORFs.detected, size = 3.2, draw.lines = F) + scale_fill_gradient2(low = 'blue', mid = 'white', high = 'red', midpoint = 0)
p11.average <- DoHeatmap(islkh.averages2, features = nfkb.genes, size = 3.2, draw.lines = F) + scale_fill_gradient2(low = 'blue', mid = 'white', high = 'red', midpoint = 0)
p12.average <- DoHeatmap(islkdox.averages2, features = nfkb.genes, size = 3.2, draw.lines = F) + scale_fill_gradient2(low = 'blue', mid = 'white', high = 'red', midpoint = 0)
p13.average <- DoHeatmap(islkIDN.averages2, features = nfkb.genes, size = 3.2, draw.lines = F) + scale_fill_gradient2(low = 'blue', mid = 'white', high = 'red', midpoint = 0)
p14.average <- DoHeatmap(islkIFN.averages2, features = nfkb.genes, size = 3.2, draw.lines = F) + scale_fill_gradient2(low = 'blue', mid = 'white', high = 'red', midpoint = 0)



postscript(file = paste('~/Documents/Genomics/reanalyze_with_insertion/IFN_differential_expression/islkh_ORFs_average.eps', sep = ''), width = 8, height = 6)
print(p1.average)
dev.off()
postscript(file = paste('~/Documents/Genomics/reanalyze_with_insertion/IFN_differential_expression/islkdox_ORFs_average.eps', sep = ''), width = 8, height = 6)
print(p2.average)
dev.off()
postscript(file = paste('~/Documents/Genomics/reanalyze_with_insertion/IFN_differential_expression/islkIDN_ORFs_average.eps', sep = ''), width = 8, height = 6)
print(p3.average)
dev.off()
postscript(file = paste('~/Documents/Genomics/reanalyze_with_insertion/IFN_differential_expression/islkIFN_ORFs_average.eps', sep = ''), width = 8, height = 6)
print(p4.average)
dev.off()
postscript(file = paste('~/Documents/Genomics/reanalyze_with_insertion/IFN_differential_expression/islkIDN_NFKB_average.eps', sep = ''), width = 8, height = 6)
print(p5.average)
dev.off()
postscript(file = paste('~/Documents/Genomics/reanalyze_with_insertion/IFN_differential_expression/islkIFN_NFKB_average.eps', sep = ''), width = 8, height = 6)
print(p6.average)
dev.off()
postscript(file = paste('~/Documents/Genomics/reanalyze_with_insertion/IFN_differential_expression/islkh_ORFs_average2.eps', sep = ''), width = 8, height = 6)
print(p7.average)
dev.off()
postscript(file = paste('~/Documents/Genomics/reanalyze_with_insertion/IFN_differential_expression/islkdox_ORFs_average2.eps', sep = ''), width = 8, height = 6)
print(p8.average)
dev.off()
postscript(file = paste('~/Documents/Genomics/reanalyze_with_insertion/IFN_differential_expression/islkIDN_ORFs_average2.eps', sep = ''), width = 8, height = 6)
print(p9.average)
dev.off()
postscript(file = paste('~/Documents/Genomics/reanalyze_with_insertion/IFN_differential_expression/islkIFN_ORFs_average2.eps', sep = ''), width = 8, height = 6)
print(p10.average)
dev.off()
postscript(file = paste('~/Documents/Genomics/reanalyze_with_insertion/IFN_differential_expression/islkh_nfkb_average2.eps', sep = ''), width = 8, height = 6)
print(p11.average)
dev.off()
postscript(file = paste('~/Documents/Genomics/reanalyze_with_insertion/IFN_differential_expression/islkdox_nfkb_average2.eps', sep = ''), width = 8, height = 6)
print(p12.average)
dev.off()
postscript(file = paste('~/Documents/Genomics/reanalyze_with_insertion/IFN_differential_expression/islkIDN_nfkb_average2.eps', sep = ''), width = 8, height = 6)
print(p13.average)
dev.off()
postscript(file = paste('~/Documents/Genomics/reanalyze_with_insertion/IFN_differential_expression/islkIFN_nfkb_average2.eps', sep = ''), width = 8, height = 6)
print(p14.average)
dev.off()

p6 <- DoHeatmap(cell.averages2, features = KSHVORFs.detected, size = 3.2, draw.lines = T, group.by = 'sample', group.bar = T, cells = cell.averages2.names) + scale_fill_gradient2(low = 'blue', mid = 'white', high = 'red', midpoint = 0)

i0 <- DoHeatmap(cell.combined, features = IFNgenes.list, size = 3.2, draw.lines = T) + scale_fill_gradient2(low = 'blue', mid = 'white', high = 'red', midpoint = 0)
i1 <- DoHeatmap(cell.islkh, features = IFNgenes.list, size = 3.2, draw.lines = T) + scale_fill_gradient2(low = 'blue', mid = 'white', high = 'red', midpoint = 0)
i2 <- DoHeatmap(cell.islkdox, features = IFNgenes.list, size = 3.2, draw.lines = T) + scale_fill_gradient2(low = 'blue', mid = 'white', high = 'red', midpoint = 0)
i3 <- DoHeatmap(cell.islkIDN, features = IFNgenes.list, size = 3.2, draw.lines = T) + scale_fill_gradient2(low = 'blue', mid = 'white', high = 'red', midpoint = 0)
i4 <- DoHeatmap(cell.islkIFN, features = IFNgenes.list, size = 3.2, draw.lines = T) + scale_fill_gradient2(low = 'blue', mid = 'white', high = 'red', midpoint = 0)
i5 <- DoHeatmap(cell.combined, features = IFNgenes.list, size = 3.2, draw.lines = T, group.by = 'sample') + scale_fill_gradient2(low = 'blue', mid = 'white', high = 'red', midpoint = 0)
i6 <- DoHeatmap(cell.averages2, features = IFNgenes.list, size = 3.2, draw.lines = T, group.by = 'sample', group.bar = T, cells = cell.averages2.names) + scale_fill_gradient2(low = 'blue', mid = 'white', high = 'red', midpoint = 0)


# Apoptotic gene markers from iSLKIFN sample
apoptoticgenes <- c('IFNB1','CFLAR','DRAM1','FAS','JAK2','NFKBIA','PRELID1','RRAGC','SIX1','TNFAIP3','TRAF1','TNFRSF10B','TNFRSF9','TAX1BP1',
                    'XAF1','BIRC3','CASP4','CTSC','C19orf12','GADD45A','HIP1R','IRF1','MAP1S','NR3C1','PMAIP1','PLSCR1','PML','PPP1R15A',
                    'RIPK2','RFK','SQSTM1','TMEM173','TNFSF10','UBE2Z')

a0 <- DoHeatmap(cell.combined, features = apoptoticgenes, size = 3.2, draw.lines = T) + scale_fill_gradient2(low = 'blue', mid = 'white', high = 'red', midpoint = 0)
a1 <- DoHeatmap(cell.islkh, features = apoptoticgenes, size = 3.2, draw.lines = T) + scale_fill_gradient2(low = 'blue', mid = 'white', high = 'red', midpoint = 0)
a2 <- DoHeatmap(cell.islkdox, features = apoptoticgenes, size = 3.2, draw.lines = T) + scale_fill_gradient2(low = 'blue', mid = 'white', high = 'red', midpoint = 0)
a3 <- DoHeatmap(cell.islkIDN, features = apoptoticgenes, size = 3.2, draw.lines = T) + scale_fill_gradient2(low = 'blue', mid = 'white', high = 'red', midpoint = 0)
a4 <- DoHeatmap(cell.islkIFN, features = apoptoticgenes, size = 3.2, draw.lines = T) + scale_fill_gradient2(low = 'blue', mid = 'white', high = 'red', midpoint = 0)
a5 <- DoHeatmap(cell.combined, features = apoptoticgenes, size = 3.2, draw.lines = T, group.by = 'sample') + scale_fill_gradient2(low = 'blue', mid = 'white', high = 'red', midpoint = 0)
a6 <- DoHeatmap(cell.averages2, features = apoptoticgenes, size = 3.2, draw.lines = T, group.by = 'sample', group.bar = T, cells = cell.averages2.names) + scale_fill_gradient2(low = 'blue', mid = 'white', high = 'red', midpoint = 0)

cell.IFNcluster <- subset(cell.islkIDN, subset = (cell.islkIDN$seurat_clusters == 9))


#Graph % cells in each cluster across treatment conditions

cell.combined.percent <- 100*prop.table(table(Idents(cell.combined), cell.combined$sample), margin = 2) %>% as.data.frame.matrix()
cell.combined.counts <- table(Idents(cell.combined), cell.combined$sample)
cluster.percentage.plot <- ggplot(cell.combined.percent, aes(fill=..Var1, x=..Var2, y=..Freq)) + geom_bar(position="stack", stat="identity") + theme_bw() +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
postscript(file = paste('~/Documents/Genomics/reanalyze_with_insertion/Cluster_Counts/cluster_percentage.eps', sep = ''), width = 6, height = 4)
print(cluster.percentage.plot)
dev.off()

