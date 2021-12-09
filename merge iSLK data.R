library(dplyr)
library(Seurat)
library(cowplot)
library(ggplot2)

# Read 10X data
print('Reading 10X data...')
islkh.data <- Read10X(data.dir = '~/Documents/Genomics/iSLKH_iSLK219_count/outs/filtered_feature_bc_matrix/')
islkdox.data <- Read10X(data.dir = '~/Documents/Genomics/iSLK219_dox_count/outs/filtered_feature_bc_matrix/')
islk_IDN.data <- Read10X(data.dir = '~/Documents/Genomics/iSLK219_IDN_count/outs/filtered_feature_bc_matrix/')
islk_IFN.data <- Read10X(data.dir = '~/Documents/Genomics/iSLK219_IDN_IFN_count/outs/filtered_feature_bc_matrix/')

print('Creating seurat objects...')
islkh <- CreateSeuratObject(counts = islkh.data, project = 'islkh')
islkdox <- CreateSeuratObject(counts = islkdox.data, project = 'islkdox')
islkIDN <- CreateSeuratObject(counts = islk_IDN.data, project = 'islkIDN')
islkIFN <- CreateSeuratObject(counts = islk_IFN.data, project = 'islkIFN')

print('merging seurat objects together...')
islk <- merge(islkh, y = c(islkdox, islkIDN, islkIFN),
              add.cell.ids = c('iSLKH_iSLK219','iSLK219+dox','iSLK219+dox+IDN','iSLK219+dox+IDN+IFN'),
              project = 'expt_283')

# Add sample ID to metadata
islk$sample <- 'temp'
islk$sample[islk$orig.ident=='islkh'] <- 'A'
islk$sample[islk$orig.ident=='islkdox'] <- 'B'
islk$sample[islk$orig.ident=='islkIDN'] <- 'C'
islk$sample[islk$orig.ident=='islkIFN'] <- 'D'

# Normalize data and perfrom dimensional reduction analysis
print('Normalizing data...')
islk[['percent.mt']] <- PercentageFeatureSet(islk, pattern = '^MT-')
islk <- subset(islk, subset = percent.mt < 60)
islk <- NormalizeData(islk)
islk <- FindVariableFeatures(islk, selection.method = 'vst', nfeatures = 2000)

print('Performing dimensional reduction and finding clusters...')
all.genes <- rownames(islk)
islk <- ScaleData(islk, features = all.genes)
islk <- RunPCA(islk, features = VariableFeatures(object = islk))
islk <- FindNeighbors(islk, dims = 1:30)
islk <- FindClusters(islk, resolution = 0.5)
islk <- RunUMAP(islk, dims = 1:30)
# Visualization
p1 <- DimPlot(islk, reduction = "umap", label = TRUE, repel = T) + NoLegend()
p2 <- DimPlot(islk, reduction = "umap", group.by = 'sample') + NoLegend()
p3 <- plot_grid(p1, p2)

print('Finding markers for each cluster...')
islk.markers <- FindAllMarkers(islk, only.pos = T, min.pct = 0.25, logfc.threshold = 0.25)
islk.top5 <- islk.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)
islk.heatmap <- DoHeatmap(islk, features = islk.top5$gene)

print('Saving to RDS file...')
saveRDS(islk, file = '~/Documents/Genomics/islk-merged.rds')
