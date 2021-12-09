library(dplyr)
library(Seurat)
library(cowplot)
library(ggpubr)

# Loading data into monocle object
# print('Loading 10X data...')
# islkhm.data <- load_cellranger_data('~/Documents/Genomics/iSLKH_iSLK219_count/')
# islkhm.data$sampleID <- 'A'
# islkdoxm.data <- load_cellranger_data('~/Documents/Genomics/iSLK219_dox_count/')
# islkdoxm.data$sampleID <- 'B'
# islk_IDNm.data <- load_cellranger_data('~/Documents/Genomics/iSLK219_IDN_count/')
# islk_IDNm.data$sampleID <- 'C'
# islk_IFNm.data <- load_cellranger_data('~/Documents/Genomics/iSLK219_IDN_IFN_count/')
# islk_IFNm.data$sampleID <- 'D'
# cellm.combined <- combine_cds(list(islkhm.data, islkdoxm.data, islk_IDNm.data, islk_IFNm.data))

s.genes <- cc.genes.updated.2019$s.genes
g2m.genes <- cc.genes.updated.2019$g2m.genes

cell.combined.cc <- CellCycleScoring(cell.combined, s.features = s.genes, g2m.features = g2m.genes, set.ident = T)
#rm(cell.combined)
cell.islkh.cc <- subset(cell.combined.cc, orig.ident == 'islkh')
cell.islkdox.cc <- subset(cell.combined.cc, orig.ident == 'islkdox')
cell.islkIDN.cc <- subset(cell.combined.cc, orig.ident == 'islkIDN')
cell.islkIFN.cc <- subset(cell.combined.cc, orig.ident == 'islkIFN')


gene.list <- c('GFP', 'DSRED', 'ORF73','ORFK4.1','ORF34','ORF42','ISG15','GAPDH','IFNB1')
phase.list <- c('G1','G2M','S')
#phase.list <- rep(c('G1','G2M','S'),c(rep(length(gene.list),3)))
#gene.list <- rep(gene.list, 3)
for (y in phase.list) {
  lapply(gene.list, function(x) {
    cc1 <- FeaturePlot(subset(cell.islkh.cc, Phase==y), features = c(x), max.cutoff = 'q99', pt.size = 0.4, order = T, col = c('lightblue','red')) + NoAxes()
    cc2 <- FeaturePlot(subset(cell.islkdox.cc, Phase==y), features = c(x), max.cutoff = 'q99', pt.size = 0.4, order = T, col = c('lightblue','red')) + NoAxes()
    cc3 <- FeaturePlot(subset(cell.islkIDN.cc, Phase==y), features = c(x), max.cutoff = 'q99', pt.size = 0.4, order = T, col = c('lightblue','red')) + NoAxes()
    cc4 <- FeaturePlot(subset(cell.islkIFN.cc, Phase==y), features = c(x), max.cutoff = 'q99', pt.size = 0.4, order = T, col = c('lightblue','red')) + NoAxes()
    #cc5 <- DimPlot(cell.combined, split.by = 'sample', label = T, repel = T, pt.size = 0.8) + NoAxes() + NoLegend()
    postscript(file = paste('~/Documents/Genomics/reanalyze_with_insertion/cellcycle/integratedscale', x, y, '.eps', sep = ''), width = 8, height = 3)
    print(plot_grid(cc1, cc2, cc3, cc4, nrow = 1))
    dev.off()
  })
}
cc5 <- DimPlot(cell.combined.cc, split.by = 'sample', label = F, repel = T, pt.size = 0.8) + NoAxes()

for (y in 0:13) {
  cc6 <- DimPlot(subset(cell.combined.cc, seurat_clusters == y), split.by = 'sample', label = F, repel = T, pt.size = 0.8) + NoAxes()
  postscript(file = paste('~/Documents/Genomics/reanalyze_with_insertion/cellcycle/cluster', y, '.eps', sep = ''), width = 8, height = 3)
  print(cc6)
  dev.off()
}
print(phase.list)
