library(dplyr)
library(Seurat)
library(cowplot)
library(ggpubr)

#gene.list <- c('GFP','DSRED','PAC','ZBP1','ORF73','ORFK4.1','ORF34','ORF42','PAN','ISG15','GAPDH','IFNB1','DUSP11','ACTB','G6PD','PGK1',
#               'POLR2A','ATF1','EIF2A','TMEM173','CGAS','IRF3','MAVS','TICAM1','DDX58','RHEBL1','PPM1G','NFKB1','NFKB2','RELA','RELB')
gene.list <- c('CCL5', 'TRAF1', 'PMAIP1', 'TNFRSF9', 'SPRY2', 'BBC3', 'NEURL3', 'FZD4', 'ZFP36L2', 'IL27', 'TLR3', 'TLR4')
cell.islkh <- subset(cell.combined, orig.ident == 'islkh')
cell.islkdox <- subset(cell.combined, orig.ident == 'islkdox')
cell.islkIDN <- subset(cell.combined, orig.ident == 'islkIDN')
cell.islkIFN <- subset(cell.combined, orig.ident == 'islkIFN')
# cell.islkh <- subset(islk, orig.ident == 'islkh')
# cell.islkdox <- subset(islk, orig.ident == 'islkdox')
# cell.islkIDN <- subset(islk, orig.ident == 'islkIDN')
# cell.islkIFN <- subset(islk, orig.ident == 'islkIFN')

#IFNb positive cells in integrated analysis
# findGenePercentages <- function(genes, celldata) {
#   populations <- unique(celldata$orig.ident)
#   subset.list <-   append(lapply(populations, FUN = function(x) { subset(celldata, orig.ident == x)}))
#   genes.data <- do.calllapply(genes, FUN = function(x) {
#   genes.sample <- append(lapply(subset.list, FUN = function(x) {x$orig.ident[[1]]}), 'total')
#   genes.count <- append(lapply(subset.list, FUN = function(x) {sum(GetAssayData(x, slot = 'data')['IFNB1',]>0)}),
#                         sum(GetAssayData(cell.combined, slot = 'data')['IFNB1',]>0))
#   genes.percent <- append(lapply(subset.list, FUN = function(x) {sum(GetAssayData(x, slot = 'data')['IFNB1',]>0)/nrow(x@meta.data)*100}),
#                          sum(GetAssayData(cell.combined, slot = 'data')['IFNB1',]>0)/nrow(cell.combined@meta.data)*100)
#   
#   })
# }
# subset.list <- c(islkh, islkdox, islkIDN, islkIFN)
# IFNb.sample <- append(lapply(subset.list, FUN = function(x) {x$orig.ident[[1]]}), 'total')
# IFNb.count <- append(lapply(subset.list, FUN = function(x) {sum(GetAssayData(x, slot = 'data')['IFNB1',]>0)}),
#               sum(GetAssayData(islk, slot = 'data')['IFNB1',]>0))
# IFNb.percent <- append(lapply(subset.list, FUN = function(x) {sum(GetAssayData(x, slot = 'data')['IFNB1',]>0)/nrow(x@meta.data)*100}),
#                 sum(GetAssayData(islk, slot = 'data')['IFNB1',]>0)/nrow(islk@meta.data)*100)
# IFNB.data.merged <- do.call(rbind.data.frame, Map('c', IFNb.sample, IFNb.count, IFNb.percent))
# colnames(IFNB.data.merged) <- c('sample', 'number IFNB1+ cells', 'percentage IFNB1+')
# write.csv(IFNB.data.merged, file = '~/Documents/Genomics/MergedAnalysis/IFNB1cells.csv', quote = F, row.names = F)

lapply(gene.list, FUN = function(x) {
  #FeaturePlot(cell.combined, features = c(x), split.by = 'sample', min.cutoff = 'q50', max.cutoff = 'q90', order=T, col = c('gray', 'red'))+NoAxes()
  #print(plot_grid(FeaturePlot(cell.combined, features = c(x), split.by = 'sample', min.cutoff = 'q15', max.cutoff = 'q90', order=T, col = c('gray','red'))+NoAxes(),
  #          DimPlot(cell.combined, group.by = 'sample')+NoLegend()+NoAxes(), nrow = 1, rel_widths = c(4.35,1)))
  # p1 <- ggpubr::annotate_figure(
  #   p = FeaturePlot(cell.islkh, features = c(x), max.cutoff = 'q99', pt.size = 0.9, order = T, col = c('light blue','red'))+NoAxes()+NoLegend(),
  #   top = ggpubr::text_grob(label = 'iSLKH', face = 'bold')
  # )
  # p2 <- ggpubr::annotate_figure(
  #   p = FeaturePlot(cell.islkdox, features = c(x), max.cutoff = 'q99', pt.size = 0.9, order = T, col = c('light blue','red'))+NoAxes()+NoLegend(),
  #   top = ggpubr::text_grob(label = 'iSLKdox', face = 'bold')
  # )
  # p3 <- ggpubr::annotate_figure(
  #   p = FeaturePlot(cell.islkIDN, features = c(x), max.cutoff = 'q99', pt.size = 0.9, order = T, col = c('light blue','red'))+NoAxes()+NoLegend(),
  #   top = ggpubr::text_grob(label = 'iSLKIDN', face = 'bold')
  # )
  # p4 <- ggpubr::annotate_figure(
  #   p = FeaturePlot(cell.islkIFN, features = c(x), max.cutoff = 'q99', pt.size = 0.9, order = T, col = c('light blue','red'))+NoAxes()+NoLegend(),
  #   top = ggpubr::text_grob(label = 'iSLKIFN', face = 'bold')
  # )
  p1 <- FeaturePlot(cell.islkh, features = c(x), max.cutoff = 'q99', pt.size = 0.4, order = T, col = c('lightblue','red')) + NoAxes()
  p2 <- FeaturePlot(cell.islkdox, features = c(x), max.cutoff = 'q99', pt.size = 0.4, order = T, col = c('lightblue','red')) + NoAxes()
  p3 <- FeaturePlot(cell.islkIDN, features = c(x), max.cutoff = 'q99', pt.size = 0.4, order = T, col = c('lightblue','red')) + NoAxes()
  p4 <- FeaturePlot(cell.islkIFN, features = c(x), max.cutoff = 'q99', pt.size = 0.4, order = T, col = c('lightblue','red')) + NoAxes()
  # p5 <- DimPlot(cell.combined, split.by = 'sample', label = T, repel = T, pt.size = 0.8) + NoAxes() + NoLegend()
  postscript(file = paste('~/Documents/Genomics/reanalyze_with_insertion/feature_plots/integratedscale', x,'.eps', sep = ''), width = 8, height = 3)
  print(plot_grid(p1, p2, p3, p4, nrow = 1))
  dev.off()
})


# p1 <- DimPlot(cell.islkh, label = T, label.size = 8, repel = T, pt.size = 0.8) + NoAxes() + NoLegend()
# p2 <- DimPlot(cell.islkdox, label = T, label.size = 8, repel = T, pt.size = 0.8) + NoAxes() + NoLegend()
# p3 <- DimPlot(cell.islkIDN, label = T, label.size = 8, repel = T, pt.size = 0.8) + NoAxes() + NoLegend()
# p4 <- DimPlot(cell.islkIFN, label = T, label.size = 8, repel = T, pt.size = 0.8) + NoAxes() + NoLegend()
p1 <- ggpubr::annotate_figure(
  p = DimPlot(cell.islkh, label = T, label.size = 4, repel = T, pt.size = 0.2) + NoAxes() + NoLegend(),
  top = ggpubr::text_grob(label = 'iSLKH', face = 'bold', size = 24)
)
p2 <- ggpubr::annotate_figure(
  p = DimPlot(cell.islkdox, label = T, label.size = 4, repel = T, pt.size = 0.2) + NoAxes() + NoLegend(),
  top = ggpubr::text_grob(label = 'iSLKdox', face = 'bold', size = 24)
)
p3 <- ggpubr::annotate_figure(
  p = DimPlot(cell.islkIDN, label = T, label.size = 4, repel = T, pt.size = 0.2) + NoAxes() + NoLegend(),
  top = ggpubr::text_grob(label = 'iSLKIDN', face = 'bold', size = 24)
)
p4 <- ggpubr::annotate_figure(
  p = DimPlot(cell.islkIFN, label = T, label.size = 4, repel = T, pt.size = 0.2) + NoAxes() + NoLegend(),
  top = ggpubr::text_grob(label = 'iSLKIFN', fac = 'bold', size = 24)
)
postscript(file = '~/Documents/Genomics/reanalyze_with_insertion/feature_plots/integrated_clusters.eps', width = 6, height = 3)
print(plot_grid(p1, p2, p3, p4, nrow = 1))
dev.off()


cell.islkh.DUSP11 <- FetchData(cell.islkh, 'DUSP11')
cell.islkdox.DUSP11 <- FetchData(cell.islkdox, 'DUSP11')
cell.islkIDN.DUSP11 <- FetchData(cell.islkIDN, 'DUSP11')
cell.islkIFN.DUSP11 <- FetchData(cell.islkIFN, 'DUSP11')
mean_DUSP11 <- c('Uninfected/Latent','Lytic+vehicle','Lytic+IDN-6556','Lytic+IDN-6556+anti-IFNs')
mean_DUSP11 <- rbind(mean_DUSP11, c(mean(cell.islkh.DUSP11[,1]),mean(cell.islkdox.DUSP11[,1]),mean(cell.islkIDN.DUSP11[,1]),mean(cell.islkIFN.DUSP11[,1])))
write.table(mean_DUSP11, file = '~/Documents/Genomics/reanalyze_with_insertion/feature_plots/DUSP11.txt', row.names = F, col.names = F, quote = F,sep = ',')

# RHEBL1 and IFNB1 correlation plots
# p1 <- FeaturePlot(cell.islkh, features = c('RHEBL1', 'IFNB1'), max.cutoff = 'q99', pt.size = 1, blend = T, order = T, col = c('lightblue','blue', 'green'))
# p2 <- FeaturePlot(cell.islkdox, features = c('RHEBL1', 'IFNB1'), max.cutoff = 'q99', pt.size = 1, blend = T, order = T, col = c('lightblue','blue', 'green'))
p3 <- FeaturePlot(cell.islkIDN, features = c('RHEBL1', 'IFNB1'), max.cutoff = 'q99', pt.size = 1, blend = T, blend.threshold = 0.5, order = T, col = c('lightblue','red', 'blue'))
# p4 <- FeaturePlot(cell.islkIFN, features = c('RHEBL1', 'IFNB1'), max.cutoff = 'q99', pt.size = 1, blend = T, order = T, col = c('lightblue','blue', 'yellow'))
# p5 <- FeaturePlot(subset(cell.islkIDN, subset = ((seurat_clusters == 9 | seurat_clusters == 11))), features = c('RHEBL1', 'IFNB1'),
#                   max.cutoff = 'q99', pt.size = 1.5, blend = T, order = T, col = c('lightblue','blue', 'yellow'))
postscript(file = paste('~/Documents/Genomics/reanalyze_with_insertion/IFN_differential_expression/rhebl1_ifnb.eps', sep = ''), width = 8, height = 3)
print(p3)
dev.off()

