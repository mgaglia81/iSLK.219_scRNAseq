library(dplyr)
library(Seurat)
library(cowplot)
library(ggpubr)
library(ggplot2)
library(tidyr)
library(tibble)
library(scales)

#Make sure cell.combined is loaded from RDS file
#cell.combined <- readRDS('~/Documents/Genomics/reanalyze_with_insertion/cell-combined.rds')
#cell.islkh <- subset(cell.combined, orig.ident == 'islkh')
#cell.islkdox <- subset(cell.combined, orig.ident == 'islkdox')
#cell.islkIDN <- subset(cell.combined, orig.ident == 'islkIDN')
#cell.islkIFN <- subset(cell.combined, orig.ident == 'islkIFN')

# Manually order clusters to follow pie graphs from figure 3F
cell.combined$seurat_clusters <- factor(cell.combined$seurat_clusters, levels = c('1', '10', '13', '2', '6', '4', '3', '0', '5', '8', '7', '12', '9', '11'))
cell.islkh$seurat_clusters <- factor(cell.islkh$seurat_clusters, levels = c('1', '10', '13', '2', '6', '4', '3', '0', '5', '8', '7', '12', '9', '11'))
cell.islkdox$seurat_clusters <- factor(cell.islkdox$seurat_clusters, levels = c('1', '10', '13', '2', '6', '4', '3', '0', '5', '8', '7', '12', '9', '11'))
cell.islkIDN$seurat_clusters <- factor(cell.islkIDN$seurat_clusters, levels = c('1', '10', '13', '2', '6', '4', '3', '0', '5', '8', '7', '12', '9', '11'))
cell.islkIFN$Seurat_clusters <- factor(cell.islkIFN$seurat_clusters, levels = c('1', '10', '13', '2', '6', '4', '3', '0', '5', '8', '7', '12', '9', '11'))


KSHVORFs.detected <- c('GFP', 'DSRED', 'PAC', 'vIL6', 'PAN', 'ORFK4.1', 'ORFK12', 'ORF73', 'ORF50', 'ORF57', 'ORFK1', 'ORF6', '1.4kb', 'ORFK5', 'ORFK6', 'ORF16',
                       'ORF34', 'ORF39', 'ORF8', 'ORF17', 'ORF21', 'ORF25', 'ORF28', 'ORF42', 'ORF53', 'ORF55', 'ORF62', 'ORF65', 'ORF69', 'ORF75')

KSHV.latent <- c('GFP', 'PAC', 'ORF73', 'ORFK12')
KSHV.early <- c('vIL6', 'ORF50','ORF57','ORFK1','ORFK4.1','ORF6','ORFK5', 'ORFK6', 'ORF16', 'ORF34', 'ORF39', 'ORF8', 'ORF53', 'ORF55', 'ORF69', '1.4kb', 'DSRED', 'PAN')
#KSHV.intermediate <- c('ORFK4.1', 'ORF42', 'ORF75')
KSHV.late <- c('ORF42', 'ORF17', 'ORF21', 'ORF25', 'ORF28', 'ORF62', 'ORF65', 'ORF75')

colors <- c('1' = '#D6D6D5', '10' = '#2935E5', '13' = '#307BFE', '2' = '#76D5FE', '6' = '#33E866', '4' = '#2EAD1F',
            '3' = '#87E12C', '0' = '#FCFA43', '5' = '#EEC85F', '8' = '#FE9705', '7' = '#CA5915', '12' = '#FD2513',
            '9' = '#B04CF2', '11' = '#8336FF')

KSHV.reporters <- c('GFP', 'DSRED', 'PAC')
KSHV.latent2 <- c('ORFK12', 'ORF73')
KSHV.early2 <- c('1.4kb', 'PAN', 'vIL6', 'ORFK1', 'ORFK4.1', 'ORFK5', 'ORFK6', 'ORF6', 'ORF8', 'ORF16', 'ORF34', 'ORF39', 'ORF50', 'ORF53', 'ORF55', 'ORF57', 'ORF69')
KSHV.late2 <- c('ORF17', 'ORF21', 'ORF25', 'ORF28', 'ORF42', 'ORF62', 'ORF65', 'ORF75')
KSHVORFs.detected <- c(KSHV.reporters, KSHV.latent2, KSHV.early2, KSHV.late2)

cell.combined[['percent.KSHV']] <- PercentageFeatureSet(object = cell.combined, features = KSHVORFs.detected)
cell.combined[['percent.latent']] <- PercentageFeatureSet(object = cell.combined, features = KSHV.latent)
cell.combined[['percent.early']] <- PercentageFeatureSet(object = cell.combined, features = KSHV.early)
#cell.combined[['percent.intermediate']] <- PercentageFeatureSet(object = cell.combined, features = KSHV.intermediate)
cell.combined[['percent.late']] <- PercentageFeatureSet(object = cell.combined, features = KSHV.late)
cell.islkh <- subset(cell.combined, orig.ident == 'islkh')
cell.islkdox <- subset(cell.combined, orig.ident == 'islkdox')
cell.islkIDN <- subset(cell.combined, orig.ident == 'islkIDN')
cell.islkIFN <- subset(cell.combined, orig.ident == 'islkIFN')
num.clusters <- max(as.numeric(levels(cell.combined)))

viral.genes <- VlnPlot(cell.combined, features = c('percent.KSHV'), pt.size = 0, group.by = 'sample') + NoLegend() +
  stat_summary(fun.y = median, geom = 'point', size = 8, color = 'black', shape = 95)

viral.genes + scale_y_log10(limits = c(0.001,100), breaks = c(0.001, 0.01, 0.1, 1, 10, 100), labels = scales::comma)
p1 <- FeaturePlot(cell.islkh, features = 'percent.KSHV', max.cutoff = 110, pt.size = 0.4, order = T) + NoAxes() + scale_color_gradientn( colors = c('lightblue', 'red'),  limits = c(0, 100))
p2 <- FeaturePlot(cell.islkdox, features = 'percent.KSHV', max.cutoff = 110, pt.size = 0.4, order = T, col = c('lightblue','red')) + NoAxes() + scale_color_gradientn( colors = c('lightblue', 'red'),  limits = c(0, 100))
p3 <- FeaturePlot(cell.islkIDN, features = 'percent.KSHV', max.cutoff = 110, pt.size = 0.4, order = T, col = c('lightblue','red')) + NoAxes() + scale_color_gradientn( colors = c('lightblue', 'red'),  limits = c(0, 100))
p4 <- FeaturePlot(cell.islkIFN, features = 'percent.KSHV', max.cutoff = 110, pt.size = 0.4, order = T, col = c('lightblue','red')) + NoAxes() + scale_color_gradientn( colors = c('lightblue', 'red'),  limits = c(0, 100))
postscript(file = paste('~/Documents/Genomics/reanalyze_with_insertion/viral_genes/KSHVgenes100.eps', sep = ''), width = 4, height = 3)
print(plot_grid(p1, p2, p3, p4, nrow = 1))
dev.off()

for (i in 0:num.clusters) {
  print(paste0('Writing txt for cluster ', as.character(i)))
  islkh <- c(subset(cell.islkh, seurat_clusters == i)$percent.latent,
             subset(cell.islkh, seurat_clusters == i)$percent.early,
             #subset(cell.islkh, seurat_clusters == i)$percent.intermediate,
             subset(cell.islkh, seurat_clusters == i)$percent.late)
  islkdox <- c(subset(cell.islkdox, seurat_clusters == i)$percent.latent,
             subset(cell.islkdox, seurat_clusters == i)$percent.early,
             #subset(cell.islkdox, seurat_clusters == i)$percent.intermediate,
             subset(cell.islkdox, seurat_clusters == i)$percent.late)
  islkIDN <- c(subset(cell.islkIDN, seurat_clusters == i)$percent.latent,
             subset(cell.islkIDN, seurat_clusters == i)$percent.early,
             #subset(cell.islkIDN, seurat_clusters == i)$percent.intermediate,
             subset(cell.islkIDN, seurat_clusters == i)$percent.late)
  islkIFN <- c(subset(cell.islkIFN, seurat_clusters == i)$percent.latent,
             subset(cell.islkIFN, seurat_clusters == i)$percent.early,
             #subset(cell.islkIFN, seurat_clusters == i)$percent.intermediate,
             subset(cell.islkIFN, seurat_clusters == i)$percent.late)
  islkh.mean <- c(mean(subset(cell.islkh, seurat_clusters == i)$percent.latent),
                  mean(subset(cell.islkh, seurat_clusters == i)$percent.early),
                  #mean(subset(cell.islkh, seurat_clusters == i)$percent.intermediate),
                  mean(subset(cell.islkh, seurat_clusters == i)$percent.late))
  islkdox.mean <- c(mean(subset(cell.islkdox, seurat_clusters == i)$percent.latent),
                  mean(subset(cell.islkdox, seurat_clusters == i)$percent.early),
                  #mean(subset(cell.islkdox, seurat_clusters == i)$percent.intermediate),
                  mean(subset(cell.islkdox, seurat_clusters == i)$percent.late))
  islkIDN.mean <- c(mean(subset(cell.islkIDN, seurat_clusters == i)$percent.latent),
                  mean(subset(cell.islkIDN, seurat_clusters == i)$percent.early),
                  #mean(subset(cell.islkIDN, seurat_clusters == i)$percent.intermediate),
                  mean(subset(cell.islkIDN, seurat_clusters == i)$percent.late))
  islkIFN.mean <- c(mean(subset(cell.islkIFN, seurat_clusters == i)$percent.latent),
                  mean(subset(cell.islkIFN, seurat_clusters == i)$percent.early),
                  #mean(subset(cell.islkIFN, seurat_clusters == i)$percent.intermediate),
                  mean(subset(cell.islkIFN, seurat_clusters == i)$percent.late))
  islkh <- matrix(islkh, ncol = 4)
  colnames(islkh) <- c('latent','early','intermediate','late')
  islkdox <- matrix(islkdox, ncol = 4)
  colnames(islkdox) <- c('latent','early','intermediate','late')
  islkIDN <- matrix(islkIDN, ncol = 4)
  colnames(islkIDN) <- c('latent','early','intermediate','late')
  islkIFN <- matrix(islkIFN, ncol = 4)
  colnames(islkIFN) <- c('latent','early','intermediate','late')
  islkh.mean <- matrix(islkh.mean, ncol = 4)
  colnames(islkh.mean) <- c('latent','early','intermediate','late')
  islkdox.mean <- matrix(islkdox.mean, ncol = 4)
  colnames(islkdox.mean) <- c('latent','early','intermediate','late')
  islkIDN.mean <- matrix(islkIDN.mean, ncol = 4)
  colnames(islkIDN.mean) <- c('latent','early','intermediate','late')
  islkIFN.mean <- matrix(islkIFN.mean, ncol = 4)
  colnames(islkIFN.mean) <- c('latent','early','intermediate','late')
  
  write.table(islkh,
              file = paste('~/Documents/Genomics/reanalyze_with_insertion/viral_genes/viral_stage/islkh_', if(i < 10) { '0' }, as.character(i),'.txt', sep = ''),
              row.names = F, col.names = T, quote = F)
  write.table(islkh.mean,
              file = paste('~/Documents/Genomics/reanalyze_with_insertion/viral_genes/viral_stage/islkh_mean_', if(i < 10) { '0' }, as.character(i),'.txt', sep = ''),
              row.names = F, col.names = T, quote = F)
  write.table(islkdox,
              file = paste('~/Documents/Genomics/reanalyze_with_insertion/viral_genes/viral_stage/islkdox_', if(i < 10) { '0' }, as.character(i),'.txt', sep = ''),
              row.names = F, col.names = T, quote = F)
  write.table(islkdox.mean,
              file = paste('~/Documents/Genomics/reanalyze_with_insertion/viral_genes/viral_stage/islkdox_mean_', if(i < 10) { '0' }, as.character(i),'.txt', sep = ''),
              row.names = F, col.names = T, quote = F)
  write.table(islkIDN,
              file = paste('~/Documents/Genomics/reanalyze_with_insertion/viral_genes/viral_stage/islkIDN_', if(i < 10) { '0' }, as.character(i),'.txt', sep = ''),
              row.names = F, col.names = T, quote = F)
  write.table(islkIDN.mean,
              file = paste('~/Documents/Genomics/reanalyze_with_insertion/viral_genes/viral_stage/islkIDN_mean_', if(i < 10) { '0' }, as.character(i),'.txt', sep = ''),
              row.names = F, col.names = T, quote = F)
  write.table(islkIFN,
              file = paste('~/Documents/Genomics/reanalyze_with_insertion/viral_genes/viral_stage/islkIFN_', if(i < 10) { '0' }, as.character(i),'.txt', sep = ''),
              row.names = F, col.names = T, quote = F)
  write.table(islkIFN.mean,
              file = paste('~/Documents/Genomics/reanalyze_with_insertion/viral_genes/viral_stage/islkIFN_mean_', if(i < 10) { '0' }, as.character(i),'.txt', sep = ''),
              row.names = F, col.names = T, quote = F)
}


samples.list <- c('islkh', 'islkdox', 'islkIDN', 'islkIFN')
lapply(samples.list, FUN = function(obj) {
  print(paste0('Working on ', obj, ' library'))
  clusters <- rep(0:num.clusters, each = 3)
  stage <- rep(c('latent', 'early', 'late'), num.clusters+1)
  value <- c()
  mean.KSHV <- c()
  deviation.KSHV <- c()
  if (obj == 'islkh') {
    cell.subset <- cell.islkh
  } else if (obj == 'islkdox') {
    cell.subset <- cell.islkdox
  } else if (obj == 'islkIDN') {
    cell.subset <- cell.islkIDN
  } else {
    cell.subset <- cell.islkIFN
  }
  for (i in 0:num.clusters) {
    print(paste0('Getting values for cluster ', as.character(i)))
    subset <- c(mean(subset(cell.subset$percent.latent, cell.subset$seurat_clusters == i)),
                mean(subset(cell.subset$percent.early, cell.subset$seurat_clusters == i)),
                #mean(subset(cell.subset$percent.intermediate, cell.subset$seurat_clusters == i)),
                mean(subset(cell.subset$percent.late, cell.subset$seurat_clusters == i)))
    mean.KSHV <- c(mean.KSHV, rep(mean(subset(cell.subset$percent.KSHV, cell.subset$seurat_clusters == i)), 3))
    deviation <- c(sd(subset(cell.subset$percent.latent, cell.subset$seurat_clusters == i)),
                   sd(subset(cell.subset$percent.early, cell.subset$seurat_clusters == i)),
                   #sd(subset(cell.subset$percent.intermediate, cell.subset$seurat_clusters == i)),
                   sd(subset(cell.subset$percent.late, cell.subset$seurat_clusters == i)))
    deviation.KSHV <- c(deviation.KSHV, deviation)
    value <- c(value, subset)
  }
  data <- data.frame(clusters,stage,value,mean.KSHV,deviation.KSHV)
  write.csv(data, file = paste('~/Documents/Genomics/reanalyze_with_insertion/viral_genes/viral_stage_ordered/', obj, '.csv', sep = ''),
              row.names = F, col.names = T, quote = F)
  p1 <- ggplot(data, aes(fill=stage, y=value, x=reorder(clusters, mean.KSHV))) + geom_bar(position="stack", stat="identity") + theme_bw() +
    geom_point(aes(y=mean.KSHV, x=reorder(clusters, mean.KSHV))) +
    theme(panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
  postscript(file = paste('~/Documents/Genomics/reanalyze_with_insertion/viral_genes/viral_stage_ordered/', obj,'_withPAN.eps', sep = ''), width = 6, height = 4)
  print(p1)
  dev.off()
})

percent.cluster.KSHV <- c()
barcodes <- c()
for (i in 0:num.clusters) {
  subset.cluster <- subset(cell.islkdox$percent.KSHV, cell.islkdox$seurat_clusters == i)
  barcodes <- c(barcodes, rownames(subset(cell.islkdox@meta.data, seurat_clusters == i)))
  num.cells <- length(subset.cluster)
  percent.cluster.KSHV <- c(percent.cluster.KSHV, rep(mean(subset.cluster),num.cells))
}
df.percent.cluster <- data.frame(percent.cluster.KSHV, row.names = barcodes)
row.names(percent.cluster.KSHV) <- barcodes
cell.islkdox <- AddMetaData(cell.islkdox, df.percent.cluster, col.name = 'mean.KSHV')

data <- cell.islkdox@meta.data[order(cell.islkdox$percent.KSHV),] %>% rownames_to_column('cells') %>% pivot_longer(c('percent.latent','percent.early','percent.late'), names_to = 'stage', values_to = 'percent')
viral.plot1 <- ggplot(data, aes(fill=stage, y=percent, x=reorder(cells, percent.KSHV))) + geom_bar(position="stack", stat="identity") +
  geom_point(aes(y = percent.KSHV, x=reorder(cells, percent.KSHV))) + theme_bw() +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

postscript(file = paste('~/Documents/Genomics/reanalyze_with_insertion/viral_genes/viral_stage_ordered/iSLKdox_all_cells_stack.eps', sep = ''), width = 6, height = 4)
print(viral.plot1)
dev.off()

