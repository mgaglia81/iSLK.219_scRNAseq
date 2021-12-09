library(monocle3)
library(ggplot2)
library(dplyr)

# Loading data into monocle object
print('Loading 10X data...')
islkhm.data <- load_cellranger_data('~/Documents/Genomics/iSLKH_iSLK219_count/')
islkhm.data$sampleID <- 'A'
islkdoxm.data <- load_cellranger_data('~/Documents/Genomics/iSLK219_dox_count/')
islkdoxm.data$sampleID <- 'B'
islk_IDNm.data <- load_cellranger_data('~/Documents/Genomics/iSLK219_IDN_count/')
islk_IDNm.data$sampleID <- 'C'
islk_IFNm.data <- load_cellranger_data('~/Documents/Genomics/iSLK219_IDN_IFN_count/')
islk_IFNm.data$sampleID <- 'D'
cellm.combined <- combine_cds(list(islkhm.data, islkdoxm.data, islk_IDNm.data, islk_IFNm.data))

print('Processing 10X data...')
cellm.combined <- preprocess_cds(cellm.combined, num_dim = 30)
cellm.pc.plot <- plot_pc_variance_explained(cellm.combined)

print('Removing batch effects...')
cellm.combined <- align_cds(cellm.combined, alignment_group = 'sampleID')

print('Performing dimensional reduction and clustering...')
cellm.combined <- reduce_dimension(cellm.combined)
cellm.combined <- cluster_cells(cellm.combined)
print(plot_cells(cellm.combined, color_cells_by = 'sampleID'))

print('Constructing single cell trajectories...')
cellm.combined <- learn_graph(cellm.combined)
cellm.combined <- order_cells(cellm.combined)
plot_cells(cellm.combined)

saveRDS(cellm.combined, file = '~/Documents/Genomics/Monocle/cellm-combined.rds')

print('Performing differential gene expression analysis...')
cellm.graph.test <- graph_test(cellm.combined, neighbor_graph = 'principal_graph',cores=4)
cellm.deg.ids <- row.names(subset(cellm.graph.test, q_value < 0.05))
head(cellm.graph.test)
head(cellm.deg.ids)