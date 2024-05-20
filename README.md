# Impact of Benzo[a]pyrene and Sex Hormones on Human ESC-Derived Cells Using Single Cell Transcriptomics


## Packages

```r
library(tidyverse)
library(magrittr)
library(Seurat)
library(scCustomize)
library(ggraph)
library(igraph)
library(HGNChelper)
library(tuple)
library(ggvenn)
library(scater)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(data.tree)
library(ReactomePA)
library(topGO)
library(Rgraphviz)
library(readxl)
library(dplyr)
library(gridExtra)
library(reshape2)
library(UpSetR)
library(shadowtext)
```r

---

## Data Pre-processing


```r
############################ Create Seurat object ##############################

FDR = 0.05
`%notin%` = Negate(`%in%`)

# Load the dataset
Untr_Ctrl.data <- Read10X(data.dir = 'data/Untr_Ctrl/raw_feature_bc_matrix/')
Untr_Ctrl.data <- CreateSeuratObject(counts = Untr_Ctrl.data, project = 'Untr_Ctrl')
Untr_Ctrl.data

Untr_BaP.data <- Read10X(data.dir = 'data/Untr_BaP/raw_feature_bc_matrix/')
Untr_BaP.data <- CreateSeuratObject(counts = Untr_BaP.data, project = 'Untr_BaP')
Untr_BaP.data

Male_Ctrl.data <- Read10X(data.dir = 'data/Male_Ctrl/raw_feature_bc_matrix/')
Male_Ctrl.data <- CreateSeuratObject(counts = Male_Ctrl.data, project = 'Male_Ctrl')
Male_Ctrl.data

Male_BaP.data <- Read10X(data.dir = 'data/Male_BaP/raw_feature_bc_matrix/')
Male_BaP.data <- CreateSeuratObject(counts = Male_BaP.data, project = 'Male_BaP')
Male_BaP.data

Female_Ctrl.data <- Read10X(data.dir = 'data/Female_Ctrl/raw_feature_bc_matrix/')
Female_Ctrl.data <- CreateSeuratObject(counts = Female_Ctrl.data, project = 'Female_Ctrl')
Female_Ctrl.data

Female_BaP.data <- Read10X(data.dir = 'data/Female_BaP/raw_feature_bc_matrix/')
Female_BaP.data <- CreateSeuratObject(counts = Female_BaP.data, project = 'Female_BaP')
Female_BaP.data

# sample names
samples <- c('Untr_Ctrl', 'Untr_BaP', 
            'Male_Ctrl','Male_BaP', 
            'Female_Ctrl', 'Female_BaP')

# initialize the Seurat object with the raw (non-normalized data).
sobj.before <- merge(Female_Ctrl.data, y = c(Female_BaP.data, Untr_Ctrl.data, Untr_BaP.data, Male_Ctrl.data, Male_BaP.data), 
             add.cell.ids = c('Female_Ctrl', 'Female_BaP', 'Untr_Ctrl', 'Untr_BaP', 'Male_Ctrl', 'Male_BaP'), 
             project = 'human_esc')

# add sample info to metadata
DefaultAssay(sobj.before) <- ASSAY = 'RNA'
sobj.before@meta.data %<>% dplyr::rename(sample = 'orig.ident') %>% 
  separate(col = 'sample', into = c('hormones', 'EDC'), sep = '_', remove = F)

# calculate % of reads mapping to mitochondrial DNA 
sobj.before[['percent.mt']] <- PercentageFeatureSet(sobj.before, pattern = '^MT-')


################################################################################
# Filter for number of features (i.e. genes), count (i.e. read count per cell) #
########## and percentage of reads mapping to the mitochondrial genome #########


sobj <- subset(sobj.before, subset = nFeature_RNA > 2000 & nFeature_RNA < 10000 & nCount_RNA > 300 & nCount_RNA < 1e+05 & percent.mt < 12)

p1 <- VlnPlot(sobj, features = 'nFeature_RNA', raster = FALSE, assay = 'RNA', pt.size = 0, group.by = 'sample') + 
scale_fill_manual(values = c('Female_BaP' = '#2775b6', 'Female_Ctrl' = '#93d5dc', 
                               'Male_BaP' = '#41b349', 'Male_Ctrl' = '#add5a2', 
                               'Untr_BaP' = '#fcc307', 'Untr_Ctrl' = '#f7de98')) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
p2 <- VlnPlot(sobj, features = 'nCount_RNA', raster = FALSE, assay = 'RNA', pt.size = 0, group.by = 'sample') + 
scale_fill_manual(values = c('Female_BaP' = '#2775b6', 'Female_Ctrl' = '#93d5dc', 
                               'Male_BaP' = '#41b349', 'Male_Ctrl' = '#add5a2', 
                               'Untr_BaP' = '#fcc307', 'Untr_Ctrl' = '#f7de98')) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
p3 <- VlnPlot(sobj, features = 'percent.mt', raster = FALSE, assay = 'RNA', pt.size = 0, group.by = 'sample') + 
      scale_fill_manual(values = c('Female_BaP' = '#2775b6', 'Female_Ctrl' = '#93d5dc', 
                               'Male_BaP' = '#41b349', 'Male_Ctrl' = '#add5a2', 
                               'Untr_BaP' = '#fcc307', 'Untr_Ctrl' = '#f7de98')) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
patchwork::wrap_plots(p1, p2, p3, ncol = 3)
ggsave('plots/after_filter_vln.png', width = 16, height = 8, unit = 'in', dpi = 300)


# Cell numbers after filtering
Sample <- unique(sobj$sample)
cell_count <- data.frame(Sample = character(), Count = numeric())
for (n in Sample) {
  count2 <- subset(sobj, subset = sample == n)
  count1 <- length(Cells(count2))
  cell_count <- rbind(cell_count, data.frame(Sample = n, Count = count1))
}
print(cell_count)
ggplot(data = cell_count, aes(x = Sample, y = Count, fill = Sample)) + 
  geom_bar(stat = 'identity', color = "black", size = 0.5) +  
  labs(title = "Cell Count") +  
  theme_bw() +
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_text(hjust = 1, vjust = .5, angle = 90),
        legend.key.size = unit(0.5, "cm"),
        panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5)) +  
  scale_fill_manual(values = c('Female_BaP' = '#2775b6', 'Female_Ctrl' = '#93d5dc',  
                               'Male_BaP' = '#41b349', 'Male_Ctrl' = '#add5a2',
                               'Untr_Ctrl' = '#f7de98', 'Untr_BaP' = '#fcc307')) +
  guides(fill = guide_legend(override.aes = list(color = "black", size = 1)))
ggsave('plots/cells_after_filter.png', width = 10, height = 16, units = 'cm', dpi = 300)


################################################################################
################ Normalize data, run PCA and divide into clusters ##############

# normalize data
sobj[['percent.mt']] = PercentageFeatureSet(sobj, pattern = '^MT-')
sobj <- NormalizeData(sobj, normalization.method = 'LogNormalize', scale.factor = 10000)
sobj <- FindVariableFeatures(sobj, selection.method = 'vst', nfeatures = 2000)
sobj@meta.data


# scale and run PCA
sobj <- ScaleData(sobj, features = rownames(sobj))
sobj <- RunPCA(sobj, features = VariableFeatures(object = sobj))

# investigate genes influencing the loading scores and visualise PCA
print(sobj[['pca']], dims = 1:5, nfeatures = 5)
VizDimLoadings(sobj, dims = 1:3, reduction = 'pca')
DimPlot(sobj, reduction = 'pca')
ggsave('/plots/PCA_plot.png', width = 8, height = 8, unit = 'in', dpi = 300)


# first construct a KNN graph based on the euclidean distance in PCA space, and refine the edge weights between any two cells based on the shared overlap in their local neighborhoods (Jaccard similarity).
# input: selected amount of PCs
sobj <- FindNeighbors(sobj, dims = 1:18)

# modularity optimization techniques: Louvain algorithm (default) or SLM
# setting the resolution between 0.4-1.2 typically returns good results for single-cell datasets of around 3K cells. Optimal resolution often increases for larger datasets.
sobj <- FindClusters(sobj, resolution = 0.6) 
# determines the amount of clusters to use more precisely then reading the elbowplot. Resolution = between 0.4 (least resolution) and 1.2 (max resolution)
clusters <- Idents(sobj)

sobj <- RunUMAP(sobj, dims = 1:length(levels(clusters)))
DimPlot(sobj, reduction = 'umap', label = T, seed = 42, group.by = 'seurat_clusters')
```r


## Data Pre-processing












