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
```

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
################ Normalise data, run PCA and divide into clusters ##############

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
```

---

## Cluster Annotation

```r
############################## Annotate clusters ###############################

# load gene set preparation function
source("gene_sets_prepare.R")
# load cell type annotation function
source("sctype_score_.R")

# DB file
db <- 'ScTypeDB_full - thyroid.xlsx'
tissue <- c('Thyroid')

# prepare gene sets
gs_list <- gene_sets_prepare(db,tissue)

# check Seurat object version (scRNA-seq matrix extracted differently in Seurat v4/v5)
#seurat_package_v5 <- isFALSE('counts' %in% names(attributes(sobj[["RNA"]])));
#print(sprintf("Seurat object %s is used", ifelse(seurat_package_v5, "v5", "v4")))

# extract scaled scRNA-seq matrix
scRNAseqData_scaled <- if (seurat_package_v5) as.matrix(sobj[["RNA"]]$scale.data) else as.matrix(sobj[["RNA"]]@scale.data)

# run ScType
es.max <- sctype_score(scRNAseqData = scRNAseqData_scaled, scaled = TRUE, gs = gs_list$gs_positive, gs2 = gs_list$gs_negative)

# merge by cluster
cL_resutls = do.call('rbind', lapply(unique(sobj@meta.data$seurat_clusters), function(cl){
  es.max.cl = sort(rowSums(es.max[ , rownames(sobj@meta.data[sobj@meta.data$seurat_clusters == cl, ])]), decreasing = !0)
  head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(sobj@meta.data$seurat_clusters == cl)), 10)
}))
sctype_scores = cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores)

# set low-confident (low ScType score) clusters to 'unknown'
sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] = 'Unknown'
print(sctype_scores[, 1:3])
# Overlay the annotation on the TSNE
sobj@meta.data$customclassif = ''
for(j in unique(sctype_scores$cluster)){
  cl_type = sctype_scores[sctype_scores$cluster == j,]; 
  sobj@meta.data$customclassif[sobj@meta.data$seurat_clusters == j] = as.character(cl_type$type[1])
}
DimPlot(sobj, reduction = 'umap', label = TRUE, seed = 42, repel = T, group.by = c('customclassif')) 
ggsave('/plots/UMAP_classification.png', width = 9, height = 8, unit = 'in', dpi = 300)


# Numbers each cell type
Cell <- unique(sobj@meta.data$customclassif)
cell_type <- data.frame(Cell = character(), Count = numeric())
for (n in Cell) {
  count1 <- sum(sobj@meta.data$customclassif == n)
  cell_type <- rbind(cell_type, data.frame(Cell = n, Count = count1))
}
print(cell_type)
ggplot(data = cell_type, aes(x = Cell, y = Count, fill = Cell)) + 
  geom_bar(stat = 'identity') +  
  labs(title = "Cell Type") +  
  theme_bw() +
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_text(hjust = 1, vjust = .5, angle = 90),
        legend.key.size = unit(0.5, "cm"),
        panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5))+
  guides(fill = guide_legend(override.aes = list(shape = 5)))
ggsave('plots/numbers_each__cell_type.png', width = 10, height = 16, units = 'cm', dpi = 300)
```

---

## Differential Gene Expression

```r
# Add a col to the @meta.data that concatenates cell type and hormone treatment and EDC treatment
sobj$cell_hormone_edc = paste(sobj@meta.data[, 'customclassif'],
                              sobj@meta.data[, 'hormones'],
                              sobj@meta.data[, 'EDC'], sep = '_')

sobj$cell_hormone = paste(sobj@meta.data[, 'customclassif'],
                          sobj@meta.data[, 'hormones'], sep = '_')

Idents(sobj) = sobj$customclassif


######################### Airway Progenitor Cells ##############################

Airway_Progenitor_Cells_sobj <- subset(sobj, subset = customclassif == 'Airway Progenitor Cells') 
Airway_Progenitor_Cells_sobj <- JoinLayers(Airway_Progenitor_Cells_sobj, what = c("RNA", "scale.data"))
test_comparisons = list(
  male = c('Airway progenitor cells_Male_Ctrl', 'Airway progenitor cells_Untr_Ctrl'),
  female = c('Airway progenitor cells_Female_Ctrl', 'Airway progenitor cells_Untr_Ctrl'),
  male_female = c('Airway progenitor cells_Male_Ctrl', 'Airway progenitor cells_Female_Ctrl'),
  bapResponse_affectedByMale = c('Airway progenitor cells_Male_BaP', 'Airway progenitor cells_Untr_BaP'),
  bapResponse_affectedByFemale = c('Airway progenitor cells_Female_BaP', 'Airway progenitor cells_Untr_BaP'),
  bapResponse = c('Airway progenitor cells_Untr_BaP', 'Airway progenitor cells_Untr_Ctrl'),
  bapResponse_inMaleEnv = c('Airway progenitor cells_Male_BaP', 'Airway progenitor cells_Male_Ctrl'),
  bapResponse_inFemaleEnv = c('Airway progenitor cells_Female_BaP', 'Airway progenitor cells_Female_Ctrl')
)
# where to store airway results
degs_airway = list() # degs
genes_airway = list() # airway genes
for (x in 1:length(test_comparisons)) {
  
  print(paste0('Comparing ', test_comparisons[[x]][1], ' VS ', test_comparisons[[x]][2]))
  
  # find DEGs
  airway_genes = FindMarkers(Airway_Progenitor_Cells_sobj, ident.1 = test_comparisons[[x]][1], ident.2 = test_comparisons[[x]][2], group.by = 'cell_hormone_edc')

  # assign DEGs results to list that collects them all
  degs_airway[[x]] = filter(airway_genes, p_val_adj < FDR)
  genes_airway[[x]] = airway_genes
  names(degs_airway)[[x]] = paste0(test_comparisons[[x]][1], 'VS', test_comparisons[[x]][2])
  names(genes_airway)[[x]] = paste0(test_comparisons[[x]][1], 'VS', test_comparisons[[x]][2])
  
  # also, save results as a table
  # airway genes
  write.csv(genes_airway[[x]], paste0(test_comparisons[[x]][1], '_', test_comparisons[[x]][2], 'All_genes.csv'))

  # DEGs
  write.csv(degs_airway[[x]], paste0(test_comparisons[[x]][1], '_', test_comparisons[[x]][2], 'DEGs_FDR', FDR, '.csv'))

  # how many genes have been tested and how many DEGs (FDR < 0.05) there are
  num_genes = nrow(genes_airway[[x]])
  num_DEGs = nrow(degs_airway[[x]])
  DEG_percent = (num_DEGs/num_genes)*100
  
  # give some info on number of tested genes and DEGs
  print(paste0('Found ', num_DEGs, ' unique DEGs (', num_genes, ' genes tested) when comparing ', test_comparisons[[x]][1], ' - ', test_comparisons[[x]][2], '; ', round(DEG_percent, digits = 1), '% is differentially expressed.'))}

# save all tested airway genes and DEGs in two excel files
degs_airway_save <- degs_airway
for (x in 1:length(degs_airway)) {
  names(degs_airway_save)[x] <- paste0(map_chr(str_split(test_comparisons[[x]][1], 'Airway progenitor cells_'), 2), 
                                    'VS',
                                    map_chr(str_split(test_comparisons[[x]][2], 'Airway progenitor cells_'), 2))
  
  degs_airway_save[[x]] <- degs_airway_save[[x]] %>% rownames_to_column(var = 'gene')
}
writexl::write_xlsx(degs_airway_save, 'DEGs_airway_comparisons.xlsx')
genes_airway_save = genes_airway
for (x in 1:length(genes_airway)) {
  names(genes_airway_save)[x] = paste0(map_chr(str_split(test_comparisons[[x]][1], 'Airway progenitor cells_'), 2), 
                                    'VS',
                                    map_chr(str_split(test_comparisons[[x]][2], 'Airway progenitor cells_'), 2))
  
  genes_airway_save[[x]] %<>% rownames_to_column(var = 'gene')
}
writexl::write_xlsx(genes_airway_save, all_genes_airway_comparisons.xlsx')


################################ Cardiovascular Cells ##################################

Cardiovascular_Cells_sobj = subset(sobj, subset = customclassif == 'Cardiovascular Cells') 
Cardiovascular_Cells_sobj <- JoinLayers(Cardiovascular_Cells_sobj, what = c("RNA", "scale.data"))
test_comparisons = list(
  male = c('Cardiovascular cells_Male_Ctrl', 'Cardiovascular cells_Untr_Ctrl'),
  female = c('Cardiovascular cells_Female_Ctrl', 'Cardiovascular cells_Untr_Ctrl'),
  male_female = c('Cardiovascular cells_Male_Ctrl', 'Cardiovascular cells_Female_Ctrl'),
  bapResponse_affectedByMale = c('Cardiovascular cells_Male_BaP', 'Cardiovascular cells_Untr_BaP'),
  bapResponse_affectedByFemale = c('Cardiovascular cells_Female_BaP', 'Cardiovascular cells_Untr_BaP'),
  bapResponse = c('Cardiovascular cells_Untr_BaP', 'Cardiovascular cells_Untr_Ctrl'),
  bapResponse_inMaleEnv = c('Cardiovascular cells_Male_BaP', 'Cardiovascular cells_Male_Ctrl'),
  bapResponse_inFemaleEnv = c('Cardiovascular cells_Female_BaP', 'Cardiovascular cells_Female_Ctrl')
)
# where to store all results
degs_cardiovascular = list() # degs
genes_cardiovascular = list() # cardiovascular genes
for (x in 1:length(test_comparisons)) {
  
  print(paste0('Comparing ', test_comparisons[[x]][1], ' VS ', test_comparisons[[x]][2]))
  
  # find DEGs
  cardiovascular_genes = FindMarkers(Cardiovascular_Cells_sobj, ident.1 = test_comparisons[[x]][1], ident.2 = test_comparisons[[x]][2], group.by = 'cell_hormone_edc')

  # assign DEGs results to list that collects them all
  degs_cardiovascular[[x]] = filter(cardiovascular_genes, p_val_adj < FDR)
  genes_cardiovascular[[x]] = cardiovascular_genes
  names(degs_cardiovascular)[[x]] = paste0(test_comparisons[[x]][1], 'VS', test_comparisons[[x]][2])
  names(genes_cardiovascular)[[x]] = paste0(test_comparisons[[x]][1], 'VS', test_comparisons[[x]][2])
  
  # also, save results as a table
  # all genes
  write.csv(degs_cardiovascular[[x]], paste0(test_comparisons[[x]][1], '_', test_comparisons[[x]][2], 'All_genes.csv'))
  
  # DEGs
  write.csv(degs_cardiovascular[[x]], paste0(test_comparisons[[x]][1], '_', test_comparisons[[x]][2], 'DEGs_FDR', FDR, '.csv'))
  
  # how many genes have been tested and how many DEGs (FDR < 0.05) there are
  num_genes = nrow(genes_cardiovascular[[x]])
  num_DEGs = nrow(degs_cardiovascular[[x]])
  DEG_percent = (num_DEGs/num_genes)*100
  
  # give some info on number of tested genes and DEGs
  print(paste0('Found ', num_DEGs, ' unique DEGs (', num_genes, ' genes tested) when comparing ', test_comparisons[[x]][1], ' - ', test_comparisons[[x]][2], '; ', round(DEG_percent, digits = 1), '% is differentially expressed.'))}

# save all tested cardiovascular genes and DEGs in two excel files
degs_cardiovascular_save <- degs_cardiovascular
for (x in 1:length(degs_cardiovascular)) {
  names(degs_cardiovascular)[x] <- paste0(map_chr(str_split(test_comparisons[[x]][1], 'Cardiovascular cells_'), 2), 
                                       'VS',
                                       map_chr(str_split(test_comparisons[[x]][2], 'Cardiovascular cells_'), 2))
  
  degs_cardiovascular[[x]] <- degs_cardiovascular[[x]] %>% rownames_to_column(var = 'gene')
}
writexl::write_xlsx(degs_cardiovascular, 'DEGs_cardiovascular_comparisons.xlsx')
genes_cardiovascular_save = genes_cardiovascular
for (x in 1:length(genes_cardiovascular)) {
  names(genes_cardiovascular_save)[x] = paste0(map_chr(str_split(test_comparisons[[x]][1], 'Cardiovascular cells_'), 2), 
                                       'VS',
                                       map_chr(str_split(test_comparisons[[x]][2], 'Cardiovascular cells_'), 2))
  
  genes_cardiovascular_save[[x]] %<>% rownames_to_column(var = 'gene')
}
writexl::write_xlsx(genes_cardiovascular_save, 'all_genes_cardiovascular_comparisons.xlsx')


############################# Embryonic Stem Cells #############################

Embryonic_Stem_Cells_sobj = subset(sobj, subset = customclassif == 'Embryonic Stem Cells') 
Embryonic_Stem_Cells_sobj <- JoinLayers(Embryonic_Stem_Cells_sobj, what = c("RNA", "scale.data"))
test_comparisons = list(
  male = c('Embryonic stem cells_Male_Ctrl', 'Embryonic stem cells_Untr_Ctrl'),
  female = c('Embryonic stem cells_Female_Ctrl', 'Embryonic stem cells_Untr_Ctrl'),
  male_female = c('Embryonic stem cells_Male_Ctrl', 'Embryonic stem cells_Female_Ctrl'),
  bapResponse_affectedByMale = c('Embryonic stem cells_Male_BaP', 'Embryonic stem cells_Untr_BaP'),
  bapResponse_affectedByFemale = c('Embryonic stem cells_Female_BaP', 'Embryonic stem cells_Untr_BaP'),
  bapResponse = c('Embryonic stem cells_Untr_BaP', 'Embryonic stem cells_Untr_Ctrl'),
  bapResponse_inMaleEnv = c('Embryonic stem cells_Male_BaP', 'Embryonic stem cells_Male_Ctrl'),
  bapResponse_inFemaleEnv = c('Embryonic stem cells_Female_BaP', 'Embryonic stem cells_Female_Ctrl')
)
# where to store all results
degs_embryonic = list() # degs
genes_embryonic = list() # embryonic genes
for (x in 1:length(test_comparisons)) {
  
  print(paste0('Comparing ', test_comparisons[[x]][1], ' VS ', test_comparisons[[x]][2]))
  
  # find DEGs
  embryonic_genes = FindMarkers(Embryonic_Stem_Cells_sobj, ident.1 = test_comparisons[[x]][1], ident.2 = test_comparisons[[x]][2], group.by = 'cell_hormone_edc')
  
  # assign DEGs results to list that collects them all
  degs_embryonic[[x]] = filter(embryonic_genes, p_val_adj < FDR)
  genes_embryonic[[x]] = embryonic_genes
  names(degs_embryonic)[[x]] = paste0(test_comparisons[[x]][1], 'VS', test_comparisons[[x]][2])
  names(genes_embryonic)[[x]] = paste0(test_comparisons[[x]][1], 'VS', test_comparisons[[x]][2])
  
  # also, save results as a table
  # embryonic genes
  write.csv(genes_embryonic[[x]], paste0(test_comparisons[[x]][1], '_', test_comparisons[[x]][2], 'All_genes.csv'))
  
  # DEGs
  write.csv(degs_embryonic[[x]], paste0(test_comparisons[[x]][1], '_', test_comparisons[[x]][2], 'DEGs_FDR', FDR, '.csv'))
  
  # how many genes have been tested and how many DEGs (FDR < 0.05) there are
  num_genes = nrow(genes_embryonic[[x]])
  num_DEGs = nrow(degs_embryonic[[x]])
  DEG_percent = (num_DEGs/num_genes)*100
  
  # give some info on number of tested genes and DEGs
  print(paste0('Found ', num_DEGs, ' unique DEGs (', num_genes, ' genes tested) when comparing ', test_comparisons[[x]][1], ' - ', test_comparisons[[x]][2], '; ', round(DEG_percent, digits = 1), '% is differentially expressed.'))}

# save all tested genes and DEGs in two excel files
degs_embryonic_save <- degs_embryonic
for (x in 1:length(degs_embryonic)) {
  names(degs_embryonic_save)[x] <- paste0(map_chr(str_split(test_comparisons[[x]][1], 'Embryonic stem cells_'), 2), 
                                    'VS',
                                    map_chr(str_split(test_comparisons[[x]][2], 'Embryonic stem cells_'), 2))
  
  degs_embryonic_save[[x]] <- degs_embryonic_save[[x]] %>% rownames_to_column(var = 'gene')
}
writexl::write_xlsx(degs_embryonic_save, 'DEGs_embryonic_comparisons.xlsx')
genes_embryonic_save = genes_embryonic
for (x in 1:length(genes_embryonic)) {
  names(genes_embryonic_save)[x] = paste0(map_chr(str_split(test_comparisons[[x]][1], 'Embryonic stem cells_'), 2), 
                                    'VS',
                                    map_chr(str_split(test_comparisons[[x]][2], 'Embryonic stem cells_'), 2))
  
  genes_embryonic_save[[x]] %<>% rownames_to_column(var = 'gene')
}
writexl::write_xlsx(genes_embryonic_save, 'all_genes_embryonic_comparisons.xlsx')


############################### Epithelial Cells ###############################
Epithelial_Cells_sobj = subset(sobj, subset = customclassif == 'Epithelial Cells') 
Epithelial_Cells_sobj <- JoinLayers(Epithelial_Cells_sobj, what = c("RNA", "scale.data"))
test_comparisons = list(
  male = c('Epithelial cells_Male_Ctrl', 'Epithelial cells_Untr_Ctrl'),
  female = c('Epithelial cells_Female_Ctrl', 'Epithelial cells_Untr_Ctrl'),
  male_female = c('Epithelial cells_Male_Ctrl', 'Epithelial cells_Female_Ctrl'),
  bapResponse_affectedByMale = c('Epithelial cells_Male_BaP', 'Epithelial cells_Untr_BaP'),
  bapResponse_affectedByFemale = c('Epithelial cells_Female_BaP', 'Epithelial cells_Untr_BaP'),
  bapResponse = c('Epithelial cells_Untr_BaP', 'Epithelial cells_Untr_Ctrl'),
  bapResponse_inMaleEnv = c('Epithelial cells_Male_BaP', 'Epithelial cells_Male_Ctrl'),
  bapResponse_inFemaleEnv = c('Epithelial cells_Female_BaP', 'Epithelial cells_Female_Ctrl')
)
# where to store epithelial results
degs_epithelial = list() # degs
genes_epithelial = list() # epithelial genes
for (x in 1:length(test_comparisons)) {
  
  print(paste0('Comparing ', test_comparisons[[x]][1], ' VS ', test_comparisons[[x]][2]))
  
  # find DEGs
  epithelial_genes = FindMarkers(Epithelial_Cells_sobj, ident.1 = test_comparisons[[x]][1], ident.2 = test_comparisons[[x]][2], group.by = 'cell_hormone_edc')
  
  # assign DEGs results to list that collects them all
  degs_epithelial[[x]] = filter(epithelial_genes, p_val_adj < FDR)
  genes_epithelial[[x]] = epithelial_genes
  names(degs_epithelial)[[x]] = paste0(test_comparisons[[x]][1], 'VS', test_comparisons[[x]][2])
  names(genes_epithelial)[[x]] = paste0(test_comparisons[[x]][1], 'VS', test_comparisons[[x]][2])
  
  # also, save results as a table
  # epithelial genes
  write.csv(genes_epithelial[[x]], paste0(test_comparisons[[x]][1], '_', test_comparisons[[x]][2], 'All_genes.csv'))
  
  # DEGs
  write.csv(degs_epithelial[[x]], paste0(test_comparisons[[x]][1], '_', test_comparisons[[x]][2], 'DEGs_FDR', FDR, '.csv'))
  
  # how many genes have been tested and how many DEGs (FDR < 0.05) there are
  num_genes = nrow(genes_epithelial[[x]])
  num_DEGs = nrow(degs_epithelial[[x]])
  DEG_percent = (num_DEGs/num_genes)*100
  
  # give some info on number of tested genes and DEGs
  print(paste0('Found ', num_DEGs, ' unique DEGs (', num_genes, ' genes tested) when comparing ', test_comparisons[[x]][1], ' - ', test_comparisons[[x]][2], '; ', round(DEG_percent, digits = 1), '% is differentially expressed.'))}

# save all tested epithelial genes and DEGs in two excel files
degs_epithelial_save <- degs_epithelial
for (x in 1:length(degs_epithelial)) {
  names(degs_epithelial_save)[x] <- paste0(map_chr(str_split(test_comparisons[[x]][1], 'Epithelial cells_'), 2), 
                                    'VS',
                                    map_chr(str_split(test_comparisons[[x]][2], 'Epithelial cells_'), 2))
  
  degs_epithelial_save[[x]] <- degs_epithelial_save[[x]] %>% rownames_to_column(var = 'gene')
}
writexl::write_xlsx(degs_epithelial_save, 'DEGs_epithelial_comparisons.xlsx')
genes_epithelial_save = genes_epithelial
for (x in 1:length(genes_epithelial)) {
  names(genes_epithelial_save)[x] = paste0(map_chr(str_split(test_comparisons[[x]][1], 'Epithelial cells_'), 2), 
                                    'VS',
                                    map_chr(str_split(test_comparisons[[x]][2], 'Epithelial cells_'), 2))
  
  genes_epithelial_save[[x]] %<>% rownames_to_column(var = 'gene')
}
writexl::write_xlsx(genes_epithelial_save, 'all_genes_epithelial_comparisons.xlsx')


################################## Fibroblasts #################################

Fibroblasts_sobj = subset(sobj, subset = customclassif == 'Fibroblasts') 
Fibroblasts_sobj <- JoinLayers(Fibroblasts_sobj, what = c("RNA", "scale.data"))
test_comparisons = list(
  male = c('Fibroblasts_Male_Ctrl', 'Fibroblasts_Untr_Ctrl'),
  female = c('Fibroblasts_Female_Ctrl', 'Fibroblasts_Untr_Ctrl'),
  male_female = c('Fibroblasts_Male_Ctrl', 'Fibroblasts_Female_Ctrl'),
  bapResponse_affectedByMale = c('Fibroblasts_Male_BaP', 'Fibroblasts_Untr_BaP'),
  bapResponse_affectedByFemale = c('Fibroblasts_Female_BaP', 'Fibroblasts_Untr_BaP'),
  bapResponse = c('Fibroblasts_Untr_BaP', 'Fibroblasts_Untr_Ctrl'),
  bapResponse_inMaleEnv = c('Fibroblasts_Male_BaP', 'Fibroblasts_Male_Ctrl'),
  bapResponse_inFemaleEnv = c('Fibroblasts_Female_BaP', 'Fibroblasts_Female_Ctrl')
)
# where to store all results
degs_fibroblast = list() # degs
genes_fibroblast = list() # fibroblast genes
for (x in 1:length(test_comparisons)) {
  
  print(paste0('Comparing ', test_comparisons[[x]][1], ' VS ', test_comparisons[[x]][2]))
  
  # find DEGs
  fibroblast_genes = FindMarkers(Fibroblasts_sobj, ident.1 = test_comparisons[[x]][1], ident.2 = test_comparisons[[x]][2], group.by = 'cell_hormone_edc')
  
  # assign DEGs results to list that collects them all
  degs_fibroblast[[x]] = filter(fibroblast_genes, p_val_adj < FDR)
  genes_fibroblast[[x]] = fibroblast_genes
  names(degs_fibroblast)[[x]] = paste0(test_comparisons[[x]][1], 'VS', test_comparisons[[x]][2])
  names(genes_fibroblast)[[x]] = paste0(test_comparisons[[x]][1], 'VS', test_comparisons[[x]][2])
  
  # also, save results as a table
  # fibroblast genes
  write.csv(genes_fibroblast[[x]], paste0(test_comparisons[[x]][1], '_', test_comparisons[[x]][2], 'All_genes.csv'))
  
  # DEGs
  write.csv(degs_fibroblast[[x]], paste0(test_comparisons[[x]][1], '_', test_comparisons[[x]][2], 'DEGs_FDR', FDR, '.csv'))
  
  # how many genes have been tested and how many DEGs (FDR < 0.05) there are
  num_genes = nrow(genes_fibroblast[[x]])
  num_DEGs = nrow(degs_fibroblast[[x]])
  DEG_percent = (num_DEGs/num_genes)*100
  
  # give some info on number of tested genes and DEGs
  print(paste0('Found ', num_DEGs, ' unique DEGs (', num_genes, ' genes tested) when comparing ', test_comparisons[[x]][1], ' - ', test_comparisons[[x]][2], '; ', round(DEG_percent, digits = 1), '% is differentially expressed.'))}

# save all tested fibroblast genes and DEGs in two excel files
degs_fibroblast_save <- degs_fibroblast
for (x in 1:length(degs_fibroblast)) {
  names(degs_fibroblast_save)[x] <- paste0(map_chr(str_split(test_comparisons[[x]][1], 'Fibroblasts_'), 2), 
                                    'VS',
                                    map_chr(str_split(test_comparisons[[x]][2], 'Fibroblasts_'), 2))
  
  degs_fibroblast_save[[x]] <- degs_fibroblast_save[[x]] %>% rownames_to_column(var = 'gene')
}
writexl::write_xlsx(degs_fibroblast_save, 'DEGs_fibroblasts_comparisons.xlsx')
genes_fibroblast_save = genes_fibroblast
for (x in 1:length(genes_fibroblast)) {
  names(genes_fibroblast_save)[x] = paste0(map_chr(str_split(test_comparisons[[x]][1], 'Fibroblasts_'), 2), 
                                    'VS',
                                    map_chr(str_split(test_comparisons[[x]][2], 'Fibroblasts_'), 2))
  
  genes_fibroblast_save[[x]] %<>% rownames_to_column(var = 'gene')
}
writexl::write_xlsx(genes_fibroblast_save, 'all_genes_fibroblasts_comparisons.xlsx')


################################ Goblet Cells ##################################

Goblet_Cells_sobj = subset(sobj, subset = customclassif == 'Goblet Cells') 
Goblet_Cells_sobj <- JoinLayers(Goblet_Cells_sobj, what = c("RNA", "scale.data"))
test_comparisons = list(
  male = c('Goblet cells_Male_Ctrl', 'Goblet cells_Untr_Ctrl'),
  female = c('Goblet cells_Female_Ctrl', 'Goblet cells_Untr_Ctrl'),
  male_female = c('Goblet cells_Male_Ctrl', 'Goblet cells_Female_Ctrl'),
  bapResponse_affectedByMale = c('Goblet cells_Male_BaP', 'Goblet cells_Untr_BaP'),
  bapResponse_affectedByFemale = c('Goblet cells_Female_BaP', 'Goblet cells_Untr_BaP'),
  bapResponse = c('Goblet cells_Untr_BaP', 'Goblet cells_Untr_Ctrl'),
  bapResponse_inMaleEnv = c('Goblet cells_Male_BaP', 'Goblet cells_Male_Ctrl'),
  bapResponse_inFemaleEnv = c('Goblet cells_Female_BaP', 'Goblet cells_Female_Ctrl')
)
# where to store all results
degs_goblet = list() # degs
genes_goblet = list() # goblet genes
for (x in 1:length(test_comparisons)) {
  
  print(paste0('Comparing ', test_comparisons[[x]][1], ' VS ', test_comparisons[[x]][2]))
  
  # find DEGs
  goblet_genes = FindMarkers(Goblet_Cells_sobj, ident.1 = test_comparisons[[x]][1], ident.2 = test_comparisons[[x]][2], group.by = 'cell_hormone_edc')
  
  # assign DEGs results to list that collects them all
  degs_goblet[[x]] = filter(goblet_genes, p_val_adj < FDR)
  genes_goblet[[x]] = goblet_genes
  names(degs_goblet)[[x]] = paste0(test_comparisons[[x]][1], 'VS', test_comparisons[[x]][2])
  names(genes_goblet)[[x]] = paste0(test_comparisons[[x]][1], 'VS', test_comparisons[[x]][2])
  
  # also, save results as a table
  # all genes
  write.csv(genes_goblet[[x]], paste0(test_comparisons[[x]][1], '_', test_comparisons[[x]][2], 'All_genes.csv'))
  
  # DEGs
  write.csv(degs_goblet[[x]], paste0(test_comparisons[[x]][1], '_', test_comparisons[[x]][2], 'DEGs_FDR', FDR, '.csv'))
  
  # how many genes have been tested and how many DEGs (FDR < 0.05) there are
  num_genes = nrow(genes_goblet[[x]])
  num_DEGs = nrow(degs_goblet[[x]])
  DEG_percent = (num_DEGs/num_genes)*100
  
  # give some info on number of tested genes and DEGs
  print(paste0('Found ', num_DEGs, ' unique DEGs (', num_genes, ' genes tested) when comparing ', test_comparisons[[x]][1], ' - ', test_comparisons[[x]][2], '; ', round(DEG_percent, digits = 1), '% is differentially expressed.'))}

# save all tested goblet genes and DEGs in two excel files
degs_goblet_save <- degs_goblet
for (x in 1:length(degs_goblet)) {
  names(degs_goblet_save)[x] <- paste0(map_chr(str_split(test_comparisons[[x]][1], 'Goblet cells_'), 2), 
                                    'VS',
                                    map_chr(str_split(test_comparisons[[x]][2], 'Goblet cells_'), 2))
  
  degs_goblet_save[[x]] <- degs_goblet_save[[x]] %>% rownames_to_column(var = 'gene')
}
writexl::write_xlsx(degs_goblet_save, 'DEGs_goblet_comparisons.xlsx')
genes_goblet_save = genes_goblet
for (x in 1:length(genes_goblet)) {
  names(genes_goblet_save)[x] = paste0(map_chr(str_split(test_comparisons[[x]][1], 'Goblet cells_'), 2), 
                                    'VS',
                                    map_chr(str_split(test_comparisons[[x]][2], 'Goblet cells_'), 2))
  
  genes_goblet_save[[x]] %<>% rownames_to_column(var = 'gene')
}
writexl::write_xlsx(genes_goblet_save, 'all_genes_goblet_comparisons.xlsx')
```

---

## Gene Ontology (GO) analysis

```r
ontology = 'BP'
go_analysis = function(DEG_LIST, i, ...) {
  
  gene_list = bitr(rownames(DEG_LIST), fromType = 'SYMBOL', toType = 'ENTREZID', OrgDb = 'org.Hs.eg.db')
  
  # if some genes failed to map, report it
  print('These genes failed to map: ')
  print(DEG_LIST[rownames(DEG_LIST) %notin% gene_list$SYMBOL, 1:2])
  
  x = enrichGO(gene = gene_list$ENTREZ, 
               OrgDb = org.Hs.eg.db,
               qvalueCutoff = 0.01,
               pvalueCutoff = 0.01,
               ont = ontology, 
               readable = TRUE)  
  
  x = clusterProfiler::simplify(x, by = 'qvalue') 
  
  npathway = x@result%>% 
    as.data.frame %>% 
    filter(qvalue < 0.01) %>% 
    nrow() 
  
  print(paste0('Enriched GO pathways: ', npathway))
  
  # also, save results as a table 
  x@result %>% 
    as.data.frame %>% 
    filter(qvalue < 0.01) %>% 
    write.csv(., paste('GO/', ontology, 'simplified', test_comparisons[[i]][1], test_comparisons[[i]][2], ..., 'qvalue_0.01.csv', sep = '_'))
  
  return(x)
  
}


#################################### Airway ####################################
degs_airway_pct0.2 <- lapply(degs_airway, function(df) {subset(df, pct.1 >= 0.2 & pct.2 >= 0.2)})
test_comparisons = list(
  male = c('Airway progenitor cells_Male_Ctrl', 'Airway progenitor cells_Untr_Ctrl'),
  female = c('Airway progenitor cells_Female_Ctrl', 'Airway progenitor cells_Untr_Ctrl'),
  male_female = c('Airway progenitor cells_Male_Ctrl', 'Airway progenitor cells_Female_Ctrl'),
  bapResponse_affectedByMale = c('Airway progenitor cells_Male_BaP', 'Airway progenitor cells_Untr_BaP'),
  bapResponse_affectedByFemale = c('Airway progenitor cells_Female_BaP', 'Airway progenitor cells_Untr_BaP'),
  bapResponse = c('Airway progenitor cells_Untr_BaP', 'Airway progenitor cells_Untr_Ctrl'),
  bapResponse_inMaleEnv = c('Airway progenitor cells_Male_BaP', 'Airway progenitor cells_Male_Ctrl'),
  bapResponse_inFemaleEnv = c('Airway progenitor cells_Female_BaP', 'Airway progenitor cells_Female_Ctrl')
)
go_a_results <- list(
  go_a_MC_vs_UC <- go_analysis(degs_airway_pct0.2[[1]], 1, 'BP'),
  go_a_FC_vs_UC <- go_analysis(degs_airway_pct0.2[[2]], 2, 'BP'),
  go_a_MC_vs_FC <- go_analysis(degs_airway_pct0.2[[3]], 3, 'BP'),
  go_a_MB_vs_UB <- go_analysis(degs_airway_pct0.2[[4]], 4, 'BP'),
  go_a_FB_vs_UB <- go_analysis(degs_airway_pct0.2[[5]], 5, 'BP'),
  go_a_UB_vs_UC <- go_analysis(degs_airway_pct0.2[[6]], 6, 'BP'),
  go_a_MB_vs_MC <- go_analysis(degs_airway_pct0.2[[7]], 7, 'BP'),
  go_a_FB_vs_FC <- go_analysis(degs_airway_pct0.2[[8]], 8, 'BP'))
plot_a_list <- list()
for (i in seq_along(go_a_results)) {
  plot <- dotplot(go_a_results[[i]],
                  x = "GeneRatio/BgRatio",
                  color = "qvalue",
                  showCategory = 20,
                  size = NULL,
                  title = "CC_dotplot") +
                  scale_color_gradient(low = "blue", high = "grey90", name = "qvalue") +
                  ggtitle(paste(titles[i], sep = "")) +  
                  theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
                  labs(x = "Fold enrichment")
  
  plot_a_list[[i]] <- plot
}
combined_a_plot <- do.call(grid.arrange, c(plot_a_list, ncol = 4))
ggsave('plots/go_a_foldenrichment.png', combined_a_plot, width = 25, height = 21, unit = 'in', dpi = 300)
plot_a_list <- list()
for (i in seq_along(go_a_results)) {
  plot <- dotplot(go_a_results[[i]],
                  x = "GeneRatio",
                  color = "qvalue",
                  showCategory = 20,
                  size = NULL,
                  title = "CC_dotplot") +
    scale_color_gradient(low = "blue", high = "grey90", name = "qvalue") +
    ggtitle(paste(titles[i], sep = "")) +  
    theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
    labs(x = "GeneRatio")
  
  plot_a_list[[i]] <- plot
}
combined_a_plot <- do.call(grid.arrange, c(plot_a_list, ncol = 4))
ggsave('plots/go_a_generatio.png', combined_a_plot, width = 25, height = 21, unit = 'in', dpi = 300)


################################## Epithelial ##################################
degs_epithelial_pct0.2 <- lapply(degs_epithelial, function(df) {subset(df, pct.1 >= 0.2 & pct.2 >= 0.2)})
test_comparisons = list(
  male = c('Embryonic stem cells_Male_Ctrl', 'Embryonic stem cells_Untr_Ctrl'),
  female = c('Embryonic stem cells_Female_Ctrl', 'Embryonic stem cells_Untr_Ctrl'),
  male_female = c('Embryonic stem cells_Male_Ctrl', 'Embryonic stem cells_Female_Ctrl'),
  bapResponse_affectedByMale = c('Embryonic stem cells_Male_BaP', 'Embryonic stem cells_Untr_BaP'),
  bapResponse_affectedByFemale = c('Embryonic stem cells_Female_BaP', 'Embryonic stem cells_Untr_BaP'),
  bapResponse = c('Embryonic stem cells_Untr_BaP', 'Embryonic stem cells_Untr_Ctrl'),
  bapResponse_inMaleEnv = c('Embryonic stem cells_Male_BaP', 'Embryonic stem cells_Male_Ctrl'),
  bapResponse_inFemaleEnv = c('Embryonic stem cells_Female_BaP', 'Embryonic stem cells_Female_Ctrl')
)
go_e_results <- list(
  go_e_MC_vs_UC <- go_analysis(degs_epithelial_pct0.2[[1]], 1, 'BP'),
  go_e_FC_vs_UC <- go_analysis(degs_epithelial_pct0.2[[2]], 2, 'BP'),
  go_e_MC_vs_FC <- go_analysis(degs_epithelial_pct0.2[[3]], 3, 'BP'),
  go_e_MB_vs_UB <- go_analysis(degs_epithelial_pct0.2[[4]], 4, 'BP'),
  go_e_FB_vs_UB <- go_analysis(degs_epithelial_pct0.2[[5]], 5, 'BP'),
  go_e_UB_vs_UC <- go_analysis(degs_epithelial_pct0.2[[6]], 6, 'BP'),
  go_e_MB_vs_MC <- go_analysis(degs_epithelial_pct0.2[[7]], 7, 'BP'),
  go_e_FB_vs_FC <- go_analysis(degs_epithelial_pct0.2[[8]], 8, 'BP'))
plot_e_list <- list()
for (i in seq_along(go_e_results)) {
  plot <- dotplot(go_e_results[[i]],
                  x = "GeneRatio/BgRatio",
                  color = "qvalue",
                  showCategory = 20,
                  size = NULL,
                  title = "CC_dotplot") +
                  scale_color_gradient(low = "blue", high = "grey90", name = "qvalue") +
                  ggtitle(paste(titles[i], sep = "")) +  
                  theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
                  labs(x = "Fold enrichment")
  
  plot_e_list[[i]] <- plot
}
combined_e_plot <- do.call(grid.arrange, c(plot_e_list, ncol = 4))
ggsave('plots/go_e_foldenrichment.png', combined_e_plot, width = 25, height = 20, unit = 'in', dpi = 300)
plot_e_list <- list()
for (i in seq_along(go_e_results)) {
  plot <- dotplot(go_e_results[[i]],
                  x = "GeneRatio",
                  color = "qvalue",
                  showCategory = 20,
                  size = NULL,
                  title = "CC_dotplot") +
    scale_color_gradient(low = "blue", high = "grey90", name = "qvalue") +
    ggtitle(paste(titles[i], sep = "")) +  
    theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
    labs(x = "GeneRatio")
  
  plot_e_list[[i]] <- plot
}
combined_e_plot <- do.call(grid.arrange, c(plot_e_list, ncol = 4))
ggsave('plots/go_e_generatio.png', combined_e_plot, width = 25, height = 20, unit = 'in', dpi = 300)


################################## Fibroblasts #################################
degs_fibroblast_pct0.2 <- lapply(degs_fibroblast, function(df) {subset(df, pct.1 >= 0.2 & pct.2 >= 0.2)})
test_comparisons = list(
  male = c('Fibroblasts_Male_Ctrl', 'Fibroblasts_Untr_Ctrl'),
  female = c('Fibroblasts_Female_Ctrl', 'Fibroblasts_Untr_Ctrl'),
  male_female = c('Fibroblasts_Male_Ctrl', 'Fibroblasts_Female_Ctrl'),
  bapResponse_affectedByMale = c('Fibroblasts_Male_BaP', 'Fibroblasts_Untr_BaP'),
  bapResponse_affectedByFemale = c('Fibroblasts_Female_BaP', 'Fibroblasts_Untr_BaP'),
  bapResponse = c('Fibroblasts_Untr_BaP', 'Fibroblasts_Untr_Ctrl'),
  bapResponse_inMaleEnv = c('Fibroblasts_Male_BaP', 'Fibroblasts_Male_Ctrl'),
  bapResponse_inFemaleEnv = c('Fibroblasts_Female_BaP', 'Fibroblasts_Female_Ctrl')
)
go_f_results <- list(
  go_f_MC_vs_UC <- go_analysis(degs_fibroblast_pct0.2[[1]], 1),
  go_f_FC_vs_UC <- go_analysis(degs_fibroblast_pct0.2[[2]], 2),
  go_f_MC_vs_FC <- go_analysis(degs_fibroblast_pct0.2[[3]], 3),
  go_f_MB_vs_UB <- go_analysis(degs_fibroblast_pct0.2[[4]], 4),
  go_f_FB_vs_UB <- go_analysis(degs_fibroblast_pct0.2[[5]], 5),
  go_f_UB_vs_UC <- go_analysis(degs_fibroblast_pct0.2[[6]], 6),
  go_f_MB_vs_MC <- go_analysis(degs_fibroblast_pct0.2[[7]], 7),
  go_f_FB_vs_FC <- go_analysis(degs_fibroblast_pct0.2[[8]], 8))
plot_f_list <- list()
for (i in seq_along(go_f_results)) {
  plot <- dotplot(go_f_results[[i]],
                  x = "GeneRatio/BgRatio",
                  color = "qvalue",
                  showCategory = 20,
                  size = NULL,
                  title = "CC_dotplot") +
                  scale_color_gradient(low = "blue", high = "grey90", name = "qvalue") +
                  ggtitle(paste(titles[i], sep = "")) +  
                  theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
                  labs(x = "Fold enrichment")
  
  plot_f_list[[i]] <- plot
}
combined_f_plot <- do.call(grid.arrange, c(plot_f_list, ncol = 4))
ggsave('plots/go_f_foldenrichment.png', combined_f_plot, width = 25, height = 20, unit = 'in', dpi = 300)
plot_f_list <- list()
for (i in seq_along(go_f_results)) {
  plot <- dotplot(go_f_results[[i]],
                  x = "GeneRatio",
                  color = "qvalue",
                  showCategory = 20,
                  size = NULL,
                  title = "CC_dotplot") +
    scale_color_gradient(low = "blue", high = "grey90", name = "qvalue") +
    ggtitle(paste(titles[i], sep = "")) +  
    theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
    labs(x = "GeneRatio")
  
  plot_f_list[[i]] <- plot
}
combined_f_plot <- do.call(grid.arrange, c(plot_f_list, ncol = 4))
ggsave('plots/go_f_generatio.png', combined_f_plot, width = 25, height = 20, unit = 'in', dpi = 300)


################################## Cardiovascular cells #################################
degs_cardiovascular_pct0.2 <- lapply(degs_cardiovascular, function(df) {subset(df, pct.1 >= 0.2 & pct.2 >= 0.2)})
test_comparisons = list(
  male = c('Cardiovascular cells_Male_Ctrl', 'Cardiovascular cells_Untr_Ctrl'),
  female = c('Cardiovascular cells_Female_Ctrl', 'Cardiovascular cells_Untr_Ctrl'),
  male_female = c('Cardiovascular cells_Male_Ctrl', 'Cardiovascular cells_Female_Ctrl'),
  bapResponse_affectedByMale = c('Cardiovascular cells_Male_BaP', 'Cardiovascular cells_Untr_BaP'),
  bapResponse_affectedByFemale = c('Cardiovascular cells_Female_BaP', 'Cardiovascular cells_Untr_BaP'),
  bapResponse = c('Cardiovascular cells_Untr_BaP', 'Cardiovascular cells_Untr_Ctrl'),
  bapResponse_inMaleEnv = c('Cardiovascular cells_Male_BaP', 'Cardiovascular cells_Male_Ctrl'),
  bapResponse_inFemaleEnv = c('Cardiovascular cells_Female_BaP', 'Cardiovascular cells_Female_Ctrl')
)
go_c_results <- list(
  go_c_MC_vs_UC <- go_analysis(degs_fibroblast_pct0.2[[1]], 1),
  go_c_FC_vs_UC <- go_analysis(degs_fibroblast_pct0.2[[2]], 2),
  go_c_MC_vs_FC <- go_analysis(degs_fibroblast_pct0.2[[3]], 3),
  go_c_MB_vs_UB <- go_analysis(degs_fibroblast_pct0.2[[4]], 4),
  go_c_FB_vs_UB <- go_analysis(degs_fibroblast_pct0.2[[5]], 5),
  go_c_UB_vs_UC <- go_analysis(degs_fibroblast_pct0.2[[6]], 6),
  go_c_MB_vs_MC <- go_analysis(degs_fibroblast_pct0.2[[7]], 7),
  go_c_FB_vs_FC <- go_analysis(degs_fibroblast_pct0.2[[8]], 8))
plot_f_list <- list()
for (i in seq_along(go_f_results)) {
  plot <- dotplot(go_f_results[[i]],
                  x = "GeneRatio/BgRatio",
                  color = "qvalue",
                  showCategory = 20,
                  size = NULL,
                  title = "CC_dotplot") +
    scale_color_gradient(low = "blue", high = "grey90", name = "qvalue") +
    ggtitle(paste(titles[i], sep = "")) +  
    theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
    labs(x = "Fold enrichment")
  
  plot_f_list[[i]] <- plot
}
combined_f_plot <- do.call(grid.arrange, c(plot_f_list, ncol = 4))
ggsave('plots/go_c_foldenrichment.png', combined_f_plot, width = 25, height = 20, unit = 'in', dpi = 300)
plot_f_list <- list()
for (i in seq_along(go_f_results)) {
  plot <- dotplot(go_f_results[[i]],
                  x = "GeneRatio",
                  color = "qvalue",
                  showCategory = 20,
                  size = NULL,
                  title = "CC_dotplot") +
    scale_color_gradient(low = "blue", high = "grey90", name = "qvalue") +
    ggtitle(paste(titles[i], sep = "")) +  
    theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
    labs(x = "GeneRatio")
  
  plot_f_list[[i]] <- plot
}
combined_f_plot <- do.call(grid.arrange, c(plot_f_list, ncol = 4))
ggsave('plots/go_c_generatio.png', combined_f_plot, width = 25, height = 20, unit = 'in', dpi = 300)
```

---

## Visualisation

```r
#################### number of cells per type per sample #######################

table(sobj@meta.data$customclassif, sobj@meta.data$sample) %>% 
  as.data.frame() %>% 
  transform(percentage = Freq / tapply(Freq, Var2, sum)[Var2] * 100) %>% # find composition in percentage
  ggplot(aes(x = Var2, y = percentage, fill = Var1)) +
  geom_bar(stat = 'identity', position = 'fill', color = 'black') +
  theme_bw() +
  labs(x = 'Sample', y = 'Composition, %') + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_discrete(name = 'Cell type')


############ Effects of BAP on general differential expression##################
##### Cytoplasmic translation (GO:0006119) retrieved from amigo2
cytran <- read.table('data/cytoplasmic translation genes.txt', sep = '\t', header = F)
cytran_genes <- cytran$V2
##### Oxidative phosphorylation (GO:0006119) retrieved from amigo2
oxphos <- read.table('data/oxidative phosphorylation genes.txt', sep = '\t', header = F)
oxphos_genes <- oxphos$V2

############# Overview of DEGs caused by BaP
p1 <- ggplot(data = degs_airway_pct0.2$`Airway progenitor cells_Female_BaP_VS_Airway progenitor cells_Female_Ctrl`, 
             aes(x = avg_log2FC, y = -log10(p_val_adj), 
                 color = factor(ifelse(rownames(degs_airway_pct0.2$`Airway progenitor cells_Female_BaP_VS_Airway progenitor cells_Female_Ctrl`) %in% oxphos_genes, "oxphos", 
                                       ifelse(rownames(degs_airway_pct0.2$`Airway progenitor cells_Female_BaP_VS_Airway progenitor cells_Female_Ctrl`) %in% cytran_genes, "cytran", "other"))))) + 
      geom_point(size = 1) +
      theme_bw() +
      labs(x = "log2 (fold change)", y = "-log10 (q-value)", title = "Female BaP VS Female Ctrl") +
      geom_hline(yintercept = -log10(0.05), linetype = 'dotted', col = 'grey25') +
      geom_vline(xintercept = 0, linetype = 'solid', col = 'grey25') +
      xlim(c(-3.27, 3.27)) +
      theme(plot.title = element_text(hjust = 0.5)) +
      scale_color_manual(name = "GO pathway", values = c("oxphos" = "#41b349", "cytran" = "#fcc307", "other" = "black"))
p2 <- ggplot(data = degs_airway_pct0.2$`Airway progenitor cells_Male_BaP_VS_Airway progenitor cells_Male_Ctrl`,
             aes(x = avg_log2FC, y = -log10(p_val_adj), 
                 color = factor(ifelse(rownames(degs_airway_pct0.2$`Airway progenitor cells_Male_BaP_VS_Airway progenitor cells_Male_Ctrl`) %in% oxphos_genes, "oxphos", 
                                       ifelse(rownames(degs_airway_pct0.2$`Airway progenitor cells_Male_BaP_VS_Airway progenitor cells_Male_Ctrl`) %in% cytran_genes, "cytran", "other"))))) + 
      geom_point(size = 1) +
      theme_bw() +
      labs(x = "log2 (fold change)", y = "-log10 (q-value)", title = "Male BaP VS Male Ctrl") +
      geom_hline(yintercept = -log10(0.05), linetype = 'dotted', col = 'grey25') +
      geom_vline(xintercept = 0, linetype = 'solid', col = 'grey25') +
      xlim(c(-4.11, 4.11)) +
      theme(plot.title = element_text(hjust = 0.5)) +
      scale_color_manual(name = "GO pathway", values = c("oxphos" = "#41b349", "cytran" = "#fcc307", "other" = "black"))
p3 <- ggplot(data = degs_airway_pct0.2$`Airway progenitor cells_Untr_BaP_VS_Airway progenitor cells_Untr_Ctrl`, 
             aes(x = avg_log2FC, y = -log10(p_val_adj), 
                 color = factor(ifelse(rownames(degs_airway_pct0.2$`Airway progenitor cells_Untr_BaP_VS_Airway progenitor cells_Untr_Ctrl`) %in% oxphos_genes, "oxphos", 
                                       ifelse(rownames(degs_airway_pct0.2$`Airway progenitor cells_Untr_BaP_VS_Airway progenitor cells_Untr_Ctrl`) %in% cytran_genes, "cytran", "other"))))) + 
      geom_point(size = 1) +
      theme_bw() +
      labs(x = "log2 (fold change)", y = "-log10 (q-value)", title = "Untr BaP VS Untr Ctrl") +
      geom_hline(yintercept = -log10(0.05), linetype = 'dotted', col = 'grey25') +
      geom_vline(xintercept = 0, linetype = 'solid', col = 'grey25') +
      xlim(c(-1.81, 1.81)) +
      theme(plot.title = element_text(hjust = 0.5)) +
      scale_color_manual(name = "GO pathway", values = c("oxphos" = "#41b349", "cytran" = "#fcc307", "other" = "black"))
p4 <- ggplot(data = degs_epithelial_pct0.2$`Epithelial cells_Female_BaP_VS_Epithelial cells_Female_Ctrl`, 
             aes(x = avg_log2FC, y = -log10(p_val_adj), 
                 color = factor(ifelse(rownames(degs_epithelial_pct0.2$`Epithelial cells_Female_BaP_VS_Epithelial cells_Female_Ctrl`) %in% oxphos_genes, "oxphos", 
                                       ifelse(rownames(degs_epithelial_pct0.2$`Epithelial cells_Female_BaP_VS_Epithelial cells_Female_Ctrl`) %in% cytran_genes, "cytran", "other"))))) + 
      geom_point(size = 1) +
      theme_bw() +
      labs(x = "log2 (fold change)", y = "-log10 (q-value)", title = "Female BaP VS Female Ctrl") +
      geom_hline(yintercept = -log10(0.05), linetype = 'dotted', col = 'grey25') +
      geom_vline(xintercept = 0, linetype = 'solid', col = 'grey25') +
      xlim(c(-3.81, 3.81)) +
      theme(plot.title = element_text(hjust = 0.5)) +
  scale_color_manual(name = "GO pathway", values = c("oxphos" = "#41b349", "cytran" = "#fcc307", "other" = "black"))
p5 <- ggplot(data = degs_epithelial_pct0.2$`Epithelial cells_Male_BaP_VS_Epithelial cells_Male_Ctrl`, 
             aes(x = avg_log2FC, y = -log10(p_val_adj), 
                 color = factor(ifelse(rownames(degs_epithelial_pct0.2$`Epithelial cells_Male_BaP_VS_Epithelial cells_Male_Ctrl`) %in% oxphos_genes, "oxphos", 
                                       ifelse(rownames(degs_epithelial_pct0.2$`Epithelial cells_Male_BaP_VS_Epithelial cells_Male_Ctrl`) %in% cytran_genes, "cytran", "other"))))) + 
      geom_point(size = 1) +
      theme_bw() +
      labs(x = "log2 (fold change)", y = "-log10 (q-value)", title = "Male BaP VS Male Ctrl") +
      geom_hline(yintercept = -log10(0.05), linetype = 'dotted', col = 'grey25') +
      geom_vline(xintercept = 0, linetype = 'solid', col = 'grey25') +
      xlim(c(-2.82, 2.82)) +
      theme(plot.title = element_text(hjust = 0.5)) +
      scale_color_manual(name = "GO pathway", values = c("oxphos" = "#41b349", "cytran" = "#fcc307", "other" = "black"))
p6 <- ggplot(data = degs_epithelial_pct0.2$`Epithelial cells_Untr_BaP_VS_Epithelial cells_Untr_Ctrl`, 
             aes(x = avg_log2FC, y = -log10(p_val_adj), 
                 color = factor(ifelse(rownames(degs_epithelial_pct0.2$`Epithelial cells_Untr_BaP_VS_Epithelial cells_Untr_Ctrl`) %in% oxphos_genes, "oxphos", 
                                       ifelse(rownames(degs_epithelial_pct0.2$`Epithelial cells_Untr_BaP_VS_Epithelial cells_Untr_Ctrl`) %in% cytran_genes, "cytran", "other"))))) + 
      geom_point(size = 1) +
      theme_bw() +
      labs(x = "log2 (fold change)", y = "-log10 (q-value)", title = "Untr BaP VS Untr Ctrl") +
      geom_hline(yintercept = -log10(0.05), linetype = 'dotted', col = 'grey25') +
      geom_vline(xintercept = 0, linetype = 'solid', col = 'grey25') +
      xlim(c(-2.42, 2.42)) +
      theme(plot.title = element_text(hjust = 0.5)) +
      scale_color_manual(name = "GO pathway", values = c("oxphos" = "#41b349", "cytran" = "#fcc307", "other" = "black"))
p7 <- ggplot(data = degs_fibroblast_pct0.2$Fibroblasts_Female_BaP_VS_Fibroblasts_Female_Ctrl, 
             aes(x = avg_log2FC, y = -log10(p_val_adj), 
                 color = factor(ifelse(rownames(degs_fibroblast_pct0.2$Fibroblasts_Female_BaP_VS_Fibroblasts_Female_Ctrl) %in% oxphos_genes, "oxphos", 
                                       ifelse(rownames(degs_fibroblast_pct0.2$Fibroblasts_Female_BaP_VS_Fibroblasts_Female_Ctrl) %in% cytran_genes, "cytran", "other"))))) + 
      geom_point(size = 1) +
      theme_bw() +
      labs(x = "log2 (fold change)", y = "-log10 (q-value)", title = "Female BaP VS Female Ctrl") +
      geom_hline(yintercept = -log10(0.05), linetype = 'dotted', col = 'grey25') +
      geom_vline(xintercept = 0, linetype = 'solid', col = 'grey25') +
      xlim(c(-1.69, 1.69)) +
      theme(plot.title = element_text(hjust = 0.5)) +
      scale_color_manual(name = "GO pathway", values = c("oxphos" = "#41b349", "cytran" = "#fcc307", "other" = "black"))
p8 <- ggplot(data = degs_fibroblast_pct0.2$Fibroblasts_Male_BaP_VS_Fibroblasts_Male_Ctrl, 
             aes(x = avg_log2FC, y = -log10(p_val_adj), 
                 color = factor(ifelse(rownames(degs_fibroblast_pct0.2$Fibroblasts_Male_BaP_VS_Fibroblasts_Male_Ctrl) %in% oxphos_genes, "oxphos", 
                                       ifelse(rownames(degs_fibroblast_pct0.2$Fibroblasts_Male_BaP_VS_Fibroblasts_Male_Ctrl) %in% cytran_genes, "cytran", "other"))))) + 
      geom_point(size = 1) +
      theme_bw() +
      labs(x = "log2 (fold change)", y = "-log10 (q-value)", title = "Male BaP VS Male Ctrl") +
      geom_hline(yintercept = -log10(0.05), linetype = 'dotted', col = 'grey25') +
      geom_vline(xintercept = 0, linetype = 'solid', col = 'grey25') +
      xlim(c(-2.39, 2.39)) +
      theme(plot.title = element_text(hjust = 0.5)) +
      scale_color_manual(name = "GO pathway", values = c("oxphos" = "#41b349", "cytran" = "#fcc307", "other" = "black"))
p9 <- ggplot(data = degs_fibroblast_pct0.2$Fibroblasts_Untr_BaP_VS_Fibroblasts_Untr_Ctrl, 
             aes(x = avg_log2FC, y = -log10(p_val_adj), 
                 color = factor(ifelse(rownames(degs_fibroblast_pct0.2$Fibroblasts_Untr_BaP_VS_Fibroblasts_Untr_Ctrl) %in% oxphos_genes, "oxphos", 
                                       ifelse(rownames(degs_fibroblast_pct0.2$Fibroblasts_Untr_BaP_VS_Fibroblasts_Untr_Ctrl) %in% cytran_genes, "cytran", "other"))))) + 
      geom_point(size = 1) +
      theme_bw() +
      labs(x = "log2 (fold change)", y = "-log10 (q-value)", title = "Untr BaP VS Untr Ctrl") +
      geom_hline(yintercept = -log10(0.05), linetype = 'dotted', col = 'grey25') +
      geom_vline(xintercept = 0, linetype = 'solid', col = 'grey25') +
      xlim(c(-1.84, 1.84)) +
      theme(plot.title = element_text(hjust = 0.5)) +
  scale_color_manual(name = "GO pathway", values = c("oxphos" = "#41b349", "cytran" = "#fcc307", "other" = "black"))
patchwork::wrap_plots(p1, p2, p3, p4, p5, p6, p7, p8, p9, ncol = 3)
ggsave('plots/general_degs_vol.png', width = 36, height = 30, units = 'cm', dpi = 300)


########## Effects of BAP on genetic expression of AHR target genes ############
ahr_genes <- c('AHR', 'AHRR', 'TIPARP', 'NQO1', 'ALDH3A1', 'ALDH3A2', 'ALDH3B1', 'ALDH3B2', 'CYP1A1', 'CYP1A2', 'CYP1B1')
expression_data <- FetchData(object = sobj, vars = ahr_genes)
average_expression <- colMeans(expression_data)
average_expression_df <- data.frame(Gene = names(average_expression), AverageExpression = average_expression)
sorted_expression <- average_expression_df %>%
  arrange(desc(AverageExpression))
print(sorted_expression)

FeaturePlot(sobj, reduction = 'umap', features = ahr_genes, ncol = 4) 
ggsave('plots/ahr_umap.png', width = 18, height = 12, unit = 'in', dpi = 300)

# Airway
for (dataset_name in names(degs_airway_pct0.2)) {
  gene_names <- rownames(degs_airway_pct0.2[[dataset_name]])
  matched_genes <- grep(paste(ahr_genes, collapse = "|"), gene_names, value = TRUE)
  cat("Matched genes in", dataset_name, ":", "\n")
  print(matched_genes)
}

# Epithelial
for (dataset_name in names(degs_epithelial_pct0.2)) {
  gene_names <- rownames(degs_epithelial_pct0.2[[dataset_name]])
  matched_genes <- grep(paste(ahr_genes, collapse = "|"), gene_names, value = TRUE)
  cat("Matched genes in", dataset_name, ":", "\n")
  print(matched_genes)
}

# Fibroblasts
for (dataset_name in names(degs_fibroblast_pct0.2)) {
  gene_names <- rownames(degs_fibroblast_pct0.2[[dataset_name]])
  matched_genes <- grep(paste(ahr_genes, collapse = "|"), gene_names, value = TRUE)
  cat("Matched genes in", dataset_name, ":", "\n")
  print(matched_genes)
}


############# DEGs on AHR genes caused by BaP
p1 <- ggplot(data = degs_airway_pct0.2$`Airway progenitor cells_Female_BaP_VS_Airway progenitor cells_Female_Ctrl`, 
             aes(x = avg_log2FC, y = -log10(p_val_adj), 
                 color = factor(ifelse(rownames(degs_airway_pct0.2$`Airway progenitor cells_Female_BaP_VS_Airway progenitor cells_Female_Ctrl`) %in% ahr_genes, "ahr", "other")))) + 
  geom_point(size = 1) +
  geom_shadowtext(data = degs_airway_pct0.2$`Airway progenitor cells_Female_BaP_VS_Airway progenitor cells_Female_Ctrl`[rownames(degs_airway_pct0.2$`Airway progenitor cells_Female_BaP_VS_Airway progenitor cells_Female_Ctrl`) %in% ahr_genes, ],
            aes(label = rownames(degs_airway_pct0.2$`Airway progenitor cells_Female_BaP_VS_Airway progenitor cells_Female_Ctrl`[rownames(degs_airway_pct0.2$`Airway progenitor cells_Female_BaP_VS_Airway progenitor cells_Female_Ctrl`) %in% ahr_genes, ])), 
            vjust = -0.5, hjust = 0.5, color = "#fcc307", size = 4, fontface = "bold") + 
  theme_bw() +
  labs(x = "log2 (fold change)", y = "-log10 (q-value)", title = "Female BaP VS Female Ctrl") +
  geom_hline(yintercept = -log10(0.05), linetype = 'dotted', col = 'grey25') +
  geom_vline(xintercept = 0, linetype = 'solid', col = 'grey25') +
  xlim(c(-3.69, 3.69)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_color_manual(name = "Cell types", values = c("ahr" = "#fcc307", "other" = "black"))
p2 <- ggplot(data = degs_airway_pct0.2$`Airway progenitor cells_Male_BaP_VS_Airway progenitor cells_Male_Ctrl`,
             aes(x = avg_log2FC, y = -log10(p_val_adj), 
                 color = factor(ifelse(rownames(degs_airway_pct0.2$`Airway progenitor cells_Male_BaP_VS_Airway progenitor cells_Male_Ctrl`) %in% ahr_genes, "ahr", "other")))) +
  geom_point(size = 1) +
  geom_shadowtext(data = degs_airway_pct0.2$`Airway progenitor cells_Male_BaP_VS_Airway progenitor cells_Male_Ctrl`[rownames(degs_airway_pct0.2$`Airway progenitor cells_Male_BaP_VS_Airway progenitor cells_Male_Ctrl`) %in% ahr_genes, ],
            aes(label = rownames(degs_airway_pct0.2$`Airway progenitor cells_Male_BaP_VS_Airway progenitor cells_Male_Ctrl`[rownames(degs_airway_pct0.2$`Airway progenitor cells_Male_BaP_VS_Airway progenitor cells_Male_Ctrl`) %in% ahr_genes, ])), 
            vjust = -0.5, hjust = 0.5, color = "#fcc307", size = 4, fontface = "bold") + 
  theme_bw() +
  labs(x = "log2 (fold change)", y = "-log10 (q-value)", title = "Male BaP VS Male Ctrl") +
  geom_hline(yintercept = -log10(0.05), linetype = 'dotted', col = 'grey25') +
  geom_vline(xintercept = 0, linetype = 'solid', col = 'grey25') +
  xlim(c(-4.11, 4.11)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_color_manual(name = "Cell types", values = c("ahr" = "#fcc307", "other" = "black"))
p3 <- ggplot(data = degs_airway_pct0.2$`Airway progenitor cells_Untr_BaP_VS_Airway progenitor cells_Untr_Ctrl`, 
             aes(x = avg_log2FC, y = -log10(p_val_adj), 
                 color = factor(ifelse(rownames(degs_airway_pct0.2$`Airway progenitor cells_Untr_BaP_VS_Airway progenitor cells_Untr_Ctrl`) %in% ahr_genes, "ahr", "other")))) +
  geom_point(size = 1) +
  geom_shadowtext(data = degs_airway_pct0.2$`Airway progenitor cells_Untr_BaP_VS_Airway progenitor cells_Untr_Ctrl`[rownames(degs_airway_pct0.2$`Airway progenitor cells_Untr_BaP_VS_Airway progenitor cells_Untr_Ctrl`) %in% ahr_genes, ],
            aes(label = rownames(degs_airway_pct0.2$`Airway progenitor cells_Untr_BaP_VS_Airway progenitor cells_Untr_Ctrl`[rownames(degs_airway_pct0.2$`Airway progenitor cells_Untr_BaP_VS_Airway progenitor cells_Untr_Ctrl`) %in% ahr_genes, ])), 
            vjust = -0.5, hjust = 0.5, color = "#fcc307", size = 4, fontface = "bold") + 
  theme_bw() +
  labs(x = "log2 (fold change)", y = "-log10 (q-value)", title = "Untr BaP VS Untr Ctrl") +
  geom_hline(yintercept = -log10(0.05), linetype = 'dotted', col = 'grey25') +
  geom_vline(xintercept = 0, linetype = 'solid', col = 'grey25') +
  xlim(c(-1.81, 1.81)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_color_manual(name = "Cell types", values = c("ahr" = "#fcc307", "other" = "black"))
p4 <- ggplot(data = degs_epithelial_pct0.2$`Epithelial cells_Female_BaP_VS_Epithelial cells_Female_Ctrl`, 
             aes(x = avg_log2FC, y = -log10(p_val_adj), 
                 color = factor(ifelse(rownames(degs_epithelial_pct0.2$`Epithelial cells_Female_BaP_VS_Epithelial cells_Female_Ctrl`) %in% ahr_genes, "ahr", "other")))) +
  geom_point(size = 1) +
  geom_shadowtext(data = degs_epithelial_pct0.2$`Epithelial cells_Female_BaP_VS_Epithelial cells_Female_Ctrl`[rownames(degs_epithelial_pct0.2$`Epithelial cells_Female_BaP_VS_Epithelial cells_Female_Ctrl`) %in% ahr_genes, ],
            aes(label = rownames(degs_epithelial_pct0.2$`Epithelial cells_Female_BaP_VS_Epithelial cells_Female_Ctrl`[rownames(degs_epithelial_pct0.2$`Epithelial cells_Female_BaP_VS_Epithelial cells_Female_Ctrl`) %in% ahr_genes, ])), 
            vjust = -0.5, hjust = 0.5, color = "#fcc307", size = 4, fontface = "bold") +
  theme_bw() +
  labs(x = "log2 (fold change)", y = "-log10 (q-value)", title = "Female BaP VS Female Ctrl") +
  geom_hline(yintercept = -log10(0.05), linetype = 'dotted', col = 'grey25') +
  geom_vline(xintercept = 0, linetype = 'solid', col = 'grey25') +
  xlim(c(-3.81, 3.81)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_color_manual(name = "Cell types", values = c("ahr" = "#fcc307", "other" = "black"))
p5 <- ggplot(data = degs_epithelial_pct0.2$`Epithelial cells_Male_BaP_VS_Epithelial cells_Male_Ctrl`, 
             aes(x = avg_log2FC, y = -log10(p_val_adj), 
                 color = factor(ifelse(rownames(degs_epithelial_pct0.2$`Epithelial cells_Male_BaP_VS_Epithelial cells_Male_Ctrl`) %in% ahr_genes, "ahr", "other")))) +
  geom_point(size = 1) +
  geom_shadowtext(data = degs_epithelial_pct0.2$`Epithelial cells_Male_BaP_VS_Epithelial cells_Male_Ctrl`[rownames(degs_epithelial_pct0.2$`Epithelial cells_Male_BaP_VS_Epithelial cells_Male_Ctrl`) %in% ahr_genes, ],
            aes(label = rownames(degs_epithelial_pct0.2$`Epithelial cells_Male_BaP_VS_Epithelial cells_Male_Ctrl`[rownames(degs_epithelial_pct0.2$`Epithelial cells_Male_BaP_VS_Epithelial cells_Male_Ctrl`) %in% ahr_genes, ])), 
            vjust = -0.5, hjust = 0.5, color = "#fcc307", size = 4, fontface = "bold") +
  theme_bw() +
  labs(x = "log2 (fold change)", y = "-log10 (q-value)", title = "Male BaP VS Male Ctrl") +
  geom_hline(yintercept = -log10(0.05), linetype = 'dotted', col = 'grey25') +
  geom_vline(xintercept = 0, linetype = 'solid', col = 'grey25') +
  xlim(c(-2.82, 2.82)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_color_manual(name = "Cell types", values = c("ahr" = "#fcc307", "other" = "black"))
p6 <- ggplot(data = degs_epithelial_pct0.2$`Epithelial cells_Untr_BaP_VS_Epithelial cells_Untr_Ctrl`, 
             aes(x = avg_log2FC, y = -log10(p_val_adj), 
                 color = factor(ifelse(rownames(degs_epithelial_pct0.2$`Epithelial cells_Untr_BaP_VS_Epithelial cells_Untr_Ctrl`) %in% ahr_genes, "ahr", "other")))) +
  geom_point(size = 1) +
  geom_shadowtext(data = degs_epithelial_pct0.2$`Epithelial cells_Untr_BaP_VS_Epithelial cells_Untr_Ctrl`[rownames(degs_epithelial_pct0.2$`Epithelial cells_Untr_BaP_VS_Epithelial cells_Untr_Ctrl`) %in% ahr_genes, ],
            aes(label = rownames(degs_epithelial_pct0.2$`Epithelial cells_Untr_BaP_VS_Epithelial cells_Untr_Ctrl`[rownames(degs_epithelial_pct0.2$`Epithelial cells_Untr_BaP_VS_Epithelial cells_Untr_Ctrl`) %in% ahr_genes, ])), 
            vjust = -0.5, hjust = 0.5, color = "#fcc307", size = 4, fontface = "bold") +
  theme_bw() +
  labs(x = "log2 (fold change)", y = "-log10 (q-value)", title = "Untr BaP VS Untr Ctrl") +
  geom_hline(yintercept = -log10(0.05), linetype = 'dotted', col = 'grey25') +
  geom_vline(xintercept = 0, linetype = 'solid', col = 'grey25') +
  xlim(c(-2.42, 2.42)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_color_manual(name = "Cell types", values = c("ahr" = "#fcc307", "other" = "black"))
p7 <- ggplot(data = degs_fibroblast_pct0.2$Fibroblasts_Female_BaP_VS_Fibroblasts_Female_Ctrl, 
             aes(x = avg_log2FC, y = -log10(p_val_adj), 
                 color = factor(ifelse(rownames(degs_fibroblast_pct0.2$Fibroblasts_Female_BaP_VS_Fibroblasts_Female_Ctrl) %in% ahr_genes, "ahr", "other")))) +
  geom_point(size = 1) +
  geom_shadowtext(data = degs_fibroblast_pct0.2$Fibroblasts_Female_BaP_VS_Fibroblasts_Female_Ctrl[rownames(degs_fibroblast_pct0.2$Fibroblasts_Female_BaP_VS_Fibroblasts_Female_Ctrl) %in% ahr_genes, ],
                  aes(label = rownames(degs_fibroblast_pct0.2$Fibroblasts_Female_BaP_VS_Fibroblasts_Female_Ctrl[rownames(degs_fibroblast_pct0.2$Fibroblasts_Female_BaP_VS_Fibroblasts_Female_Ctrl) %in% ahr_genes, ])), 
                  vjust = -0.5, hjust = 0.5, color = "#fcc307", size = 4, fontface = "bold", bg.color = "black") + 
  theme_bw() +
  labs(x = "log2 (fold change)", y = "-log10 (q-value)", title = "Female BaP VS Female Ctrl") +
  geom_hline(yintercept = -log10(0.05), linetype = 'dotted', col = 'grey25') +
  geom_vline(xintercept = 0, linetype = 'solid', col = 'grey25') +
  xlim(c(-1.69, 1.69)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_color_manual(name = "Cell types", values = c("ahr" = "#fcc307", "other" = "black"))
p7 <- ggplot(data = degs_fibroblast_pct0.2$Fibroblasts_Female_BaP_VS_Fibroblasts_Female_Ctrl, 
             aes(x = avg_log2FC, y = -log10(p_val_adj), 
                 color = factor(ifelse(rownames(degs_fibroblast_pct0.2$Fibroblasts_Female_BaP_VS_Fibroblasts_Female_Ctrl) %in% ahr_genes, "ahr", "other")))) +
  geom_point(size = 1) +
  geom_text(data = degs_fibroblast_pct0.2$Fibroblasts_Female_BaP_VS_Fibroblasts_Female_Ctrl[rownames(degs_fibroblast_pct0.2$Fibroblasts_Female_BaP_VS_Fibroblasts_Female_Ctrl) %in% ahr_genes, ],
            aes(label = rownames(degs_fibroblast_pct0.2$Fibroblasts_Female_BaP_VS_Fibroblasts_Female_Ctrl[rownames(degs_fibroblast_pct0.2$Fibroblasts_Female_BaP_VS_Fibroblasts_Female_Ctrl) %in% ahr_genes, ])), 
            vjust = -0.5, hjust = 0.5, color = "#fcc307", size = 4, fontface = "bold") +
  theme_bw() +
  labs(x = "log2 (fold change)", y = "-log10 (q-value)", title = "Female BaP VS Female Ctrl") +
  geom_hline(yintercept = -log10(0.05), linetype = 'dotted', col = 'grey25') +
  geom_vline(xintercept = 0, linetype = 'solid', col = 'grey25') +
  xlim(c(-1.69, 1.69)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_color_manual(name = "Cell types", values = c("ahr" = "#fcc307", "other" = "black"))
p8 <- ggplot(data = degs_fibroblast_pct0.2$Fibroblasts_Male_BaP_VS_Fibroblasts_Male_Ctrl, 
             aes(x = avg_log2FC, y = -log10(p_val_adj), 
                 color = factor(ifelse(rownames(degs_fibroblast_pct0.2$Fibroblasts_Male_BaP_VS_Fibroblasts_Male_Ctrl) %in% ahr_genes, "ahr", "other")))) +
  geom_point(size = 1) +
  geom_shadowtext(data = degs_fibroblast_pct0.2$Fibroblasts_Male_BaP_VS_Fibroblasts_Male_Ctrl[rownames(degs_fibroblast_pct0.2$Fibroblasts_Male_BaP_VS_Fibroblasts_Male_Ctrl) %in% ahr_genes, ],
            aes(label = rownames(degs_fibroblast_pct0.2$Fibroblasts_Male_BaP_VS_Fibroblasts_Male_Ctrl[rownames(degs_fibroblast_pct0.2$Fibroblasts_Male_BaP_VS_Fibroblasts_Male_Ctrl) %in% ahr_genes, ])), 
            vjust = -0.5, hjust = 0.5, color = "#fcc307", size = 4, fontface = "bold") +
  theme_bw() +
  labs(x = "log2 (fold change)", y = "-log10 (q-value)", title = "Male BaP VS Male Ctrl") +
  geom_hline(yintercept = -log10(0.05), linetype = 'dotted', col = 'grey25') +
  geom_vline(xintercept = 0, linetype = 'solid', col = 'grey25') +
  xlim(c(-2.39, 2.39)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_color_manual(name = "Cell types", values = c("ahr" = "#fcc307", "other" = "black"))
p9 <- ggplot(data = degs_fibroblast_pct0.2$Fibroblasts_Untr_BaP_VS_Fibroblasts_Untr_Ctrl, 
             aes(x = avg_log2FC, y = -log10(p_val_adj), 
                 color = factor(ifelse(rownames(degs_fibroblast_pct0.2$Fibroblasts_Untr_BaP_VS_Fibroblasts_Untr_Ctrl) %in% ahr_genes, "ahr", "other")))) +
  geom_point(size = 1) +
  geom_shadowtext(data = degs_fibroblast_pct0.2$Fibroblasts_Untr_BaP_VS_Fibroblasts_Untr_Ctrl[rownames(degs_fibroblast_pct0.2$Fibroblasts_Untr_BaP_VS_Fibroblasts_Untr_Ctrl) %in% ahr_genes, ],
            aes(label = rownames(degs_fibroblast_pct0.2$Fibroblasts_Untr_BaP_VS_Fibroblasts_Untr_Ctrl[rownames(degs_fibroblast_pct0.2$Fibroblasts_Untr_BaP_VS_Fibroblasts_Untr_Ctrl) %in% ahr_genes, ])), 
            vjust = -0.5, hjust = 0.5, color = "#fcc307", size = 4, fontface = "bold") +
  theme_bw() +
  labs(x = "log2 (fold change)", y = "-log10 (q-value)", title = "Untr BaP VS Untr Ctrl") +
  geom_hline(yintercept = -log10(0.05), linetype = 'dotted', col = 'grey25') +
  geom_vline(xintercept = 0, linetype = 'solid', col = 'grey25') +
  xlim(c(-1.84, 1.84)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_color_manual(name = "Cell types", values = c("ahr" = "#fcc307", "other" = "black"))
patchwork::wrap_plots(p1, p2, p3, p4, p5, p6, p7, p8, p9, ncol = 3)
ggsave('plots/ahr_degs_vol.png', width = 36, height = 30, units = 'cm', dpi = 300)


####### Effects of BAP on genetic expression of cell-type specific markers #####

# Airway
airway_genes = as.factor(c('ABI3BP','AQP3','DAPL1','GSTM2','HPGD','ICAM1','KRT14','KRT15','KRT5','PHLDA3','RPS18','SDC1'))
for (dataset_name in names(degs_airway_pct0.2)) {
  dataset <- degs_airway_pct0.2[[dataset_name]]
  gene_names <- rownames(dataset)
  matched_genes <- gene_names[gene_names %in% airway_genes]
  
  cat("Matched genes in", dataset_name, ":\n")
  
  for (gene in matched_genes) {
    log2fc <- dataset[gene, "avg_log2FC"]
    regulation <- ifelse(log2fc > 0, "Upregulated", "Downregulated")
    print(paste(gene, ": log2FC =", log2fc, "(", regulation, ")"))
  }
  cat("\n")
}


# Epithelial
epithelial_genes = as.factor(c('CD24','CEACAM1','ST6GAL1','ITGB4','IL1R1','PROM1','CDH1','KRT1','KRT7',
          'MUC1','ICAM1','KRT14','KRT5','ITGAL','CD2','KLK3','ITGA5','ITGA4','ITGA2',
          'KRT3','KRT16','SCNN1A','KRT15','ITGA1','KRT2','SCNN1B','SCNN1D','SCNN1G',
          'IFI16','BOK','NKD1','FZD6','DKK3','NRP2','SFRP5','RAI14','DEFB1','KLK1',
          'AGR2','APOA1','GPA33','ANPEP','CRYBA1','BMI1','BRCA1','MUC16','CEACAM5',
          'CTSE','SCGB1A1','EXO1','FOXA1','GABRP','GGT1','SFN','KRT13','LTF','SLC46A2',
          'KLK10','P2RX7','CDKN2A','TP63','CDH3','PSCA','AGER','ZFP42','SPRR1B','SI',
          'TTR','TM4SF20','TSTD1','SYCN','HBEGF','PIGR','MUC13','SELENBP1','ELF3',
          'TSPAN1','GUCA2A','PHGR1','LYPD8','LGALS4','GATA2','SEC23B','TSPAN8','DLX5',
          'DGAT2','ITPR2','THRSP','PLA2G4A','SLC25A48','PGR','FERMT1','EHF','PLEKHS1',
          'CDKL1','MECOM','MSX1','RNF128','ANLN','CKAP2','HMMR','KIF15','CKAP2L','KIF20B',
          'HIRIP3','INCENP','KIF23','PRC1','ECT2','CXCL10','CXCL8','CCL20','CXCL17','PRG4',
          'ALOX15','F5','EMILIN2','SPTSSB','FMO5','IVL','VSIG2','AQP3','PAQR5','EPCAM',
          'CLDN1','OCLN','MUC5AC'))
for (dataset_name in names(degs_epithelial_pct0.2)) {
  dataset <- degs_epithelial_pct0.2[[dataset_name]]
  gene_names <- rownames(dataset)
  matched_genes <- gene_names[gene_names %in% epithelial_genes]
  
  cat("Matched genes in", dataset_name, ":\n")
  
  for (gene in matched_genes) {
    log2fc <- dataset[gene, "avg_log2FC"]
    regulation <- ifelse(log2fc > 0, "Upregulated", "Downregulated")
    print(paste(gene, ": log2FC =", log2fc, "(", regulation, ")"))
  }
  cat("\n")
}


# Fibroblasts
fibroblast_genes <- as.factor(c('IL1R1','FAP','FLI1','CELA1','LOX','PDGFRB','P4HA1','UCP2','CCR2','ITGAL',
          'FGR','HCK','TNFRSF1B','PRKCD','ENO3','ABI3','TREML4','PIP4K2A','CD300E',
          'SERPINB10','CTHRC1','TBX18','COL15A1','GJB2','IL34','EDN3','SLC6A13','VTN',
          'ITIH5','LUM','DPT','POSTN','PENK','MMP14','COL6A2','FABP4','ASPN','ANGPTL2',
          'EFEMP1','SCARA5','IGFBP3','COPZ2','DPEP1','ADAMTS5','COL5A1','CD248','PI16',
          'PAMR1','TNXB','MMP2','COL14A1','CLEC3B','IGFBP6','COL5A2','FBN1','MFAP5',
          'FKBP10','PALLD','WIF1','SNHG18','CDH11','PTCH1','ARAP1','FBLN2','IGF1','PRRX1',
          'FKBP7','OAF','COL6A3','CTSK','DKK1','C1S','RARRES2','GREM1','SPON2','TCF21',
          'PCSK6','COL8A1','ENTPD2','CXCL8','CXCL3','IL6','CYP1B1','COL13A1','ADAMTS10',
          'CCL11','ADAM33','COL4A3','COL4A4','LAMA2','ACKR3','CD55','FBLN7','FIBIN',
          'THBS2','NOV','PTX3','MMP3','LRRK1','HGF','FRZB','COL12A1','COL7A1','MEOX1',
          'PRG4','PKD2','CCL19','NNMT','FOXF1','HAS1','CTGF','ERCC1','WISP1','TWIST2',
          'RIPK3','DDR2','ELN','FN1','HHIP','FMO2','COL1A2','COL3A1','VIM','FSTL1','GSN',
          'SPARC','S100A4','NT5E','COL1A1','MGP','NOX4','THY1','CD40','SERPINH1','CD44',
          'PDGFRA','EN1','DCN','CEBPB','EGR1','FOSB','FOSL2','HIF1A','KLF2','KLF4','KLF6',
          'KLF9','NFAT5','NFATC1','NFKB1','NR4A1','NR4A2','PBX1','RUNX1','STAT3','TCF4',
          'ZEB2','LAMC1','MEDAG','LAMB1','DKK3','TBX20','MDK','GSTM5','NGF','VEGFA','FGF2',
          'P4HTM','CKAP4','INMT','CXCL14'))
for (dataset_name in names(degs_fibroblast_pct0.2)) {
  dataset <- degs_fibroblast_pct0.2[[dataset_name]]
  gene_names <- rownames(dataset)
  matched_genes <- gene_names[gene_names %in% fibroblast_genes]
  cat("Matched genes in", dataset_name, ":\n")
  for (gene in matched_genes) {
    log2fc <- dataset[gene, "avg_log2FC"]
    regulation <- ifelse(log2fc > 0, "Upregulated", "Downregulated")
    cat(gene, ": log2FC =", log2fc, "(", regulation, ")\n")
  }
  cat("\n")
}


############# DEGs on markers caused by BaP
p1 <- ggplot(data = degs_airway_pct0.2$`Airway progenitor cells_Female_BaP_VS_Airway progenitor cells_Female_Ctrl`, 
             aes(x = avg_log2FC, y = -log10(p_val_adj), 
                 color = factor(ifelse(rownames(degs_airway_pct0.2$`Airway progenitor cells_Female_BaP_VS_Airway progenitor cells_Female_Ctrl`) %in% airway_genes, "airway", "other")))) + 
  geom_point(size = 1) +
  theme_bw() +
  labs(x = "log2 (fold change)", y = "-log10 (q-value)", title = "Female BaP VS Female Ctrl") +
  geom_hline(yintercept = -log10(0.05), linetype = 'dotted', col = 'grey25') +
  geom_vline(xintercept = 0, linetype = 'solid', col = 'grey25') +
  xlim(c(-3.27, 3.27)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_color_manual(name = "Cell types", values = c("airway" = "#2775b6", "other" = "black"))
p2 <- ggplot(data = degs_airway_pct0.2$`Airway progenitor cells_Male_BaP_VS_Airway progenitor cells_Male_Ctrl`,
             aes(x = avg_log2FC, y = -log10(p_val_adj), 
                 color = factor(ifelse(rownames(degs_airway_pct0.2$`Airway progenitor cells_Male_BaP_VS_Airway progenitor cells_Male_Ctrl`) %in% airway_genes, "airway", "other")))) + 
  geom_point(size = 1) +
  theme_bw() +
  labs(x = "log2 (fold change)", y = "-log10 (q-value)", title = "Male BaP VS Male Ctrl") +
  geom_hline(yintercept = -log10(0.05), linetype = 'dotted', col = 'grey25') +
  geom_vline(xintercept = 0, linetype = 'solid', col = 'grey25') +
  xlim(c(-4.11, 4.11)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_color_manual(name = "Cell types", values = c("airway" = "#2775b6", "other" = "black"))
p3 <- ggplot(data = degs_airway_pct0.2$`Airway progenitor cells_Untr_BaP_VS_Airway progenitor cells_Untr_Ctrl`, 
             aes(x = avg_log2FC, y = -log10(p_val_adj), 
                 color = factor(ifelse(rownames(degs_airway_pct0.2$`Airway progenitor cells_Untr_BaP_VS_Airway progenitor cells_Untr_Ctrl`) %in% airway_genes, "airway", "other")))) + 
  geom_point(size = 1) +
  theme_bw() +
  labs(x = "log2 (fold change)", y = "-log10 (q-value)", title = "Untr BaP VS Untr Ctrl") +
  geom_hline(yintercept = -log10(0.05), linetype = 'dotted', col = 'grey25') +
  geom_vline(xintercept = 0, linetype = 'solid', col = 'grey25') +
  xlim(c(-1.81, 1.81)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_color_manual(name = "Cell types", values = c("airway" = "#2775b6", "other" = "black"))
p4 <- ggplot(data = degs_epithelial_pct0.2$`Epithelial cells_Female_BaP_VS_Epithelial cells_Female_Ctrl`, 
             aes(x = avg_log2FC, y = -log10(p_val_adj), 
                 color = factor(ifelse(rownames(degs_epithelial_pct0.2$`Epithelial cells_Female_BaP_VS_Epithelial cells_Female_Ctrl`) %in% epithelial_genes, "epithelial","other")))) + 
  geom_point(size = 1) +
  theme_bw() +
  labs(x = "log2 (fold change)", y = "-log10 (q-value)", title = "Female BaP VS Female Ctrl") +
  geom_hline(yintercept = -log10(0.05), linetype = 'dotted', col = 'grey25') +
  geom_vline(xintercept = 0, linetype = 'solid', col = 'grey25') +
  xlim(c(-3.81, 3.81)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_color_manual(name = "Cell types", values = c("epithelial" = "#41b349", "other" = "black"))
p5 <- ggplot(data = degs_epithelial_pct0.2$`Epithelial cells_Male_BaP_VS_Epithelial cells_Male_Ctrl`, 
             aes(x = avg_log2FC, y = -log10(p_val_adj), 
                 color = factor(ifelse(rownames(degs_epithelial_pct0.2$`Epithelial cells_Male_BaP_VS_Epithelial cells_Male_Ctrl`) %in% epithelial_genes, "epithelial","other")))) + 
  geom_point(size = 1) +
  theme_bw() +
  labs(x = "log2 (fold change)", y = "-log10 (q-value)", title = "Male BaP VS Male Ctrl") +
  geom_hline(yintercept = -log10(0.05), linetype = 'dotted', col = 'grey25') +
  geom_vline(xintercept = 0, linetype = 'solid', col = 'grey25') +
  xlim(c(-2.82, 2.82)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_color_manual(name = "Cell types", values = c("epithelial" = "#41b349", "other" = "black"))
p6 <- ggplot(data = degs_epithelial_pct0.2$`Epithelial cells_Untr_BaP_VS_Epithelial cells_Untr_Ctrl`, 
             aes(x = avg_log2FC, y = -log10(p_val_adj), 
                 color = factor(ifelse(rownames(degs_epithelial_pct0.2$`Epithelial cells_Untr_BaP_VS_Epithelial cells_Untr_Ctrl`) %in% epithelial_genes, "epithelial", "other")))) + 
  geom_point(size = 1) +
  theme_bw() +
  labs(x = "log2 (fold change)", y = "-log10 (q-value)", title = "Untr BaP VS Untr Ctrl") +
  geom_hline(yintercept = -log10(0.05), linetype = 'dotted', col = 'grey25') +
  geom_vline(xintercept = 0, linetype = 'solid', col = 'grey25') +
  xlim(c(-2.42, 2.42)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_color_manual(name = "Cell types", values = c("epithelial" = "#41b349", "other" = "black"))
p7 <- ggplot(data = degs_fibroblast_pct0.2$Fibroblasts_Female_BaP_VS_Fibroblasts_Female_Ctrl, 
             aes(x = avg_log2FC, y = -log10(p_val_adj), 
                 color = factor(ifelse(rownames(degs_fibroblast_pct0.2$Fibroblasts_Female_BaP_VS_Fibroblasts_Female_Ctrl) %in% fibroblast_genes, "fibroblast", "other")))) + 
  geom_point(size = 1) +
  theme_bw() +
  labs(x = "log2 (fold change)", y = "-log10 (q-value)", title = "Female BaP VS Female Ctrl") +
  geom_hline(yintercept = -log10(0.05), linetype = 'dotted', col = 'grey25') +
  geom_vline(xintercept = 0, linetype = 'solid', col = 'grey25') +
  xlim(c(-1.69, 1.69)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_color_manual(name = "Cell types", values = c("fibroblast" = "#fcc307", "other" = "black"))
p8 <- ggplot(data = degs_fibroblast_pct0.2$Fibroblasts_Male_BaP_VS_Fibroblasts_Male_Ctrl, 
             aes(x = avg_log2FC, y = -log10(p_val_adj), 
                 color = factor(ifelse(rownames(degs_fibroblast_pct0.2$Fibroblasts_Male_BaP_VS_Fibroblasts_Male_Ctrl) %in% fibroblast_genes, "fibroblast","other")))) + 
  geom_point(size = 1) +
  theme_bw() +
  labs(x = "log2 (fold change)", y = "-log10 (q-value)", title = "Male BaP VS Male Ctrl") +
  geom_hline(yintercept = -log10(0.05), linetype = 'dotted', col = 'grey25') +
  geom_vline(xintercept = 0, linetype = 'solid', col = 'grey25') +
  xlim(c(-2.39, 2.39)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_color_manual(name = "Cell types", values = c("fibroblast" = "#fcc307", "other" = "black"))
p9 <- ggplot(data = degs_fibroblast_pct0.2$Fibroblasts_Untr_BaP_VS_Fibroblasts_Untr_Ctrl, 
             aes(x = avg_log2FC, y = -log10(p_val_adj), 
                 color = factor(ifelse(rownames(degs_fibroblast_pct0.2$Fibroblasts_Untr_BaP_VS_Fibroblasts_Untr_Ctrl) %in% fibroblast_genes, "fibroblast","other")))) + 
  geom_point(size = 1) +
  theme_bw() +
  labs(x = "log2 (fold change)", y = "-log10 (q-value)", title = "Untr BaP VS Untr Ctrl") +
  geom_hline(yintercept = -log10(0.05), linetype = 'dotted', col = 'grey25') +
  geom_vline(xintercept = 0, linetype = 'solid', col = 'grey25') +
  xlim(c(-1.84, 1.84)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_color_manual(name = "Cell types", values = c("fibroblast" = "#fcc307", "other" = "black"))
patchwork::wrap_plots(p1, p2, p3, p4, p5, p6, p7, p8, p9, ncol = 3)
ggsave('plots/marker_degs_vol.png', width = 36, height = 30, units = 'cm', dpi = 300)


######################## Hormone receptor genes ################################

genes = c('AR', 'ESR1', 'ESR2', 'PGR')
FeaturePlot(sobj, reduction = 'umap', features = genes, ncol = 4) 
ggsave('plots/hormone_receptors_umap.png', width = 18, height = 5, unit = 'in', dpi = 300)
p1 <- AverageExpression(sobj, features = genes, group.by = 'customclassif') %>%
  as.data.frame() %>%
  rownames_to_column(var = 'gene') %>%
  reshape2::melt(id.vars = 'gene',
                 variable.name = 'Annotated cell type',
                 value.name = 'Average expression') %>%
  mutate(`Annotated cell type` = sub('RNA\\.', '', `Annotated cell type`)) %>%
  mutate(`Annotated cell type` = gsub('\\.', ' ', `Annotated cell type`)) %>%
  ggplot(aes(x = gene, y = `Average expression`, fill = `Annotated cell type`)) +
  geom_bar(stat = 'identity', color = 'black', linewidth = .5, position = 'dodge') +
  theme_bw() +
  ggtitle('All cells') +
  theme(plot.title = element_text(hjust = .5)) 
sample_colors <- c("Female Ctrl" = '#93d5dc', "Female BaP" = '#2775b6',  
                   "Untr Ctrl" = '#f7de98', "Untr BaP" = '#fcc307',
                   "Male BaP" = '#41b349', "Male Ctrl" = '#add5a2')
p2 <- AverageExpression(Airway_Progenitor_Cells_sobj, features = genes, group.by = 'sample') %>%
  as.data.frame() %>%
  rownames_to_column(var = 'gene') %>%
  reshape2::melt(id.vars = 'gene',
                 variable.name = 'Sample',  
                 value.name = 'Average expression') %>%
  mutate(Sample = sub('RNA\\.', '', Sample)) %>%
  mutate(Sample = gsub('\\.', ' ', Sample)) %>%
  ggplot(aes(x = gene, y = `Average expression`, fill = Sample)) +  
  geom_bar(stat = 'identity', color = 'black', linewidth = .5, position = 'dodge') +
  theme_bw() +
  ggtitle('Airway Progenitor cells') +
  theme(plot.title = element_text(hjust = .5)) +
  scale_fill_manual(values = sample_colors)
p3 <- AverageExpression(Epithelial_Cells_sobj, features = genes, group.by = 'sample') %>%
  as.data.frame() %>%
  rownames_to_column(var = 'gene') %>%
  reshape2::melt(id.vars = 'gene',
                 variable.name = 'Sample',  
                 value.name = 'Average expression') %>%
  mutate(Sample = sub('RNA\\.', '', Sample)) %>%
  mutate(Sample = gsub('\\.', ' ', Sample)) %>%
  ggplot(aes(x = gene, y = `Average expression`, fill = Sample)) +  
  geom_bar(stat = 'identity', color = 'black', linewidth = .5, position = 'dodge') +
  theme_bw() +
  ggtitle('Epithelial cells') +
  theme(plot.title = element_text(hjust = .5)) +
  scale_fill_manual(values = sample_colors)
p4 <- AverageExpression(Fibroblasts_sobj, features = genes, group.by = 'sample') %>%
  as.data.frame() %>%
  rownames_to_column(var = 'gene') %>%
  reshape2::melt(id.vars = 'gene',
                 variable.name = 'Sample',  # 将变量名更改为'Sample'
                 value.name = 'Average expression') %>%
  mutate(Sample = sub('RNA\\.', '', Sample)) %>%
  mutate(Sample = gsub('\\.', ' ', Sample)) %>%
  ggplot(aes(x = gene, y = `Average expression`, fill = Sample)) +  
  geom_bar(stat = 'identity', color = 'black', linewidth = .5, position = 'dodge') +
  theme_bw() +
  ggtitle('Fibroblasts') +
  theme(plot.title = element_text(hjust = .5)) +
  scale_fill_manual(values = sample_colors)
patchwork::wrap_plots(p1, p2, p3, p4, ncol = 2)
ggsave('plots/average_expression_hormone_receptors.png', width = 24, height = 14, units = 'cm', dpi = 300)


############# Effects of BaP on Gene Expressions of Sex Hormones Receptor Genes
p1 <- ggplot(data = degs_airway_pct0.2$`Airway progenitor cells_Female_Ctrl_VS_Airway progenitor cells_Untr_Ctrl`, 
             aes(x = avg_log2FC, y = -log10(p_val_adj), 
                 color = factor(ifelse(rownames(degs_airway_pct0.2$`Airway progenitor cells_Female_Ctrl_VS_Airway progenitor cells_Untr_Ctrl`) %in% oxphos_genes, "oxphos", 
                                       ifelse(rownames(degs_airway_pct0.2$`Airway progenitor cells_Female_Ctrl_VS_Airway progenitor cells_Untr_Ctrl`) %in% cytran_genes, "cytran", "other"))))) + 
  geom_point(size = 1) +
  theme_bw() +
  labs(x = "log2 (fold change)", y = "-log10 (q-value)", title = "Female Ctrl VS Untr Ctrl") +
  geom_hline(yintercept = -log10(0.05), linetype = 'dotted', col = 'grey25') +
  geom_vline(xintercept = 0, linetype = 'solid', col = 'grey25') +
  xlim(c(-0.7, 0.7)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_color_manual(name = "GO pathway", values = c("aeres" = "#2775b6", "oxphos" = "#41b349", "cytran" = "#fcc307", "other" = "black"))
p2 <- ggplot(data = degs_airway_pct0.2$`Airway progenitor cells_Male_Ctrl_VS_Airway progenitor cells_Untr_Ctrl`,
             aes(x = avg_log2FC, y = -log10(p_val_adj), 
                 color = factor(ifelse(rownames(degs_airway_pct0.2$`Airway progenitor cells_Male_Ctrl_VS_Airway progenitor cells_Untr_Ctrl`) %in% oxphos_genes, "oxphos", 
                                       ifelse(rownames(degs_airway_pct0.2$`Airway progenitor cells_Male_Ctrl_VS_Airway progenitor cells_Untr_Ctrl`) %in% cytran_genes, "cytran", "other"))))) + 
  geom_point(size = 1) +
  theme_bw() +
  labs(x = "log2 (fold change)", y = "-log10 (q-value)", title = "Male Ctrl VS Untr Ctrl") +
  geom_hline(yintercept = -log10(0.05), linetype = 'dotted', col = 'grey25') +
  geom_vline(xintercept = 0, linetype = 'solid', col = 'grey25') +
  xlim(c(-3.04, 3.04)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_color_manual(name = "GO pathway", values = c("oxphos" = "#41b349", "cytran" = "#fcc307", "other" = "black"))
p3 <- ggplot(data = degs_airway_pct0.2$`Airway progenitor cells_Male_Ctrl_VS_Airway progenitor cells_Female_Ctrl`, 
             aes(x = avg_log2FC, y = -log10(p_val_adj), 
                 color = factor(ifelse(rownames(degs_airway_pct0.2$`Airway progenitor cells_Male_Ctrl_VS_Airway progenitor cells_Female_Ctrl`) %in% oxphos_genes, "oxphos", 
                                       ifelse(rownames(degs_airway_pct0.2$`Airway progenitor cells_Male_Ctrl_VS_Airway progenitor cells_Female_Ctrl`) %in% cytran_genes, "cytran", "other"))))) + 
  geom_point(size = 1) +
  theme_bw() +
  labs(x = "log2 (fold change)", y = "-log10 (q-value)", title = "Male Ctrl VS Female Ctrl") +
  geom_hline(yintercept = -log10(0.05), linetype = 'dotted', col = 'grey25') +
  geom_vline(xintercept = 0, linetype = 'solid', col = 'grey25') +
  xlim(c(-3.87, 3.87)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_color_manual(name = "GO pathway", values = c("oxphos" = "#41b349", "cytran" = "#fcc307", "other" = "black"))
p4 <- ggplot(data = degs_epithelial_pct0.2$`Epithelial cells_Female_Ctrl_VS_Epithelial cells_Untr_Ctrl`, 
             aes(x = avg_log2FC, y = -log10(p_val_adj), 
                 color = factor(ifelse(rownames(degs_epithelial_pct0.2$`Epithelial cells_Female_Ctrl_VS_Epithelial cells_Untr_Ctrl`) %in% oxphos_genes, "oxphos", 
                                       ifelse(rownames(degs_epithelial_pct0.2$`Epithelial cells_Female_Ctrl_VS_Epithelial cells_Untr_Ctrl`) %in% cytran_genes, "cytran", "other"))))) + 
  geom_point(size = 1) +
  theme_bw() +
  labs(x = "log2 (fold change)", y = "-log10 (q-value)", title = "Female Ctrl VS Untr Ctrl") +
  geom_hline(yintercept = -log10(0.05), linetype = 'dotted', col = 'grey25') +
  geom_vline(xintercept = 0, linetype = 'solid', col = 'grey25') +
  xlim(c(-0.85, 0.85)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_color_manual(name = "GO pathway", values = c("oxphos" = "#41b349", "cytran" = "#fcc307", "other" = "black"))
p5 <- ggplot(data = degs_epithelial_pct0.2$`Epithelial cells_Male_Ctrl_VS_Epithelial cells_Untr_Ctrl`, 
             aes(x = avg_log2FC, y = -log10(p_val_adj), 
                 color = factor(ifelse(rownames(degs_epithelial_pct0.2$`Epithelial cells_Male_Ctrl_VS_Epithelial cells_Untr_Ctrl`) %in% oxphos_genes, "oxphos", 
                                       ifelse(rownames(degs_epithelial_pct0.2$`Epithelial cells_Male_Ctrl_VS_Epithelial cells_Untr_Ctrl`) %in% cytran_genes, "cytran", "other"))))) + 
  geom_point(size = 1) +
  theme_bw() +
  labs(x = "log2 (fold change)", y = "-log10 (q-value)", title = "Male Ctrl VS Untr Ctrl") +
  geom_hline(yintercept = -log10(0.05), linetype = 'dotted', col = 'grey25') +
  geom_vline(xintercept = 0, linetype = 'solid', col = 'grey25') +
  xlim(c(-1.82, 1.82)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_color_manual(name = "GO pathway", values = c("oxphos" = "#41b349", "cytran" = "#fcc307", "other" = "black"))
p6 <- ggplot(data = degs_epithelial_pct0.2$`Epithelial cells_Male_Ctrl_VS_Epithelial cells_Female_Ctrl`, 
             aes(x = avg_log2FC, y = -log10(p_val_adj), 
                 color = factor(ifelse(rownames(degs_epithelial_pct0.2$`Epithelial cells_Male_Ctrl_VS_Epithelial cells_Female_Ctrl`) %in% oxphos_genes, "oxphos", 
                                       ifelse(rownames(degs_epithelial_pct0.2$`Epithelial cells_Male_Ctrl_VS_Epithelial cells_Female_Ctrl`) %in% cytran_genes, "cytran", "other"))))) + 
  geom_point(size = 1) +
  theme_bw() +
  labs(x = "log2 (fold change)", y = "-log10 (q-value)", title = "Male Ctrl VS Female Ctrl") +
  geom_hline(yintercept = -log10(0.05), linetype = 'dotted', col = 'grey25') +
  geom_vline(xintercept = 0, linetype = 'solid', col = 'grey25') +
  xlim(c(-2.34, 2.34)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_color_manual(name = "GO pathway", values = c("oxphos" = "#41b349", "cytran" = "#fcc307", "other" = "black"))
p7 <- ggplot(data = degs_fibroblast_pct0.2$Fibroblasts_Female_Ctrl_VS_Fibroblasts_Untr_Ctrl, 
             aes(x = avg_log2FC, y = -log10(p_val_adj), 
                 color = factor(ifelse(rownames(degs_fibroblast_pct0.2$Fibroblasts_Female_Ctrl_VS_Fibroblasts_Untr_Ctrl) %in% oxphos_genes, "oxphos", 
                                       ifelse(rownames(degs_fibroblast_pct0.2$Fibroblasts_Female_Ctrl_VS_Fibroblasts_Untr_Ctrl) %in% cytran_genes, "cytran", "other"))))) + 
  geom_point(size = 1) +
  theme_bw() +
  labs(x = "log2 (fold change)", y = "-log10 (q-value)", title = "Female Ctrl VS Untr Ctrl") +
  geom_hline(yintercept = -log10(0.05), linetype = 'dotted', col = 'grey25') +
  geom_vline(xintercept = 0, linetype = 'solid', col = 'grey25') +
  xlim(c(-0.93, 0.93)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_color_manual(name = "GO pathway", values = c("oxphos" = "#41b349", "cytran" = "#fcc307", "other" = "black"))
p8 <- ggplot(data = degs_fibroblast_pct0.2$Fibroblasts_Male_Ctrl_VS_Fibroblasts_Untr_Ctrl, 
             aes(x = avg_log2FC, y = -log10(p_val_adj), 
                 color = factor(ifelse(rownames(degs_fibroblast_pct0.2$Fibroblasts_Male_Ctrl_VS_Fibroblasts_Untr_Ctrl) %in% oxphos_genes, "oxphos", 
                                       ifelse(rownames(degs_fibroblast_pct0.2$Fibroblasts_Male_Ctrl_VS_Fibroblasts_Untr_Ctrl) %in% cytran_genes, "cytran", "other"))))) + 
  geom_point(size = 1) +
  theme_bw() +
  labs(x = "log2 (fold change)", y = "-log10 (q-value)", title = "Male Ctrl VS Untr Ctrl") +
  geom_hline(yintercept = -log10(0.05), linetype = 'dotted', col = 'grey25') +
  geom_vline(xintercept = 0, linetype = 'solid', col = 'grey25') +
  xlim(c(-1.52, 1.52)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_color_manual(name = "GO pathway", values = c("oxphos" = "#41b349", "cytran" = "#fcc307", "other" = "black"))
p9 <- ggplot(data = degs_fibroblast_pct0.2$Fibroblasts_Male_Ctrl_VS_Fibroblasts_Female_Ctrl, 
             aes(x = avg_log2FC, y = -log10(p_val_adj), 
                 color = factor(ifelse(rownames(degs_fibroblast_pct0.2$Fibroblasts_Male_Ctrl_VS_Fibroblasts_Female_Ctrl) %in% oxphos_genes, "oxphos", 
                                       ifelse(rownames(degs_fibroblast_pct0.2$Fibroblasts_Male_Ctrl_VS_Fibroblasts_Female_Ctrl) %in% cytran_genes, "cytran", "other"))))) + 
  geom_point(size = 1) +
  theme_bw() +
  labs(x = "log2 (fold change)", y = "-log10 (q-value)", title = "Male Ctrl VS Female Ctrl") +
  geom_hline(yintercept = -log10(0.05), linetype = 'dotted', col = 'grey25') +
  geom_vline(xintercept = 0, linetype = 'solid', col = 'grey25') +
  xlim(c(-1.61, 1.61)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_color_manual(name = "GO pathway", values = c("oxphos" = "#41b349", "cytran" = "#fcc307", "other" = "black"))
patchwork::wrap_plots(p1, p2, p3, p4, p5, p6, p7, p8, p9, ncol = 3)
ggsave('plots/sex_hornome_degs_vol.png', width = 36, height = 30, units = 'cm', dpi = 300)
```























