# install packages --------------------------------------------------------
BiocManager::install("openxlsx")
library(dplyr)
library(Seurat)
library(patchwork)
library(RColorBrewer)
library(tibble, quietly=T, warn.conflicts=F)
library(magrittr)
library(ggplot2)
library(gplots)
library(EnhancedVolcano)
library(ComplexHeatmap)
library(devtools)
library(presto)
library(dittoSeq)
library(SeuratWrappers)
library(monocle3)
library(TSCAN)
library(DOSE)

library(enrichplot)
library(clusterProfiler)
library(msigdbr)
library(Cairo)
library(circlize)


#presets & useful functions
palette <- c(brewer.pal(12, "Paired"), brewer.pal(8, "Set2"), brewer.pal(3, "Set3")) #color palette
readRDS(file = "file path") #read Seurat object
saveRDS(mq833.renamed, file = "file path") #save Seurat object
ggsave(file.path(figdir, 'DE_c9_WT_MQ833_vs_STAT2KO_MQ833.tiff'), plot=g, width=6.5, height=8) #save ggplot figure

# Load files --------------------------------------------------------------
#set working directory to the current folder
setwd("~/Desktop/Projects/scRNA-seq/20211126-MQ833 bilateral tumor and LN")

MQ833_I.data <- Read10X(data.dir = 'source/1_MI/filtered_feature_bc_matrix')
MQ833_NI.data <- Read10X(data.dir = 'source/2_MNI/filtered_feature_bc_matrix')
PBS.data <- Read10X(data.dir = 'source/3_P/filtered_feature_bc_matrix')
MQ833_TDLN.data <- Read10X(data.dir = 'source/4_MLN/filtered_feature_bc_matrix')
PBS_TDLN.data <- Read10X(data.dir = 'source/5_PLN/filtered_feature_bc_matrix')


MQ833_I <- CreateSeuratObject(counts = MQ833_I.data, project = "MQ833 injected tumor", min.cells = 3, min.features = 200)
MQ833_NI <- CreateSeuratObject(counts = MQ833_NI.data, project = "MQ833 non-injected tumor", min.cells = 3, min.features = 200)
MQ833_TDLN <- CreateSeuratObject(counts = MQ833_TDLN.data, project = "MQ833 TDLN", min.cells = 3, min.features = 200)
PBS <- CreateSeuratObject(counts = PBS.data, project = "PBS injected tumor", min.cells = 3, min.features = 200)
PBS_TDLN <- CreateSeuratObject(counts = PBS_TDLN.data, project = "PBS TDLN", min.cells = 3, min.features = 200)

#merge all 5 samples
mbi <-merge(MQ833_I , y = c(MQ833_NI, PBS, MQ833_TDLN, PBS_TDLN), 
            add.cell.ids = c("MQ833_I", "MQ833_NI", "PBS", "MQ833_TDLN", "PBS_TDLN"), project = "MQ833 bilateral")
mbi
head(mbi)

mbi$orig.ident <- factor(mbi$orig.ident, levels = c("MQ833 injected tumor", "MQ833 non-injected tumor", "PBS injected tumor", "MQ833 TDLN", "PBS TDLN"))


#filter mitochondrial genes
mbi[["percent.mt"]] <- PercentageFeatureSet(mbi, pattern = "^mt-")
head(mbi@meta.data, 5)
g <- VlnPlot(mbi, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0)
g <- g + theme(plot.margin = margin(1, 1, 1, 1, "cm"))
g
plotall <- FeatureScatter(mbi, feature1 ="nCount_RNA", feature2 = "nFeature_RNA", group.by = "ident")
plotall

#subsetting
mbi <- subset(mbi, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 5)

mbi <- NormalizeData(mbi, normalization.method = "LogNormalize", scale.factor = 10000) #normalize data
mbi <- FindVariableFeatures(mbi, selection.method = "vst", nfeatures = 5000) #find variable genes
mbi <- ScaleData(mbi) #scale the data, only for variable features
mbi <- RunPCA(mbi, features = VariableFeatures(object = mbi)) #PCA

ElbowPlot(mbi) #determine the dimentionality of the dataset by identifying the 'elbow' point in this graph

mbi <- FindNeighbors(mbi, dims = 1:15)
mbi <- FindClusters(mbi, resolution = 0.5)
mbi <- RunUMAP(mbi, dims = 1:15) #use the same dims as the FindNeighbors function
mbi

#plot UMAP seperated by sample
g <- DimPlot(mbi, reduction = "umap",split.by = 'orig.ident', cols = palette)
g <- g + theme(plot.margin = margin(1, 1, 1, 1, "cm"))
tiff(filename = "figures/cluster_by_sample.tiff", width=14, height=4, units="in", res = 600)
g
dev.off()

#heatmap
mbi.markers <- FindAllMarkers(mbi, only.pos = TRUE, min.pct = 0.5, logfc.threshold = 0.25, max.cells.per.ident = 200)
mbi.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
mbi.small <- mbi[, sample(colnames(mbi), size =2999, replace=F)]

tiff(filename = "figures/heatmap.tiff", width=14, height=18, units="in", res = 600)
DoHeatmap(mbi.small, features = top10$gene) + scale_fill_gradient2(low="royalblue3", high="firebrick3")
dev.off()


cluster_info <- sort(mbi$seurat_clusters)

mat <- GetAssayData(mbi, slot = "counts")
mat <- ScaleData(mat)

mat <- as.matrix(mat[top10$gene, names(cluster_info)])

palette.named <- c("0" = "#A6CEE3", "1" = "#1F78B4", "2" = "#B2DF8A", "3" = "#33A02C", 
                   "4" = "#FB9A99", "5" = "#E31A1C", "6" = "#FDBF6F", "7" = "#FF7F00", 
                   "8" = "#CAB2D6", "9" = "#6A3D9A", "10" = "#FFFF99", "11" = "#B15928", 
                   "12" = "#66C2A5", "13" = "#FC8D62", "14" = "#8DA0CB", "15" = "#E78AC3",
                   "16" = "#A6D854")
ha = HeatmapAnnotation(seurat_clusters = 1:33896, col = list(seurat_clusters = palette.named))
palette
col_fun = colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))
Heatmap(mat, cluster_rows = FALSE,
        cluster_columns = FALSE,
        show_column_names = TRUE,
        show_row_names = TRUE,
        use_raster = TRUE, 
        col = col_fun, 
        heatmap_width = unit(12, "cm"), 
        heatmap_height = unit(18, "cm"), 
        top_annotation = ha)
ncol(mat)

?dittoHeatmap()

FeaturePlot(mbi, features = c('Cd3d', 'Cd3g', 'Tcf7')) #lymphocytes
FeaturePlot(mbi, features = 'Ifng', split.by = 'orig.ident')
FeaturePlot(mbi, features = c('Cxcl2', 'S100a8', 'S100a9')) #neutrophils


# Sctype ------------------------------------------------------------------


# load libraries and functions

lapply(c("dplyr","Seurat","HGNChelper"), library, character.only = T)
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R"); source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")

# load gene set preparation function
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
# load cell type annotation function
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")
# DB file
db_ = "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx";
tissue = "Immune system" # e.g. Immune system, Liver, Pancreas, Kidney, Eye, Brain

# prepare gene sets
gs_list = gene_sets_prepare(db_, tissue)

# get cell-type by cell matrix
start_time <- Sys.time()
es.max = sctype_score(scRNAseqData = mbit[["RNA"]]@scale.data, scaled = TRUE, 
                      gs = gs_list$gs_positive, gs2 = gs_list$gs_negative)
end_time <- Sys.time()
end_time - start_time
# merge by cluster
cL_resutls = do.call("rbind", lapply(unique(mbit@meta.data$seurat_clusters), function(cl){
  es.max.cl = sort(rowSums(es.max[ ,rownames(mbit@meta.data[mbit@meta.data$seurat_clusters==cl, ])]), decreasing = !0)
  head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(mbit@meta.data$seurat_clusters==cl)), 10)
}))
sctype_scores = cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores)  

# set low-confident (low ScType score) clusters to "unknown"
sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] = "Unknown"
print(sctype_scores[,1:3])

mbit@meta.data$customclassif = ""
for(j in unique(sctype_scores$cluster)){
  cl_type = sctype_scores[sctype_scores$cluster==j,]; 
  mbit@meta.data$customclassif[mbit@meta.data$seurat_clusters == j] = as.character(cl_type$type[1])
}

DimPlot(mbit, reduction = "umap", cols = palette)        



start_time <- Sys.time()

end_time <- Sys.time()

end_time - start_time



