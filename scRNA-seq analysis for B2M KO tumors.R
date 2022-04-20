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

# Stacked bar plot function -----------------------------------------------

myStackBarCellComposition <- function(tobj, tcluster_by='ident', tcluster_order=NULL, tgroup_by='orig.ident', tgroup_order=NULL, tcells=NULL, tanncolor=NULL, ttlsize=28, ttxsize=24, tltsize=18){
  # calculate percentage of cells per cluster in each group
  tpercent <- FetchData(tobj, vars=c(tgroup_by, tcluster_by), cells=tcells) %>% dplyr::rename('Cluster'=tcluster_by, 'Group'=tgroup_by) %>%
    rownames_to_column('CellID') %>% group_by(Group, Cluster) %>% summarise_at('CellID', length) %>%
    tidyr::complete(Cluster=factor(tcluster_order), fill=list(CellID=0)) %>%
    group_by(Group) %>% mutate_at('CellID', function(x) { round(x/sum(x)*100,2) })
  colnames(tpercent) <- c("Group", "Cluster", "Percent")
  # reorder group/cluster if provided
  if (! is.null(tcluster_order)){
    tpercent$Cluster <- factor(tpercent$Cluster, levels=tcluster_order)
  }
  if (! is.null(tgroup_order)){
    tpercent$Group <- factor(tpercent$Group, levels=tgroup_order)
  }
  # plot
  tg <- ggplot(tpercent, aes(x=Group, y=Percent, fill=Cluster))
  tg <- tg + geom_bar(stat='identity')
  if (! is.null(tanncolor)){
    tg <- tg + scale_fill_manual(values=tanncolor)
  }
  tg <- tg + ylab('Percentage of cells')
  tg <- tg + theme_bw()
  tg <- tg + theme(axis.title.x=element_blank(), axis.title.y=element_text(color='black', size=ttlsize))
  tg <- tg + theme(axis.text.x=element_text(color='black', size=ttxsize, angle=60, hjust=1), axis.text.y=element_text(color='black', size=ttxsize))
  tg <- tg + theme(legend.text=element_text(size=tltsize), legend.title=element_blank())
  return(tg)
}

# Load files --------------------------------------------------------------
#set working directory to the current folder
setwd("~/Desktop/Projects/scRNA-seq/20220225-MQ833 B2MKO")

B2M_MT.data <- Read10X(data.dir = 'source/B2M_MT/filtered_feature_bc_matrix')
B2M_PT.data <- Read10X(data.dir = 'source/B2M_PT/filtered_feature_bc_matrix')
WT_MT.data <- Read10X(data.dir = '~/Desktop/Projects/scRNA-seq/20211126-MQ833 bilateral tumor and LN/source/1_MI/filtered_feature_bc_matrix')
WT_PT.data <- Read10X(data.dir = '~/Desktop/Projects/scRNA-seq/20211126-MQ833 bilateral tumor and LN/source/3_P/filtered_feature_bc_matrix')

B2M_MT <- CreateSeuratObject(counts = B2M_MT.data, project = "B2M MQ833 tumor", min.cells = 3, min.features = 200)
B2M_PT <- CreateSeuratObject(counts = B2M_PT.data, project = "B2M PBS tumor", min.cells = 3, min.features = 200)
WT_MT <- CreateSeuratObject(counts = WT_MT.data, project = "WT MQ833 tumor", min.cells = 3, min.features = 200)
WT_PT <- CreateSeuratObject(counts = WT_PT.data, project = "WT PBS tumor", min.cells = 3, min.features = 200)

b2m <- merge(B2M_MT, y= c(B2M_PT, WT_MT, WT_PT), add.cell.ids = c("B2M_MT", "B2M_PT", "WT_MT", "WT_PT"), project = "B2M KO MQ833")
b2m$orig.ident <- factor(b2m$orig.ident, levels = c("B2M_MT", "B2M_PT", "WT_MT", "WT_PT"))
#filter mitochondrial genes
b2m[["percent.mt"]] <- PercentageFeatureSet(b2m, pattern = "^mt-")
head(b2m@meta.data, 5)
g <- VlnPlot(b2m, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0)
g <- g + theme(plot.margin = margin(1, 1, 1, 1, "cm"))
g
plotall <- FeatureScatter(b2m, feature1 ="nCount_RNA", feature2 = "nFeature_RNA", group.by = "ident")
plotall

#subsetting
b2m <- subset(b2m, subset = nFeature_RNA > 200 & nFeature_RNA < 7000 & percent.mt < 5)

BiocManager::install("glmGamPoi")
b2m <- SCTransform(b2m, method = "glmGamPoi", vars.to.regress = "percent.mt", verbose = FALSE)

#b2m <- NormalizeData(b2m, normalization.method = "LogNormalize", scale.factor = 10000) #normalize data
#b2m <- FindVariableFeatures(b2m, selection.method = "vst", nfeatures = 5000) #find variable genes
#b2m <- ScaleData(b2m) #scale the data, only for variable features
b2m <- RunPCA(b2m, features = VariableFeatures(object = b2m)) #PCA

ElbowPlot(b2m) #determine the dimentionality of the dataset by identifying the 'elbow' point in this graph

b2m <- FindNeighbors(b2m, dims = 1:15)
b2m <- FindClusters(b2m, resolution = 0.5)
b2m <- RunUMAP(b2m, dims = 1:15) #use the same dims as the FindNeighbors function
b2m@meta.data

saveRDS(b2m, file = "seurat objects/b2m.rds")

#plot UMAP seperated by sample
DimPlot(b2m, reduction = "umap", split.by = "orig.ident", cols = palette, ncol = 2)

g <- g + theme(plot.margin = margin(1, 1, 1, 1, "cm"))
tiff(filename = "figures/cluster_by_sample_LN.tiff", width=14, height=4, units="in", res = 600)
g
dev.off()
levels(b2m)
levels(x = b2m.renamed) <- c("Naïve CD4+ T cells", "Naïve CD8+ T cells", "Treg", "Ccl5+ NKT cells", 
                               "Proliferating NKT cells", "Naïve B cells", "ISG+ B cells", "Non-classical Mono", "Unknown")

#heatmap
b2m.markers <- FindAllMarkers(b2m, slot = "scale.data", only.pos = TRUE, 
                                    min.pct = 0.25, logfc.threshold = 0.25, max.cells.per.ident = 200)
b2m.markers %>%
  group_by(cluster) %>%
  slice_max(prop = 0.1, order_by = avg_diff) -> filtered_markers
b2m.downsampled <- subset(b2m, downsample = 20)
mat <- GetAssayData(b2m.downsampled, slot = "scale.data", assay = "SCT")
cluster_info <- sort(b2m.downsampled@active.ident)
cluster_info
mat <- as.matrix(mat[filtered_markers$gene, names(cluster_info)])

col_fun = colorRamp2(c(-2, 0, 2), c("#bbdded", "white", "#9c322a"))

g<- Heatmap(mat,
            cluster_rows = FALSE,
            cluster_columns = FALSE,
            show_column_names = FALSE,
            show_row_names = TRUE, 
            row_names_gp = gpar(fontsize = 5),
            col=col_fun,
            column_split = cluster_info,
            column_title = c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14),
            column_title_rot = 0,
            column_title_side = "top",
            row_split = filtered_markers$cluster,
            row_gap = unit(0.5, "mm"), 
            column_gap = unit(0.5, "mm"), 
            border = FALSE,
            border_gp = gpar(col = "black", lty = "dotted", lwd = 1.2) ,
            row_title = NULL)
g <- g + theme(plot.margin = margin(1, 1, 1, 1, "cm"))
g
?gpar
column_ha = HeatmapAnnotation(foo = anno_block(gp = gpar(fill = palette, col = 0), height = unit(2, "mm")), border = FALSE)
row_ha = rowAnnotation(foo = anno_mark(at = c(3, 34, 260, 228, 45, 71, 73, 103, 47, 138, 141, 366, 369, 220, 202, 153,268, 269, 261), 
                                       labels = c("Cd8a", "Cd4", "Cd79a", "Nkg7", "Isg15", 'Tnfrsf4', 'Foxp3', 
                                                  'Ccl5', 'Ifit3', 'Mki67', 'Top2a', 'Fcer1g', 'Cst3', "Ncr1", "Eif5a", "Gzmb", "Ighm", "Ighd", "Cd79b")))
b2m.sct.renamed@active.ident <- plyr::mapvalues(x = b2m.sct.renamed@active.ident, 
                                                  from = c("11: Pre-B cells", "12: Non-classical monocytes"), 
                                                  to = c("10: Pre-B cells", "11: Non-classical monocytes"))
which(filtered_markers$gene == 'Cd79b')


tiff(filename = "figures/heatmap.tiff", width=14, height=18, units="in", res = 600)
DoHeatmap(b2m.small, features = top10_new$gene) + scale_fill_gradient2(low="royalblue3", high="firebrick3")
dev.off()

#feature plot
FeaturePlot(b2m, features = c("Ifng", 'Gzmb'))


VlnPlot(b2m, features = c('Cd8a', 'Cd4','Cd7', 'Tcf7', 'Cd40lg', 'Isg15', 'Ifit3', 'Gzmb', 'Ifng', 'Klrg1', 'Tox','Havcr2', 'Foxp3', 'Tnfrsf4', 'Mki67'), 
        stack = TRUE, pt.size = 0, cols = palette, flip = TRUE, fill.by = 'ident') + NoLegend()+ xlab('Clusters')


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
es.max = sctype_score(scRNAseqData = b2m[["SCT"]]@scale.data, scaled = TRUE, 
                      gs = gs_list$gs_positive, gs2 = gs_list$gs_negative)
end_time <- Sys.time()
end_time - start_time
# merge by cluster
cL_resutls = do.call("rbind", lapply(unique(b2m@meta.data$seurat_clusters), function(cl){
  es.max.cl = sort(rowSums(es.max[ ,rownames(b2m@meta.data[b2m@meta.data$seurat_clusters==cl, ])]), decreasing = !0)
  head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(b2m@meta.data$seurat_clusters==cl)), 10)
}))
sctype_scores = cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores)  

# set low-confident (low ScType score) clusters to "unknown"
sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] = "Unknown"
sctype_scores <- sctype_scores[order(as.numeric(as.character(sctype_scores$cluster))),]
print(sctype_scores[,1:3])

b2m@meta.data$customclassif = ""
for(j in unique(sctype_scores$cluster)){
  cl_type = sctype_scores[sctype_scores$cluster==j,]; 
  b2m@meta.data$customclassif[b2m@meta.data$seurat_clusters == j] = as.character(cl_type$type[1])
}


# Rename clusters ---------------------------------------------------------

b2m.renamed <- b2m
current.cluster.ids <- c(0,1,2,3,4,5,6,7,8,9,10,11, 12, 13, 14)
new.cluster.ids <- c("7. M1 Mφ", "11. cDC2", "6. Inflammatory Mono", "1. Tcf7+ T cells", "8. Nos2+ M1 Mφ", 
                     "5. NK cells", "2. Proliferating CD8+ T cells", "3. Effector CD8+ T cells", "10. Neutrophils", 
                     "9. M2 Mφ", "12. Migrating DCs", "13. pDCs", "14. cDC1", "4. CD8+ NKT-like cells", 
                     "15. B cells")
b2m.renamed@active.ident <- plyr::mapvalues(x = b2m.renamed@active.ident, 
                                              from = current.cluster.ids, 
                                              to = new.cluster.ids)
levels(x = b2m.renamed) <- c("1. Tcf7+ T cells",  "2. Proliferating CD8+ T cells", "3. Effector CD8+ T cells",
                             "4. CD8+ NKT-like cells", "5. NK cells", "6. Inflammatory Mono",
                             "7. M1 Mφ","8. Nos2+ M1 Mφ","9. M2 Mφ",  "10. Neutrophils",  "11. cDC2", 
                             "12. Migrating DCs", "13. pDCs", "14. cDC1", "15. B cells")

saveRDS(b2m.renamed, file = "seurat objects/b2m_renamed.rds")

# Plotting after cluster annotation ---------------------------------------

DimPlot(b2m.renamed, reduction = "umap", cols = new_palette, ncol = 2, split.by = "orig.ident")

myStackBarCellComposition(b2m.renamed, tgroup_order=c("B2M MQ833 tumor", "B2M PBS tumor", "WT MQ833 tumor", "WT PBS tumor"), 
                          tcells=NULL, tanncolor=new_palette, ttlsize=20, ttxsize=18, tltsize=16)

#heatmap
b2m.renamed.markers <- FindAllMarkers(b2m.renamed, assay = "SCT", slot = "scale.data", only.pos = TRUE, 
                              min.pct = 0.25, logfc.threshold = 0.25, max.cells.per.ident = 200)
b2m.renamed.markers %>%
  group_by(cluster) %>%
  slice_max(n = 10, order_by = avg_diff) -> filtered_markers
b2m.downsampled <- subset(b2m.renamed, downsample = 20)
mat <- GetAssayData(b2m.downsampled, slot = "scale.data", assay = "SCT")
cluster_info <- sort(b2m.downsampled@active.ident)
cluster_info
mat <- as.matrix(mat[filtered_markers$gene, names(cluster_info)])

col_fun = colorRamp2(c(-1, 0, 3), c("#bbdded", "white", "#9c322a"))

ht <- Heatmap(mat,
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        show_column_names = FALSE,
        show_row_names = FALSE, 
        row_names_gp = gpar(fontsize = 5),
        col=col_fun,
        column_split = cluster_info,
        column_title = c("1. Tcf7+ T cells",  "2. Proliferating CD8+ T cells", "3. Effector CD8+ T cells",
                         "4. CD8+ NKT-like cells", "5. NK cells", "6. Inflammatory Mono",
                         "7. M1 Mφ","8. Nos2+ M1 Mφ","9. M2 Mφ",  "10. Neutrophils",  "11. cDC2", 
                         "12. Migrating DCs", "13. pDCs", "14. cDC1", "15. B cells"),
        column_title_rot = 90, 
        column_title_side = "bottom",
        row_split = filtered_markers$cluster,
        row_gap = unit(0.5, "mm"), 
        column_gap = unit(0.5, "mm"), 
        border = FALSE,
        border_gp = gpar(col = "black", lty = "dotted", lwd = 1.2) ,
        row_title = NULL, 
        top_annotation = column_ha, 
        right_annotation = row_ha)
draw(ht, padding = unit(c(25, 2, 2, 2), "mm"))
pdf(file = "figures/heatmap_simplified.pdf", width = 6, height = 8)
warnings()
dev.off()
column_ha = HeatmapAnnotation(foo = anno_block(gp = gpar(fill = new_palette, col = 0), height = unit(2, "mm")), border = FALSE)
row_ha = rowAnnotation(foo = anno_mark(at = c(178, 186, 198, 182, 
                                              6, 12, 190, 
                                              182, 238, 106,  
                                              85, 109, 
                                              100, 55, 
                                              92, 272, 273, 
                                              252, 352, 375, 
                                              411, 419, 53,
                                              431, 434, 233, 340), 
                                       labels = c("Chil3", 'Cxcl10','Ccl12', "Fcgr1", 
                                                  'Tcf7', 'Cd3d', 'Ccl7', 
                                                  'Fcgr1', 'Il1rn', 'Gzmb',  
                                                  'Ifng', 'Nkg7', 
                                                  'Gzma', 'Mki67', 
                                                  'Ctla4', 'S100a9', 'S100a8', 
                                                  'Apoe', 'Ccr7', 'Siglech', 
                                                  'Cst3', 'Wdfy4', 'Top2a',
                                                  'Cd79a', 'Cd79b', 'Nos2', 'Cd209a') ))

which(filtered_markers$gene == 'Ifit3')

new_palette <-c("#A6CEE3", "#1F78B4", "#66C2A5", "#B2DF8A", "#33A02C","#CAB2D6", "#6A3D9A", "#f29df5", "#cc55d9", "#B15928", "#FB9A99","#E31A1C", "#FDBF6F", "#FF7F00", "#FFFF99")

#pathway analysis
immune_gene_sets = msigdbr(species = "mouse", category = "C7", subcategory = "IMMUNESIGDB")
head(immune_gene_sets)
#retrieve all mouse gene sets
m_df <- msigdbr(species = "Mus musculus")
head(m_df, 2) %>% as.data.frame
rm(m_df)
#retrieve the hallmark gene set
m_t2g <- msigdbr(species = "Mus musculus", category = "H") %>% 
  dplyr::select(gs_name, gene_symbol)

nos2m1.markers <- FindMarkers(object = b2m.renamed, ident.1 = "B2M MQ833 tumor", ident.2 = "B2M PBS tumor", 
                                   group.by = "orig.ident", subset.ident = "1. Tcf7+ T cells", min.cells.group = 20)


genes <- rownames(nos2m1.mq833.markers)
genelist <- nos2m1.mq833.markers$avg_log2FC
names(genelist) <- genes
head(genelist)
genelist <- na.omit(genelist)
genelist <- sort(genelist, decreasing = TRUE)
em <- clusterProfiler::GSEA(genelist, TERM2GENE = m_t2g, eps = 0)

#dotplot: generatio=count/set size
tiff(filename = file.path(figdir, "GSEA_c3_WT MQ833 vs WT PBS.tiff"), width=6, height=4, units="in", res = 600)
dotplot(em, showCategory = 10, title = "B2M KO tumor vs. B2M KO PBS in Nos2+ M1 Mφ" , split=".sign") + facet_grid(.~.sign)
dev.off()

gseaplot2(em, geneSetID = 3, title = em$Description[3])



#macrophage DE analysis 
tcf7t.mq833.markers <- FindMarkers(object = b2m.renamed, ident.1 = "B2M MQ833 tumor", ident.2 = "WT MQ833 tumor", 
                          group.by = "orig.ident", subset.ident = "1. Tcf7+ T cells", min.cells.group = 20)

write.csv(tcf7t.mq833.markers,"tables/tcf7t_b2m mq833 vs wt mq833_markers.csv")

keyvals <- ifelse(
  tcf7t.mq833.markers$p_val_adj > 10e-10, 'grey70',
  ifelse(tcf7t.mq833.markers$avg_log2FC < -0.5, 'royalblue',
         ifelse(tcf7t.mq833.markers$avg_log2FC > 0.5, 'deeppink2',
                'grey70')))
keyvals[is.na(keyvals)] <- 'grey70'
names(keyvals)[keyvals == 'deeppink2'] <- 'Up in B2M KO'
names(keyvals)[keyvals == 'grey70'] <- 'non-significant'
names(keyvals)[keyvals == 'royalblue'] <- 'Down in B2M KO'

g <- EnhancedVolcano(tcf7t.mq833.markers,
                     lab = rownames(tcf7t.mq833.markers),
                     x = 'avg_log2FC',
                     y = 'p_val_adj',
                     title = 'Tcf7+ T cells',
                     subtitle = bquote(italic('B2M KO MQ833 vs. WT MQ833')),
                     pCutoff = 10e-10,
                     FCcutoff = 0.5,
                    # selectLab = c('Iigp1', 'Cxcl2', 'Ccl4', 'Gbp2', 'Cxcl10', 'Saa3', 
                    #               'Ccl12', 'Cxcl9', 'Isg15', 'Fcgr4', 'Chil3', 'Cd274', 
                    #               'Apoe', 'Fn1', 'Ccl6'),
                     cutoffLineType = 'blank',
                     colCustom = keyvals,
                     pointSize = 2.0,
                     labSize = 3.0,
                     # labFace = 'bold',
                     # boxedLabels = TRUE,
                     drawConnectors = TRUE,
                     widthConnectors = 0.3,
                     arrowheads = FALSE,
                     maxoverlapsConnectors = 40,
                     gridlines.major = FALSE,
                     gridlines.minor = FALSE,
                     legendDropLevels = TRUE,
                     colAlpha = 1,
                     legendPosition = 'right',
                     legendIconSize = 2.0)
g <- g + theme(plot.margin = margin(1, 1, 1, 1, "cm"))
g

# DE genes for nos2 m1 cluster
# 'Saa3', 'Nos2', 'Ifitm1', 
# 'Gzmb', 'Cxcl10', 'Gzma', 'Cd274', 
# 'Nkg7', 'Ifng', 'Cd3d', 'Il27', 'Apoe', 'H2-Ab1', 'H2-Aa', 'H2-Eb1', 'Fn1'

nos2m1 <- WhichCells(b2m.renamed, idents = c("8. Nos2+ M1 Mφ"))
DimPlot(b2m.renamed, reduction = 'umap', cells.highlight = nos2m1,cols.highlight = "#f29df5", cols = "lightgrey")

FeaturePlot(b2m, features = c('Itgam', 'Cd3d'), blend = TRUE, order = TRUE, cols = c('red', 'blue'), blend.threshold = 0.5, pt.size = 0.3)

#average expression and heatmap
avgexp <- AverageExpression(b2m.renamed, assay = "SCT", slot = "scale.data", features = filtered_markers$gene, group.by = "ident", return.seurat = TRUE)

DoHeatmap(avgexp, features = filtered_markers$gene) + scale_fill_gradient2(low="royalblue3", high="firebrick3")
DimPlot(b2m.renamed, reduction = 'umap', cols = new_palette)


#expression of gene signature
library(stringr)
ifng.gene <- "GBP6	SP110	MX1	RSAD2	MX2	PSMB9	PML	NUP93	LAP3	MTHFD2	IFI27	TRIM14	RBCK1	IDO1	IFITM2	IFITM3	FGL2	NOD1	UBE2L6	PSMB10	PNP	APOL6	HERC6	ST3GAL5	SLAMF7	B2M	GBP4	EPSTI1	BATF2	ST8SIA4	PLA2G4A	HLA-A	HLA-B	HLA-G	IFI44	PARP14	LATS2	PARP12	VAMP8	PLSCR1	BANK1	SERPING1	PSME2	SSPN	PSME1	VAMP5	CFB	HLA-DRB1	CFH	IFI35	HIF1A	IFIT2	IFI30	IFIT1	ICAM1	IFIT3	RNF213	MT2A	HELZ2	KLRK1	TNFSF10	CMKLR1	NCOA3	GPR18	PNPT1	IL15	IL10RA	NMI	BST2	RAPGEF6	IL18BP	CD40	BTG1	SECTM1	FPR1	SOCS3	SOCS1	RIPK2	CD38	RIPK1	HLA-DQA1	P2RY14	STAT1	STAT2	GCH1	STAT3	STAT4	PSMA2	PSMA3	CMPK2	PELI1	MARCHF1	CDKN1A	OGFR	PTGS2	IFIH1	PSMB8	NAMPT	PSMB2	METTL7B	HLA-DMA	JAK2	PTPN1	ZBP1	PTPN2	IL4R	CD74	TOR1B	TAP1	VCAM1	CXCL10	CXCL11	IFI44L	ZNFX1	TRAFD1	FAS	LY6E	CD69	PTPN6	CD86	SAMD9L	RNF31	CMTR1	DDX58	AUTS2	ARID5B	TDRD7	ISG15	NFKB1	ISG20	MYD88	SELP	XAF1	IFNAR2	CD274	LGALS3BP	SPPL2A	CSF2RB	DDX60	PIM1	XCL1	EIF2AK2	TAPBP	IRF1	NFKBIA	IL7	IL6	OAS2	IRF4	IRF5	OAS3	IRF2	IRF8	IRF9	IRF7	LCP2	EIF4E3	ARL4A	CXCL9	TNFAIP6	TNFAIP3	NLRC5	TNFAIP2	ADAR	CASP8	CASP7	CASP4	PDE4B	CASP3	ISOC1	CASP1	CIITA	IL15RA	GZMA	LYSMD2	BPGM	WARS1	SOD2	IL2RB	TXNIP	PFKP	C1R	RTP4	C1S	MVP	SRI	USP18	OASL	SAMHD1	SLC25A28	FCGR1A	CCL7	DHX58	CCL5	TRIM25	TRIM26	CCL2	ITGB7	TRIM21	UPP1"
ifng.gene <- strsplit(ifng.gene, "\\s+")[[1]]
ifng.gene
ifng.gene <- str_to_title(ifng.gene)

ifng.sig <- list(ifng.gene)

m1.sig <- c("ABCA9","ABHD1","ACSS1","ACYP2","ADA","AIF1","ALAD","ARFGEF3","ARSB","ATP6V0E2","AXL","BCKDHB","BLOC1S6","CADM1","CAP1","CAPN5","CBX6","CD59","CFH","CLBA1","CNRIP1","COLEC12","COMT","CRIM1","CXCL14","CXCR4","DST","DYNLT1","EMC1","ENO2","FAM124A","FAM135A","FAM9B","FGD2","FILIP1L","GALNT11","GATM","GDA","GJA1","GLO1","GNB4","HAUS2","HDDC3","HLA-DQA1","HMGN3","KCNJ10","LAMA3","LCORL","LYPLAL1","MAF","MALAT1","MARCKSL1","MARCO","MSR1","NAT8L","NRCAM","OCEL1","OGFRL1","P2RY13","PIANP","PIK3AP1","PLAAT3","PLBD1","PLXDC2","PPP2R5C","PTGER3","RAB10","RAPSN","RASAL2","RCBTB2","RCN1","RFX3","RPL14","SFI1","SLC35A1","SLC7A7","SLCO2B1","SRD5A3","TGFBI","TIFAB","TM7SF3","TOR3A","TTC3","TUBB2B","TXNIP","ZNF681") 
m1.sig <- str_to_title(m1.sig)
m1.sig <- list(m1.sig)

ifna.sig <- "IL4R	CD74	SP110	MX1	RSAD2	TAP1	PSMB9	CXCL10	CXCL11	IFI44L	LAP3	LPAR6	IFI27	TRAFD1	TRIM14	LY6E	IFITM2	IFITM3	SAMD9	CSF1	IFITM1	UBE2L6	SAMD9L	RNF31	HERC6	LAMP3	CCRL2	MVB12A	B2M	CMTR1	GBP2	GBP4	EPSTI1	BATF2	HLA-C	TENT5A	TDRD7	ISG15	IFI44	PARP14	MOV10	ISG20	PARP12	PLSCR1	PROCR	SELL	PSME2	PSME1	CNP	LGALS3BP	IFI35	IFIT2	IFI30	DDX60	IFIT3	NUB1	HELZ2	PNPT1	IL15	NCOA7	EIF2AK2	NMI	OAS1	BST2	IRF1	ELF1	IL7	IRF2	IRF9	IRF7	UBA7	GMPR	ADAR	CASP8	TRIM5	RIPK2	CASP1	STAT2	WARS1	PARP9	PSMA3	TXNIP	CMPK2	CD47	RTP4	C1S	OGFR	TMEM140	USP18	OASL	SLC25A28	IFIH1	PSMB8	DHX58	TRIM25	TRIM26	TRIM21"
ifna.sig <- strsplit(ifna.sig, "\\s+")[[1]]  
ifna.sig <- str_to_title(ifna.sig)
ifna.sig <- list(ifna.sig)


b2m.renamed <- AddModuleScore(object = b2m.renamed, features = ifna.sig, name = "ifna_sig_score")
FeaturePlot(b2m.renamed, features = "ifna_sig_score1", split.by = "orig.ident")







