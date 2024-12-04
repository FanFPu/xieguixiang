### Seurat clustering
library(Seurat)
library(dplyr)
library(patchwork)
library(ggplot2)
library(future)
library(Cairo)
library(harmony)
options(bitmapType = "cairo")

maxGenes <- 8000
SampleID <- "BLCA"
minUMIs <- 700
minGenes <- 200
maxPercent.mt <-10
dim.usage <- 30
res.usage <- 1.0
doublets.percentage <- 0.075
workdir <- "/jdfsbjcas1/ST_BJ/P21H28400N0232/xieguixiang/BLCA/SC/R28BC/myocyte/intersectHVG_10/"
organs <- "BLCA"

source("/jdfsbjcas1/ST_BJ/P21H28400N0232/wangdi/tumors_public/codes/scRNA_primary.R")
source("/jdfsbjcas1/ST_BJ/P21H28400N0232/wangdi/tumors_public/codes/colors.R")
source("/jdfsbjcas1/ST_BJ/P21H28400N0232/wangdi/tumors_public/codes/signatureGenes.R")
source("/jdfsbjcas1/ST_BJ/P21H28400N0232/wangdi/tumors_public/codes/0_initNEWanalysis.R")
source("/jdfsbjcas1/ST_BJ/P21H28400N0232/wangdi/tumors_public/codes/DEG_wilcox.R")

setwd(workdir)
createObjdir(workdir = workdir)

x <- load('/jdfsbjcas1/ST_BJ/P21H28400N0232/xieguixiang/BLCA/SC/R28BC/myocyte/intersectHVG_10/saveData/Muscle_hvg_10.RData')
seurat_comb_singlet <- get(x)

# 找出所有不以"KRT"开头的基因
keep_genes <- !grepl("^KRT", rownames(seurat_comb_singlet))
# 使用这些基因创建一个新的Seurat子集
seurat_comb_singlet <- subset(seurat_comb_singlet, features = keep_genes)

# # seurat_comb_singlet <- subset(seurat_comb_singlet, cellTypes_new != 'Undefined')
# seurat_comb_singlet <- subset(seurat_comb_singlet, cellTypes_new == 'Fibroblast')
# str(seurat_comb_singlet)
# seurat_comb_singlet <- readRDS('/jdfsbjcas1/ST_BJ/P21H28400N0232/xieguixiang/BLCA/TSanalysis/SCdata/tumor.rds')
# Idents(seurat_comb_singlet) <- seurat_comb_singlet@meta.data$TS

#### normalize and reduce the data ####
# seurat_comb_singlet <- subset(seurat_comb_singlet, Patients %in% c("HBCP1", "HBCP3A", "HBCP3B", "HBCP3C","HBCP5A","HBCP5B","HBCP6","HBCP7","HBCP8","HBCP9","HBCP10","HBCP11","HBCP12","HBCP13","HBCP14","HBCP15","HBCP17","HBCP18"))
# print(table(seurat_comb_singlet$Patients))
# seurat_comb_singlet <- NormalizeData(seurat_comb_singlet)
# seurat_comb_singlet <- FindVariableFeatures(seurat_comb_singlet, selection.method = "vst", nfeatures = 3000)
# all.genes <- rownames(seurat_comb_singlet)
# seurat_comb_singlet <- ScaleData(seurat_comb_singlet, features = all.genes)

#### harmony
# seurat_comb_singlet <- PercentageFeatureSet(seurat_comb_singlet, pattern = "^MT-", col.name = "percent.mt")# run sctransform
# seurat_comb_singlet <- SCTransform(seurat_comb_singlet, vars.to.regress = "percent.mt", verbose = FALSE)

#### use common genes among samples
geneset_test <- c()
for (i in c("HBCP1","HBCP2","HBCP3A","HBCP3B","HBCP3C","HBCP4A","HBCP4B","HBCP5A","HBCP5B","HBCP6","HBCP7","HBCP8","HBCP9","HBCP10","HBCP11","HBCP12","HBCP13","HBCP14","HBCP15","HBCP16","HBCP17","HBCP18","HBCP19","HBCP20","HBCP21","HBCP22")) {
    sc_test <- subset(seurat_comb_singlet, Patients == i)
    sc_test <- NormalizeData(sc_test)
    sc_test <- FindVariableFeatures(sc_test, selection.method = "vst", nfeatures = 2000)
    geneset_test <- c(geneset_test, VariableFeatures(sc_test))
}
geneset <- c()
for (i in geneset_test) {
    if (table(geneset_test)[i] >= 15) {
        geneset <- c(geneset, i)
    }
}
geneset <- unique(geneset)
print(length(geneset))
# print(geneset)

seurat_comb_singlet <- RunPCA(seurat_comb_singlet, features = geneset)
# cordata <- data.frame(seurat_comb_singlet@reductions[["pca"]]@cell.embeddings)
# cordata <- cbind(cordata, seurat_comb_singlet$Patients)
# colnames(cordata)[colnames(cordata) == 'seurat_comb_singlet$Patients'] <- 'Patients'
# print(colnames(cordata))
# cordata$Patients <- as.numeric(as.factor(cordata$Patients))
# for (i in 1:30) {
#     print(colnames(cordata[i]))
#     print(cor(cordata[i], cordata[51], method="pearson"))
# }

# seurat_comb_singlet <- RunPCA(seurat_comb_singlet)
# seurat_comb_singlet <- RunHarmony(seurat_comb_singlet, reduction='pca', group.by.vars='orig.ident', reduction.save='harmony')
# pdf(file = paste0(workdir, "Results/Figures/2_pca_reduce.pdf"), width = 8, height = 6)
# # print(DimPlot(seurat_merge4, reduction = "umap", pt.size = 0.2, group.by = "cellTypes", label = TRUE,
#  #cells.highlight = list("HBCP4B-tumor"=rownames(seurat_Epithe@meta.data[seurat_Epithe@meta.data$cellSubtype_infercnv %in% c("CNV-1","CNV-2",) , ]), 
#  #cols.highlight = list("HBCP4B-tumor"= 'black'))))
# print(DimPlot(seurat_comb_singlet, reduction = "pca", group.by = "orig.ident"))
# print(DimPlot(seurat_comb_singlet, reduction = "pca", group.by = "group"))
# print(DimPlot(seurat_comb_singlet, reduction = "pca", group.by = "Patients"))
# print(DimHeatmap(seurat_comb_singlet, dims = 1:10, cells = 500, balanced = TRUE))
# dev.off()

## umap and tsne reduction
seurat_comb_singlet <- RunUMAP(seurat_comb_singlet, dims = 1:dim.usage, reduction.use='pca')
# seurat_comb_singlet <- RunUMAP(seurat_comb_singlet, reduction.use='pca', dims=c(3,4,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30))
# pdf(file = paste0(workdir, "Results/Figures/pca_reduce.pdf"), width = 8, height = 6)
# print(FeaturePlot(seurat_comb_singlet, features='PC_1'))
# print(FeaturePlot(seurat_comb_singlet, features='PC_2'))
# print(FeaturePlot(seurat_comb_singlet, features='PC_3'))
# print(FeaturePlot(seurat_comb_singlet, features='PC_4'))
# print(FeaturePlot(seurat_comb_singlet, features='PC_5'))
# print(FeaturePlot(seurat_comb_singlet, features='PC_6'))
# print(FeaturePlot(seurat_comb_singlet, features='PC_7'))
# print(FeaturePlot(seurat_comb_singlet, features='PC_8'))
# print(FeaturePlot(seurat_comb_singlet, features='PC_9'))
# print(FeaturePlot(seurat_comb_singlet, features='PC_10'))
# dev.off()

# pdf(file = paste0(workdir, "Results/Figures/3_RNA_featureplots.pdf"), width = 11, height = 5)
# FeaturePlot(seurat_comb_singlet, reduction = "umap", features = c("nFeature_RNA", "nCount_RNA"))
# dev.off()
pdf(file = paste0(workdir, "Results/Figures/3_umap.pdf"), width = 8, height = 4)
print(DimPlot(seurat_comb_singlet, reduction = "umap", pt.size = 0.2, group.by = "orig.ident", label = TRUE))
print(DimPlot(seurat_comb_singlet, reduction = "umap", pt.size = 0.2, group.by = "group", label = TRUE))
print(DimPlot(seurat_comb_singlet, reduction = "umap", pt.size = 0.05, group.by = "Patients", label = TRUE))
dev.off()
# save(seurat_comb_singlet, file = paste0(workdir, "saveData/seurat_", SampleID, "_umap.RData"))
# seurat_comb_singlet <- RunTSNE(seurat_comb_singlet, reduction='harmony', dims = 1:dim.usage, check_duplicates = FALSE, reduction.name='tsne')
# pdf(file = paste0(workdir, "Results/Figures/3_tsne.pdf"), width = 5, height = 4)
# # print(DimPlot(seurat_comb_singlet, reduction = "tsne", pt.size = 0.2, group.by = "orig.ident", label = TRUE))
# print(DimPlot(seurat_comb_singlet, reduction = "tsne", pt.size = 0.05, group.by = "Patients", label = FALSE))
# dev.off()

# #### cluster all the cells  ####
seurat_comb_singlet <- FindNeighbors(seurat_comb_singlet, reduction.use='pca')
seurat_comb_singlet <- FindClusters(seurat_comb_singlet, resolution = 0.2)
seurat_comb_singlet <- FindClusters(seurat_comb_singlet, resolution = 0.4)
seurat_comb_singlet <- FindClusters(seurat_comb_singlet, resolution = 0.5)
seurat_comb_singlet <- FindClusters(seurat_comb_singlet, resolution = 0.6)
seurat_comb_singlet <- FindClusters(seurat_comb_singlet, resolution = 0.8)
# seurat_comb_singlet <- FindClusters(seurat_comb_singlet, resolution = 1)
# seurat_comb_singlet <- FindClusters(seurat_comb_singlet, resolution = 1.2)

pdf(file = paste0(workdir, "Results/Figures/3_umap_cluster.pdf"), width = 7, height = 6)
print(DimPlot(seurat_comb_singlet, reduction = "umap", pt.size = 0.2, group.by = "RNA_snn_res.0.2", label = TRUE))
print(DimPlot(seurat_comb_singlet, reduction = "umap", pt.size = 0.2, group.by = "RNA_snn_res.0.4", label = TRUE))
print(DimPlot(seurat_comb_singlet, reduction = "umap", pt.size = 0.2, group.by = "RNA_snn_res.0.5", label = TRUE))
print(DimPlot(seurat_comb_singlet, reduction = "umap", pt.size = 0.2, group.by = "RNA_snn_res.0.6", label = TRUE))
print(DimPlot(seurat_comb_singlet, reduction = "umap", pt.size = 0.2, group.by = "RNA_snn_res.0.8", label = TRUE))
# print(DimPlot(seurat_comb_singlet, reduction = "umap", pt.size = 0.2, group.by = "RNA_snn_res.1", label = TRUE))
# print(DimPlot(seurat_comb_singlet, reduction = "umap", pt.size = 0.2, group.by = "RNA_snn_res.1.2", label = TRUE))
dev.off()

save(seurat_comb_singlet, file = paste0(workdir, "saveData/seurat_", SampleID, "_cluster.RData"))

pdf(file = paste0(workdir, "/Results/Figures/umap_all_immune_cluster_singlet.pdf"), width = 7, height = 6)
print(DimPlot(seurat_comb_singlet, reduction = "umap", pt.size = 0.2, group.by = 'cellTypes_new'))
dev.off()

### cluster find DEGs, and marker genes plots  ####
# Idents(seurat_comb_singlet) <- seurat_comb_singlet$RNA_snn_res.0.8
# seurat_comb_singlet_markers <- FindAllMarkers(seurat_comb_singlet, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# save(seurat_comb_singlet_markers, file = paste0(workdir, "saveData/seurat_comb_markers.RData"))
# seurat_comb_singlet_markers %>%
#     group_by(cluster) %>%
#     top_n(n = 30, wt = avg_log2FC) -> top30

# height <- 0.8 + nrow(top30) * 0.11
# pdf(file = paste0(workdir, "Results/Figures/4_ClusterMarkers_heatmap.pdf"), width = 25, height = height)
# print(DoHeatmap(seurat_comb_singlet, features = top30$gene) + NoLegend())
# dev.off()

# heig <- ceiling(length(top30$gene) / 4) * 3
# pdf(file = paste0(workdir, "Results/Figures/4_ClusterMarkers_featureplots.pdf"), width = 25, height = heig)
# print(FeaturePlot(seurat_comb_singlet, features = unique(top30$gene), reduction = "umap", ncol = 4))
# dev.off()

#### correlation of cluter
# table(seurat_comb_singlet$RNA_snn_res.1)
# av <- AverageExpression(seurat_comb_singlet, group.by = 'RNA_snn_res.1', assays = 'RNA')
# av <- av[[1]]
# head(av)

# #选出标准差最大的1000个基因
# cg <- names(tail(sort(apply(av,1,sd)),1000))
# #查看这1000个基因在各细胞群中的表达矩阵
# View(av[cg,])
# #查看细胞群的相关性矩阵
# View(cor(av[cg,],method='spearman'))
# # heatmap
# pdf(file = paste0(workdir, "Results/Figures/cluster_correlation.pdf"), width = 10, height = 10)
# pheatmap::pheatmap(cor(av[cg,],method='spearman'))
# dev.off()
