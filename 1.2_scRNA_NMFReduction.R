### Seurat clustering
library(Seurat)
library(tidyverse)
library(patchwork)
library(ggplot2)
library(future)
library(Cairo)
library(harmony)
library(NMF)
options(bitmapType = "cairo")

maxGenes <- 8000
SampleID <- "BLCA"
minUMIs <- 700
minGenes <- 200
maxPercent.mt <-10
dim.usage <- 20
res.usage <- 1.0
doublets.percentage <- 0.075
workdir <- "/jdfsbjcas1/ST_BJ/P21H28400N0232/fanfengpu/BC/SingleCell/R28BC/Muscle/Results/NMFreduction_8/"
organs <- "BLCA"

source("/jdfsbjcas1/ST_BJ/P21H28400N0232/wangdi/tumors_public/codes/scRNA_primary.R")
source("/jdfsbjcas1/ST_BJ/P21H28400N0232/wangdi/tumors_public/codes/colors.R")
source("/jdfsbjcas1/ST_BJ/P21H28400N0232/wangdi/tumors_public/codes/signatureGenes.R")
source("/jdfsbjcas1/ST_BJ/P21H28400N0232/wangdi/tumors_public/codes/0_initNEWanalysis.R")
source("/jdfsbjcas1/ST_BJ/P21H28400N0232/wangdi/tumors_public/codes/DEG_wilcox.R")

setwd(workdir)
createObjdir(workdir = workdir)

x <- load('/jdfsbjcas1/ST_BJ/P21H28400N0232/fanfengpu/BC/SingleCell/R28BC/Data/Mesenchymal/seurat_merge_Mesenchymal.RData')
seurat_comb_singlet <- get(x)
# seurat_comb_singlet <- subset(seurat_comb_singlet, cellTypes_new != 'Undefined')
seurat_comb_singlet <- subset(seurat_comb_singlet, cellTypes_new == 'Myocyte')
str(seurat_comb_singlet)

# seurat_comb_singlet <- readRDS('/jdfsbjcas1/ST_BJ/P21H28400N0232/xieguixiang/BLCA/TSanalysis/SCdata/tumor.rds')
# Idents(seurat_comb_singlet) <- seurat_comb_singlet@meta.data$TS
#### normalize and reduce the data ####
# seurat_comb_singlet <- subset(seurat_comb_singlet, Patients %in% c("HBCP1", "HBCP3A", "HBCP3B", "HBCP3C","HBCP5A","HBCP5B","HBCP6","HBCP7","HBCP8","HBCP9","HBCP10","HBCP11","HBCP12","HBCP13","HBCP14","HBCP15","HBCP17","HBCP18"))
# print(table(seurat_comb_singlet$Patients))
seurat_comb_singlet <- NormalizeData(seurat_comb_singlet) %>% FindVariableFeatures() %>% ScaleData(do.center = F)
vm <- seurat_comb_singlet@assays$RNA@scale.data
saveRDS(vm, file = 'vm.rds')
# vm <- readRDS("/jdfsbjcas1/ST_BJ/P21H28400N0232/xieguixiang/BLCA/SC/R28BC/fibroblast/Reduction_NMF/vm.rds")
res <- nmf(vm, 8, method='snmf/r')
save(res, file='nmf_res.rda')
# load(file='/jdfsbjcas1/ST_BJ/P21H28400N0232/xieguixiang/BLCA/SC/R28BC/fibroblast/Reduction_NMF/nmf_res.rda')

##### Extract cluster top loading features
# 每个因子提取30个
fs <- extractFeatures(res, 30L)
fs <- lapply(fs, function(x) rownames(res)[x])
fs <- do.call('rbind', fs)
rownames(fs) <- paste0('cluster', 1:8)
write.csv(t(fs), 'c_NMF_TopGenes.csv')
DT::datatable(t(fs))


#### 选择用于后续分析的因子，使用NMF运行的结果进行降维和聚类
s.f <- 1:8
# 降维
cell1 <- colnames(seurat_comb_singlet)
cell2 <- colnames(coef(res))
cells <- intersect(cell1, cell2)
seurat_comb_singlet <- seurat_comb_singlet[,cells]
seurat_comb_singlet <- RunPCA(seurat_comb_singlet, verbose=F)
seurat_comb_singlet@reductions$nmf <- seurat_comb_singlet@reductions$pca
seurat_comb_singlet@reductions$nmf@cell.embeddings <- t(coef(res)[,cells])
seurat_comb_singlet@reductions$nmf@feature.loadings <- basis(res)


# seurat_comb_singlet <- PercentageFeatureSet(seurat_comb_singlet, pattern = "^MT-", col.name = "percent.mt")# run sctransform
# seurat_comb_singlet <- SCTransform(seurat_comb_singlet, vars.to.regress = "percent.mt", verbose = FALSE)


# seurat_comb_singlet <- RunPCA(seurat_comb_singlet, features = VariableFeatures(object = seurat_comb_singlet))
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
seurat_comb_singlet <- RunUMAP(seurat_comb_singlet, dims = s.f, reduction='nmf')
save(seurat_comb_singlet, file = paste0(workdir, "saveData/seurat_", SampleID, "_umap.RData"))

# pdf(file = paste0(workdir, "Results/Figures/3_RNA_featureplots.pdf"), width = 11, height = 5)
# FeaturePlot(seurat_comb_singlet, reduction = "umap", features = c("nFeature_RNA", "nCount_RNA"))
# dev.off()
pdf(file = paste0(workdir, "Results/Figures/3_umap.pdf"), width = 8, height = 4)
# print(DimPlot(seurat_comb_singlet, reduction = "umap", pt.size = 0.2, group.by = "orig.ident", label = TRUE))
# print(DimPlot(seurat_comb_singlet, reduction = "umap", pt.size = 0.2, group.by = "group", label = TRUE))
print(DimPlot(seurat_comb_singlet, reduction = "umap", pt.size = 0.05, group.by = "Patients", label = TRUE))
dev.off()

# seurat_comb_singlet <- RunTSNE(seurat_comb_singlet, reduction='harmony', dims = 1:dim.usage, check_duplicates = FALSE, reduction.name='tsne')
# pdf(file = paste0(workdir, "Results/Figures/3_tsne.pdf"), width = 5, height = 4)
# # print(DimPlot(seurat_comb_singlet, reduction = "tsne", pt.size = 0.2, group.by = "orig.ident", label = TRUE))
# print(DimPlot(seurat_comb_singlet, reduction = "tsne", pt.size = 0.05, group.by = "Patients", label = FALSE))
# dev.off()

#### 基于NMF降维矩阵的聚类
seurat_comb_singlet <- FindNeighbors(seurat_comb_singlet, reduction='nmf', dims = s.f) %>% FindClusters()
#### 基于因子最大载荷分类
seurat_comb_singlet$cluster <- apply(NMF::coefficients(res)[s.f,], 2, which.max) 

save(seurat_comb_singlet, file = paste0(workdir, "saveData/seurat_", SampleID, "_cluster.RData"))


##### 降维聚类结果可视化
pdf(file = paste0(workdir, "Results/Figures/NMF_cluster.pdf"), width = 6, height = 5)
DimPlot(seurat_comb_singlet, group.by = 'cluster', label = T) + ggtitle('Clustered by max loading')
dev.off()

# pdf(file = paste0(workdir, "Results/Figures/featurePlot.pdf"), width = 13, height = 20)
# FeaturePlot(seurat_comb_singlet, features=c('CD34','ACTA2','ACTG2','TPM1','SULF1','CD74','IL6','LIF','GSN','CXCL12',
#                                             'MFAP5','TRPA1','CALML5','MMP19','NR4A2','CXCL14','WNT5A','KIF26B','KCNIP1','ERBB2'), ncol=3)
# dev.off()
#### cluster find DEGs, and marker genes plots  ####
# Idents(seurat_comb_singlet) <- seurat_comb_singlet$RNA_snn_res.0.2
# seurat_comb_singlet_markers <- FindAllMarkers(seurat_comb_singlet, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# save(seurat_comb_singlet_markers, file = paste0(workdir, "saveData/seurat_comb_markers.RData"))
# seurat_comb_singlet_markers %>%
#     group_by(cluster) %>%
#     top_n(n = 10, wt = avg_log2FC) -> top10

# height <- 0.8 + nrow(top10) * 0.11
# pdf(file = paste0(workdir, "Results/Figures/4_ClusterMarkers_heatmap.pdf"), width = 10, height = height)
# print(DoHeatmap(seurat_comb_singlet, features = top10$gene) + NoLegend())
# dev.off()

# heig <- ceiling(length(top10$gene) / 4) * 3
# pdf(file = paste0(workdir, "Results/Figures/4_ClusterMarkers_featureplots.pdf"), width = 14.5, height = heig)
# print(FeaturePlot(seurat_comb_singlet, features = unique(top10$gene), reduction = "umap", ncol = 4))
# dev.off()