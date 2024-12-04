### Seurat clustering
library(Seurat)
library(dplyr)
library(patchwork)
library(DoubletFinder)
library(ggplot2)
library(future)

library(Cairo)
options(bitmapType = "cairo")

maxGenes <- 8000
SampleID <- "BRCA"
minUMIs <- 700
minGenes <- 200
maxPercent.mt <-10
dim.usage <- 20
res.usage <- 1.0
doublets.percentage <- 0.075
workdir <- "/jdfsbjcas1/ST_BJ/P21H28400N0232/wangdi/xieguixiang/BRCAsc_Bcell_plasmocyte/"
organs <- "BRCA"

source("/jdfsbjcas1/ST_BJ/P21H28400N0232/wangdi/tumors_public/codes/scRNA_primary.R")
source("/jdfsbjcas1/ST_BJ/P21H28400N0232/wangdi/tumors_public/codes/colors.R")
source("/jdfsbjcas1/ST_BJ/P21H28400N0232/wangdi/tumors_public/codes/signatureGenes.R")
source("/jdfsbjcas1/ST_BJ/P21H28400N0232/wangdi/tumors_public/codes/0_initNEWanalysis.R")
source("/jdfsbjcas1/ST_BJ/P21H28400N0232/wangdi/tumors_public/codes/DEG_wilcox.R")

setwd(workdir)
createObjdir(workdir = workdir)

plan()
plan("multiprocess", workers = 10)
plan()
options(future.globals.maxSize = 429496729600)

sc1 <- load("/jdfsbjcas1/ST_BJ/P21H28400N0232/wangdi/xieguixiang/BRCAsc_Bcell_plasmocyte/1.3.1_seurat_Bcell_batch.RData")
bcell <- get(sc1)
str(bcell)
sc2 <- load("/jdfsbjcas1/ST_BJ/P21H28400N0232/wangdi/xieguixiang/BRCAsc_Bcell_plasmocyte/1.3.1_seurat_plasmocyte_batch.RData")
plasmocyte <- get(sc2)
str(plasmocyte)

seurat_comb <- merge(bcell, y=plasmocyte, add.cell.ids = c("Bcell", "plasmocyte"))
str(seurat_comb)

#### plot the qc plots ####
VlnPlot(seurat_comb, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, group.by = "orig.ident")
plot1 <- FeatureScatter(seurat_comb, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(seurat_comb, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

plot3 <- FeatureScatter(seurat_comb, feature1 = "nCount_RNA", feature2 = "percent.mt", group.by = "group")
plot4 <- FeatureScatter(seurat_comb, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = "group")

if (!is.null(dev.list())) {
    dev.off()
}

pdf(file = paste0(workdir, "Results/Figures/1_qc.pdf"), width = 10, height = 6)
print(VlnPlot(seurat_comb, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, group.by = "group"))
print(VlnPlot(seurat_comb, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, group.by = NULL))

print(plot1 + plot2)
print(plot3 + plot4)

dev.off()

qc_stat <- rbind(summary(seurat_comb@meta.data$nCount_RNA), summary(seurat_comb@meta.data$nFeature_RNA), summary(seurat_comb@meta.data$percent.mt))
rownames(qc_stat) <- c("nUMIs", "nGenes", "percent.mt")
write.table(qc_stat, file = paste0(workdir, "/Results/", SampleID, "_raw_qc_stat.txt"), quote = F, sep = "\t", row.names = T, col.names = NA)


#### filter the data  ####
print(minUMIs)
print(minGenes)
print(maxGenes)
print(maxPercent.mt)
seurat_comb_filter <- subset(seurat_comb, subset = nCount_RNA > minUMIs & nFeature_RNA > minGenes & nFeature_RNA < maxGenes & percent.mt < maxPercent.mt)

cat("2.Filter Low quality cell matrix:", dim(seurat_comb_filter), "\n")

#### plot the plots of post plot  ####
plot1 <- FeatureScatter(seurat_comb_filter, feature1 = "nCount_RNA", feature2 = "percent.mt", group.by = "orig.ident")
plot2 <- FeatureScatter(seurat_comb_filter, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = "orig.ident")

plot3 <- FeatureScatter(seurat_comb_filter, feature1 = "nCount_RNA", feature2 = "percent.mt", group.by = "group")
plot4 <- FeatureScatter(seurat_comb_filter, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = "group")

pdf(file = paste0(workdir, "Results/Figures/1_qc_post.pdf"), width = 10, height = 6)
print(VlnPlot(seurat_comb_filter, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, group.by = "group"))
print(VlnPlot(seurat_comb_filter, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, group.by = "orig.ident"))

print(plot1 + plot2)
print(plot3 + plot4)

dev.off()

#### statistic the QC results  ####
qc <- c(
    ncol(seurat_comb), nrow(seurat_comb), median(seurat_comb$nFeature_RNA), median(seurat_comb$nCount_RNA), median(seurat_comb$percent.mt),
    mean(seurat_comb$nFeature_RNA), mean(seurat_comb$nCount_RNA), mean(seurat_comb$percent.mt)
)
qcPost <- c(
    ncol(seurat_comb_filter), nrow(seurat_comb_filter), median(seurat_comb_filter$nFeature_RNA), median(seurat_comb_filter$nCount_RNA),
    median(seurat_comb_filter$percent.mt), mean(seurat_comb_filter$nFeature_RNA), mean(seurat_comb_filter$nCount_RNA), mean(seurat_comb_filter$percent.mt)
)
stat <- data.frame(qcBefore = qc, qcPost = qcPost)
rownames(stat) <- c("cells", "genes", "medianFeatures", "medianCounts", "medianMT", "meanFeatures", "meanCounts", "meanMT")
stat$percent <- 100 * qcPost / qc
write.csv(stat, file = paste0(workdir, "Results/QCstat.csv"), quote = FALSE)

save(seurat_comb, file = paste0(workdir, "saveData/seurat_", SampleID, ".RData"))
rm(seurat_comb)

qc_stat <- rbind(summary(seurat_comb_filter@meta.data$nCount_RNA), summary(seurat_comb_filter@meta.data$nFeature_RNA), summary(seurat_comb_filter@meta.data$percent.mt))
rownames(qc_stat) <- c("nUMIs", "nGenes", "percent.mt")
write.table(qc_stat, file = paste0(workdir, "/Results/", SampleID, "post_qc_stat.txt"), quote = F, sep = "\t", row.names = T, col.names = NA)

#### normalize and reduce the data ####
seurat_comb_filter <- NormalizeData(seurat_comb_filter)
seurat_comb_filter <- FindVariableFeatures(seurat_comb_filter, selection.method = "vst", nfeatures = 2000)

top10genes <- head(VariableFeatures(seurat_comb_filter), 10)
plot1 <- VariableFeaturePlot(seurat_comb_filter)
plot2 <- LabelPoints(plot = plot1, points = top10genes, repel = TRUE)

pdf(file = paste0(workdir, "Results/Figures/2_vstgenes.pdf"), width = 18, height = 6)
print(plot1 + plot2)
dev.off()

all.genes <- rownames(seurat_comb_filter)
seurat_comb_filter <- ScaleData(seurat_comb_filter, features = all.genes)

seurat_comb_filter <- RunPCA(seurat_comb_filter, features = VariableFeatures(object = seurat_comb_filter))


pdf(file = paste0(workdir, "Results/Figures/2_pca_reduce.pdf"), width = 8, height = 6)
print(DimPlot(seurat_comb_filter, reduction = "pca", group.by = "orig.ident"))
print(DimPlot(seurat_comb_filter, reduction = "pca", group.by = "group"))
print(DimHeatmap(seurat_comb_filter, dims = 1:10, cells = 500, balanced = TRUE))
dev.off()

## umap and tsne reduction
seurat_comb_filter <- RunUMAP(seurat_comb_filter, dims = 1:dim.usage)
seurat_comb_filter <- RunTSNE(seurat_comb_filter, dims = 1:dim.usage)

pdf(file = paste0(workdir, "Results/Figures/3_umap_tsne.pdf"), width = 5, height = 4)
print(DimPlot(seurat_comb_filter, reduction = "umap", pt.size = 0.2, group.by = "orig.ident", label = TRUE))
print(DimPlot(seurat_comb_filter, reduction = "umap", pt.size = 0.2, group.by = "group", label = TRUE))
print(DimPlot(seurat_comb_filter, reduction = "tsne", pt.size = 0.2, group.by = "orig.ident", label = TRUE))
print(DimPlot(seurat_comb_filter, reduction = "tsne", pt.size = 0.2, group.by = "group", label = TRUE))
dev.off()

pdf(file = paste0(workdir, "Results/Figures/3_RNA_featureplots.pdf"), width = 11, height = 5)
FeaturePlot(seurat_comb_filter, features = c("nFeature_RNA", "nCount_RNA"))
dev.off()

cat("3. reduce thedata:", dim(seurat_comb_filter), "\n")

#### to remove the doublets  ####

pdf(paste0(workdir, "Results/Figures/5_find_doublet.pdf"))
seurat_comb_filter <- Find_doublet(seurat_comb_filter, doublets.percentage = doublets.percentage)
DimPlot(seurat_comb_filter, reduction = "umap", group.by = "doublet_info")
dev.off()

write.table(seurat_comb_filter@meta.data, paste0(workdir, "Results/", SampleID, "_doublets_info.txt"), quote = F, sep = "\t", row.names = T, col.names = NA)

cat("4. remove the doublets ", table(seurat_comb_filter@meta.data$doublet_info), "\n")
#### to analysis with the singlets ####
seurat_comb_singlet <- subset(seurat_comb_filter, subset = doublet_info == "Singlet")

save(seurat_comb_filter, file = paste0(workdir, "saveData/seurat_", SampleID, "_filter.RData"))
rm(seurat_comb_filter)

#### normalize and reduce the data ####
seurat_comb_singlet <- NormalizeData(seurat_comb_singlet)
seurat_comb_singlet <- FindVariableFeatures(seurat_comb_singlet, selection.method = "vst", nfeatures = 2000)

top10genes <- head(VariableFeatures(seurat_comb_singlet), 10)
plot1 <- VariableFeaturePlot(seurat_comb_singlet)
plot2 <- LabelPoints(plot = plot1, points = top10genes, repel = TRUE)

pdf(file = paste0(workdir, "Results/Figures/2_vstgenes_singlet.pdf"), width = 18, height = 6)
print(plot1 + plot2)
dev.off()

all.genes <- rownames(seurat_comb_singlet)
seurat_comb_singlet <- ScaleData(seurat_comb_singlet, features = all.genes)

seurat_comb_singlet <- RunPCA(seurat_comb_singlet, features = VariableFeatures(object = seurat_comb_singlet))


pdf(file = paste0(workdir, "Results/Figures/2_pca_reduce_singlet.pdf"), width = 8, height = 6)
# print(DimPlot(seurat_merge4, reduction = "umap", pt.size = 0.2, group.by = "cellTypes", label = TRUE,
 #cells.highlight = list("HBCP4B-tumor"=rownames(seurat_Epithe@meta.data[seurat_Epithe@meta.data$cellSubtype_infercnv %in% c("CNV-1","CNV-2",) , ]), 
 #cols.highlight = list("HBCP4B-tumor"= 'black'))))

print(DimPlot(seurat_comb_singlet, reduction = "pca", group.by = "orig.ident"))
print(DimPlot(seurat_comb_singlet, reduction = "pca", group.by = "group"))
print(DimHeatmap(seurat_comb_singlet, dims = 1:10, cells = 500, balanced = TRUE))
dev.off()

## umap and tsne reduction
seurat_comb_singlet <- RunUMAP(seurat_comb_singlet, dims = 1:dim.usage)
seurat_comb_singlet <- RunTSNE(seurat_comb_singlet, dims = 1:dim.usage)

pdf(file = paste0(workdir, "Results/Figures/3_RNA_featureplots_singlet.pdf"), width = 11, height = 5)
FeaturePlot(seurat_comb_singlet, reduction = "umap", features = c("nFeature_RNA", "nCount_RNA"))
dev.off()

pdf(file = paste0(workdir, "Results/Figures/3_umap_tsne_singlet.pdf"), width = 5, height = 4)
print(DimPlot(seurat_comb_singlet, reduction = "umap", pt.size = 0.2, group.by = "orig.ident", label = TRUE))
print(DimPlot(seurat_comb_singlet, reduction = "umap", pt.size = 0.2, group.by = "group", label = TRUE))
print(DimPlot(seurat_comb_singlet, reduction = "tsne", pt.size = 0.2, group.by = "orig.ident", label = TRUE))
print(DimPlot(seurat_comb_singlet, reduction = "tsne", pt.size = 0.2, group.by = "group", label = TRUE))
dev.off()

cat("5. reduce the singlets ", table(seurat_comb_singlet@meta.data$doublet_info), "\n")
#### cluster all the cells  ####
seurat_comb_singlet <- FindNeighbors(seurat_comb_singlet, dims = 1:dim.usage)

seurat_comb_singlet <- FindClusters(seurat_comb_singlet, resolution = 0.4)
seurat_comb_singlet <- FindClusters(seurat_comb_singlet, resolution = 0.5)
seurat_comb_singlet <- FindClusters(seurat_comb_singlet, resolution = 0.6)
seurat_comb_singlet <- FindClusters(seurat_comb_singlet, resolution = 0.8)
seurat_comb_singlet <- FindClusters(seurat_comb_singlet, resolution = 1.2)
seurat_comb_singlet <- FindClusters(seurat_comb_singlet, resolution = 1)

pdf(file = paste0(workdir, "Results/Figures/3_umap_cluster_singlet.pdf"), width = 7, height = 6)
print(DimPlot(seurat_comb_singlet, reduction = "umap", pt.size = 0.2, group.by = "RNA_snn_res.0.4", label = TRUE))
print(DimPlot(seurat_comb_singlet, reduction = "umap", pt.size = 0.2, group.by = "RNA_snn_res.0.5", label = TRUE))
print(DimPlot(seurat_comb_singlet, reduction = "umap", pt.size = 0.2, group.by = "RNA_snn_res.0.6", label = TRUE))
print(DimPlot(seurat_comb_singlet, reduction = "umap", pt.size = 0.2, group.by = "RNA_snn_res.0.8", label = TRUE))
print(DimPlot(seurat_comb_singlet, reduction = "umap", pt.size = 0.2, group.by = "RNA_snn_res.1", label = TRUE))
print(DimPlot(seurat_comb_singlet, reduction = "umap", pt.size = 0.2, group.by = "RNA_snn_res.1.2", label = TRUE))
dev.off()

cat("6. cluster the singlets ", table(seurat_comb_singlet@meta.data$doublet_info), "\n")
save(seurat_comb_singlet, file = paste0(workdir, "saveData/seurat_", SampleID, "_single.RData"))

#### cluster find DEGs, and marker genes plots  ####

SigGeneral_filt <- SigGeneral[SigGeneral %in% rownames(seurat_comb_singlet)]

seurat_comb_singlet_markers <- FindAllMarkers(seurat_comb_singlet, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

save(seurat_comb_singlet_markers, file = paste0(workdir, "saveData/seurat_comb_singlet_markers.RData"))
seurat_comb_singlet_markers %>%
    group_by(cluster) %>%
    top_n(n = 10, wt = avg_log2FC) -> top10
markers_PDAC <- seurat_comb_singlet_markers[seurat_comb_singlet_markers$gene %in% unlist(PDAC_marker), ]
height <- 0.8 + nrow(top10) * 0.11
pdf(file = paste0(workdir, "Results/Figures/4_ClusterMarkers_heatmap_singlet.pdf"), width = 10, height = height)
print(DoHeatmap(seurat_comb_singlet, features = top10$gene) + NoLegend())
dev.off()
pdf(file = paste0(workdir, "Results/Figures/4_ClusterMarkers_heatmap_singlet_PDAC.pdf"), width = 10, height = height)
print(DoHeatmap(seurat_comb_singlet, features = markers_PDAC$gene) + NoLegend())
dev.off()
hei <- ceiling(length(SigGeneral_filt) / 4) * 3
pdf(file = paste0(workdir, "Results/Figures/4_SigGeneral_featureplots_singlet.pdf"), width = 14.5, height = hei)
print(FeaturePlot(seurat_comb_singlet, features = SigGeneral_filt, reduction = "umap", ncol = 4))
dev.off()

heig <- ceiling(length(top10$gene) / 4) * 3
pdf(file = paste0(workdir, "Results/Figures/4_ClusterMarkers_featureplots_singlet.pdf"), width = 14.5, height = heig)
print(FeaturePlot(seurat_comb_singlet, features = unique(top10$gene), reduction = "umap", ncol = 4))
dev.off()

strommarker <- unique(SigGeneral[40:81])
strommarker <- strommarker[strommarker %in% rownames(seurat_comb_singlet)]

immunemarker <- unique(SigGeneral[1:39])
immunemarker <- immunemarker[immunemarker %in% rownames(seurat_comb_singlet)]

if (organs == "BRCA") {
    genemarker <- BRCA_marker
} else if (organs == "PDAC") {
    genemarker <- PDAC_marker
}
genemarker <- unique(unlist(genemarker))
genemarker <- genemarker[genemarker %in% rownames(seurat_comb_singlet)]

hei <- ceiling(length(strommarker) / 4) * 3
pdf(file = paste0(workdir, "Results/Figures/4_strommarker_featureplots.pdf"), width = 14.5, height = hei)
print(FeaturePlot(seurat_comb_singlet, features = strommarker, reduction = "umap", ncol = 4))
dev.off()
hei <- ceiling(length(immunemarker) / 4) * 3
pdf(file = paste0(workdir, "Results/Figures/4_immunemarker_featureplots.pdf"), width = 14.5, height = hei)
print(FeaturePlot(seurat_comb_singlet, features = immunemarker, reduction = "umap", ncol = 4))
dev.off()
hei <- ceiling(length(genemarker) / 4) * 3
pdf(file = paste0(workdir, "Results/Figures/4_genemarker_featureplots.pdf"), width = 14.5, height = hei)
print(FeaturePlot(seurat_comb_singlet, features = genemarker, reduction = "umap", ncol = 4))
dev.off()

wid <- length(unique(Idents(seurat_comb_singlet))) * 0.3 + 1

hei <- length(strommarker) * 0.3 + 1
pdf(file = paste0(workdir, "Results/Figures/4_strommarker_vlnplots.pdf"), width = wid, height = hei)
print(VlnPlot(seurat_comb_singlet, features = strommarker, log = F, stack = TRUE, flip = TRUE) + NoLegend())
dev.off()
hei <- length(immunemarker) * 0.3 + 1
pdf(file = paste0(workdir, "Results/Figures/4_immunemarker_vlnplots.pdf"), width = wid, height = hei)
print(VlnPlot(seurat_comb_singlet, features = immunemarker, log = F, stack = TRUE, flip = TRUE) + NoLegend())
dev.off()
hei <- length(genemarker) * 0.3 + 1
pdf(file = paste0(workdir, "Results/Figures/4_genemarker_vlnplots.pdf"), width = wid, height = hei)
print(VlnPlot(seurat_comb_singlet, features = genemarker, log = F, stack = TRUE, flip = TRUE) + NoLegend())
dev.off()
if (FALSE) {
    cat("7. cluster cell type enrichment", table(seurat_comb_singlet@meta.data$doublet_info), "\n")

    t1 <- Sys.time()
    seurat_comb_singlet_sigmarker_av <- Seurat_cluster_wilcoxRankEnrichment(seurat_obj = seurat_comb_singlet, genesets_len = 1, rank_meth = "averag", workdir = workdir)
    t2 <- Sys.time()
    t2 - t1
    seurat_comb_singlet_sigmarker_fc <- Seurat_cluster_wilcoxRankEnrichment(seurat_obj = seurat_comb_singlet, genesets_len = 1, rank_meth = "FCrank", workdir = workdir)
    t3 <- Sys.time()
    t3 - t2
    save(seurat_comb_singlet_sigmarker_av, file = paste0(workdir, "saveData/seurat_comb_singlet_sigmarker_av.RData"))
    save(seurat_comb_singlet_sigmarker_fc, file = paste0(workdir, "saveData/seurat_comb_singlet_sigmarker_fc.RData"))

    pdf(paste0(workdir, "/Results/Figures/6_celltype_averag_enrichment.pdf"))
    ggplot(seurat_comb_singlet_sigmarker_av$Enrichment_all, aes(clust, Description)) +
        geom_point(aes(size = NES, fill = p.Val), alpha = 0.9, pch = 21, col = "grey25") +
        theme_classic() +
        scale_fill_gradient2(low = "white", high = "red") +
        theme(
            axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, colour = "black"),
            axis.text.y = element_text(colour = "black")
        )
    dev.off()

    pdf(paste0(workdir, "/Results/Figures/6_celltype_FCrank_enrichment.pdf"))
    ggplot(seurat_comb_singlet_sigmarker_fc$Enrichment_all, aes(clust, Description)) +
        geom_point(aes(size = NES, fill = p.Val), alpha = 0.9, pch = 21, col = "grey25") +
        theme_classic() +
        scale_fill_gradient2(low = "white", high = "red") +
        theme(
            axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, colour = "black"),
            axis.text.y = element_text(colour = "black")
        )
    dev.off()


    cat("8. finished", table(seurat_comb_singlet@meta.data$doublet_info), "\n")
}
#### to enrichment the clusters
