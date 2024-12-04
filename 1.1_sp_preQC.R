################################################################################
## script to pre analysis and QC for spatial data, to chose the fit combined bins
## there are 2 model to operate the script
## 1. to run the script by Rscript, and set the parameters: full path of matrix data
##    bins, full path outdir: Rscript spatial_preQC input bins outdir
## 2. to run in command lines, set the argument on dommand lines
################################################################################

# setwd("/data/wangdi/PAAD_xiehe/ST/R4923_ST/")
library(Seurat)
library(dplyr)
library(data.table)
library(Matrix)
library(rjson)
library(ggplot2)
library(ggsci)
library(patchwork)
library(RColorBrewer)
library(pheatmap)
library(cowplot)
library(ggpubr)
library(SingleCellExperiment)
library(RcppML)
library(ggplot2)

library(Cairo)
options(bitmapType = "cairo")
options(RcppML.threads = 3)

source("/jdfsbjcas1/ST_BJ/P21H28400N0232/wangdi/tumors_public/codes/scRNA_primary.R")
source("/jdfsbjcas1/ST_BJ/P21H28400N0232/wangdi/tumors_public/codes/colors.R")
source("/jdfsbjcas1/ST_BJ/P21H28400N0232/wangdi/tumors_public/codes/signatureGenes.R")

library(getopt)
spec <- matrix(c(
    "help", "h", 0, "logical", "help document",
    "infile", "i", 1, "character", "input file,must bt gem file or gem.gz file",
    "outdir", "o", 1, "character", "output dierectory",
    "binsize", "b", 1, "integer", "BIN size for spatial rna-seq,default=50",
    "id", "d", 1, "character", "project id",
    "pc", "p", 1, "integer", "dimension usage in umap and tsne,default=30",
    "res", "r", 1, "double", "resolution of cluster,default=1.0",
    "tumor", "t", "1", "character", "datatype of sample,[PDAC,BRCA]",
    "nGenes", "n", "1", "character", "differential gene number,default=2000"
), byrow = TRUE, ncol = 5)
opt <- getopt(spec)
if (!is.null(opt$help) || is.null(opt$infile) || is.null(opt$outdir)) {
    cat(getopt(spec, usage = TRUE))
    q(status = 1)
}
infile <- opt$infile
outdir <- opt$outdir
if (is.null(opt$binsize)) {
    opt$binsize <- 50
}


if (is.null(opt$id)) {
    opt$id <- tail(unlist(strsplit(infile, "/", fixed = TRUE)), 1)
    opt$id <- gsub("mat.|.txt|.tsv|.gz|_filtered|.gem", "", opt$id)
}
if (is.null(opt$nGenes)) {
    opt$nGenes <- 2000
}
if (is.null(opt$pc)) {
    opt$pc <- 30
}
if (is.null(opt$res)) {
    opt$res <- 0.6
}
pc <- opt$pc
bs <- opt$binsize
# setwd("/data/public/gongchanghao/script/spatial_pipeline/test/")
organs <- opt$tumor

dir.create(outdir)
setwd(outdir)

pro <- opt$id

############################## 1. bin data  ##############################
dat <- fread(file = infile)
if (grep("MIDCount|MIDCounts", colnames(dat)) > 0) {
    colnames(dat) <- gsub("MIDCount|MIDCounts", "UMICount", colnames(dat))
}
out <- as.data.frame(dat)

dat$x <- trunc((dat$x - min(dat$x)) / bs + 1)
dat$y <- trunc((dat$y - min(dat$y)) / bs + 1)

out <- cbind(dat$y, dat$x, out)
colnames(out)[1:2] <- c(paste0("bin", bs, ".y"), paste0("bin", bs, ".x"))

fwrite(out, paste0(pro, "_bins", bs, "_information.txt"), col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)

dat <- dat[, sum(UMICount), by = .(geneID, x, y)]
dat$bin_ID <- max(dat$x) * (dat$y - 1) + dat$x
bin.coor <- dat[, sum(V1), by = .(x, y)]

out <- as.data.frame(cbind(paste0("BIN.", unique(dat$bin_ID)), bin.coor$y, bin.coor$x))
colnames(out) <- c(paste0("BIN.", bs), paste0("bin", bs, ".y"), paste0("bin", bs, ".x"))
rownames(out) <- out[, 1]
fwrite(out, paste0(pro, "_bin", bs, "_position.txt"),
    col.names = T, row.names = F, sep = "\t", quote = FALSE
)

##
geneID <- seq(length(unique(dat$geneID))) ## 36249 detected genes
hash.G <- data.frame(row.names = unique(dat$geneID), values = geneID)
gen <- hash.G[dat$geneID, "values"]

##
bin_ID <- unique(dat$bin_ID)
hash.B <- data.frame(row.names = sprintf("%d", bin_ID), values = bin_ID)
bin <- hash.B[sprintf("%d", dat$bin_ID), "values"]

##
cnt <- dat$V1

rm(dat)
gc()

##
tissue_lowres_image <- matrix(1, max(bin.coor$y), max(bin.coor$x))
tissue_positions_list <- data.frame(
    row.names = paste("BIN", rownames(hash.B), sep = "."),
    tissue = 1,
    row = bin.coor$y,
    col = bin.coor$x,
    imagerow = bin.coor$y,
    imagecol = bin.coor$x
)



scalefactors_json <- toJSON(list(
    fiducial_diameter_fullres = 1,
    tissue_hires_scalef = 1,
    tissue_lowres_scalef = 1
))




##
mat <- sparseMatrix(i = gen, j = bin, x = cnt)
rownames(mat) <- rownames(hash.G)
colnames(mat) <- paste("BIN", sprintf("%d", seq(max(hash.B[, "values"]))), sep = ".")


############################## 2. creat Spatial Object  ##############################
seurat_spatialObj <- CreateSeuratObject(mat, project = "Spatial", assay = "Spatial", min.cells = 1, min.features = 1)
generate_spatialObj <- function(image, scale.factors, tissue.positions, filter.matrix = TRUE) {
    if (filter.matrix) {
        tissue.positions <- tissue.positions[which(tissue.positions$tissue == 1), , drop = FALSE]
    }

    unnormalized.radius <- scale.factors$fiducial_diameter_fullres * scale.factors$tissue_lowres_scalef

    spot.radius <- unnormalized.radius / max(dim(image))

    return(new(
        Class = "VisiumV1",
        image = image,
        scale.factors = scalefactors(
            spot = scale.factors$tissue_hires_scalef,
            fiducial = scale.factors$fiducial_diameter_fullres,
            hires = scale.factors$tissue_hires_scalef,
            lowres = scale.factors$tissue_lowres_scalef
        ),
        coordinates = tissue.positions,
        spot.radius = spot.radius
    ))
}

spatialObj <- generate_spatialObj(
    image = tissue_lowres_image,
    scale.factors = fromJSON(scalefactors_json),
    tissue.positions = tissue_positions_list
)

##
spatialObj <- spatialObj[Cells(seurat_spatialObj)]
DefaultAssay(spatialObj) <- "Spatial"

seurat_spatialObj[["slice1"]] <- spatialObj

##
rm(mat)
rm(bin.coor)
rm(hash.G)
rm(hash.B)
rm(bin)
rm(gen)
rm(cnt)

dir.create("figures")

##############################  3. Spatial Analyse  ##############################
seurat_spatialObj[["percent.mt"]] <- PercentageFeatureSet(seurat_spatialObj, pattern = "^MT-")

##
Q1 <- quantile(seurat_spatialObj$nFeature_Spatial)[2]
Q3 <- quantile(seurat_spatialObj$nFeature_Spatial)[4]
upper <- as.numeric(Q3 + 1.5 * (Q3 - Q1))
lower <- as.numeric(Q1 - 1.5 * (Q3 - Q1))

save(list = ls(), file = "tmpfiles.RData")

pdf(paste0("figures/", pro, "_bin", bs, "_preQC.pdf"), width = 10, height = 8)
p1 <- VlnPlot(seurat_spatialObj,
    features = c("nFeature_Spatial", "nCount_Spatial", "percent.mt"), group.by = "orig.ident",
    ncol = 4, pt.size = 0
) +
    theme(axis.text.x = element_text(angle = 0, size = 0), axis.title.x = element_text(angle = 20, size = 8)) +
    labs(x = paste0("nGene:", dim(seurat_spatialObj)[1], "; ", "nBIN:", dim(seurat_spatialObj)[2]))
print(p1)

p2 <- ggplot(seurat_spatialObj@meta.data, aes(x = nFeature_Spatial)) +
    geom_density(colour = "black") +
    theme_classic() +
    theme(
        plot.title = element_text(hjust = 0.5, size = 18, face = "bold.italic"), legend.position = "none",
        axis.title = element_text(size = 15, face = "bold.italic"), axis.text.x = element_text(size = 12), axis.ticks.x = element_blank()
    ) +
    geom_vline(aes(xintercept = 100, colour = "#999999", linetype = "twodash")) +
    geom_vline(aes(xintercept = 200, colour = "#999999", linetype = "twodash")) +
    geom_vline(aes(xintercept = 300, colour = "#999999", linetype = "twodash")) +
    geom_vline(aes(xintercept = 500, colour = "#999999", linetype = "twodash")) +
    geom_vline(aes(xintercept = lower, colour = "#377EB8", linetype = "twodash")) +
    geom_vline(aes(xintercept = upper, colour = "#E41A1C", linetype = "twodash")) +
    xlim(min(seurat_spatialObj@meta.data$nFeature_Spatial), max(seurat_spatialObj@meta.data$nFeature_Spatial)) +
    ggtitle(paste0(pro, ".nBIN_", bs, ":", dim(seurat_spatialObj@meta.data)[1]))
print(p2)

dev.off()

pdf(paste0("figures/", pro, "_bin", bs, "_spatial_dis.pdf"))

# if(ncol(seurat_spatialObj) > 10000){
#   pt.s <- 1}else if(ncol(seurat_spatialObj) > 5000){
#     pt.s <- 1.2}else if(ncol(seurat_spatialObj) > 2000){
#       pt.s <- 1.4}else if(ncol(seurat_spatialObj) > 1000){
#         pt.s <- 2}else {pt.s <- 2.4
#         }

SpatialFeaturePlot(seurat_spatialObj, features = "nFeature_Spatial", stroke = 0) +
    theme(legend.position = "right")
SpatialFeaturePlot(seurat_spatialObj, features = "nCount_Spatial", stroke = 0) +
    theme(legend.position = "right")

dev.off()

seurat_spatialObj <- NormalizeData(seurat_spatialObj, assay = "Spatial", verbose = FALSE)
removeBiasGenes <- function(mat){

  RPgenes <- rownames(mat)[intersect(grep("^RP", rownames(mat)), grep("-", rownames(mat)))]
  RPgenes2 <- rownames(mat)[grep("^RP[SL]", rownames(mat))]
  MTgenes <- rownames(mat)[grep("^MT-", rownames(mat))]
  CTCgenes <- rownames(mat)[intersect(grep("^CTC", rownames(mat)), grep("-", rownames(mat)))]
  MIRgenes <- rownames(mat)[grep("^MIR", rownames(mat))]
  ACgenes <- rownames(mat)[intersect(grep("^AC[0-9]", rownames(mat)), grep(".", rownames(mat)))]
  CTgenes <- rownames(mat)[intersect(grep("^CT", rownames(mat)), grep("-", rownames(mat)))]
  LINCgenes <- rownames(mat)[grep("^LINC[0-9]", rownames(mat))]
  ALgenes <- rownames(mat)[intersect(grep("^AL", rownames(mat)), grep(".", rownames(mat)))]

  rmgenes <- c(RPgenes, RPgenes2, MTgenes, CTCgenes, MIRgenes, ACgenes, CTgenes, LINCgenes, ALgenes)

  # datacount <- mat[!rownames(mat)%in%rmgenes,]
  # datacount <- datacount[rowSums(datacount > 0) > 1,]
  return(rmgenes)
}
rmgenes <- removeBiasGenes(seurat_spatialObj)
filetered_genes <- rownames(seurat_spatialObj)[!rownames(seurat_spatialObj)%in%rmgenes]
seurat_spatialObj <- subset(seurat_spatialObj, features = filetered_genes)
seurat_spatialObj <- subset(seurat_spatialObj, subset = nFeature_Spatial > 100)
seurat_spatialObj <- FindVariableFeatures(seurat_spatialObj, selection.method = "vst", nfeatures = opt$nGenes)
seurat_spatialObj <- ScaleData(seurat_spatialObj)
seurat_spatialObj <- RunPCA(seurat_spatialObj, features = VariableFeatures(object = seurat_spatialObj), verbose = F)
save(seurat_spatialObj, file = "seurat_spatialObj.RData")
# seurat_spatialObj <- JackStraw(seurat_spatialObj, num.replicate = 100)
# seurat_spatialObj <- ScoreJackStraw(seurat_spatialObj, dims = 1:50)
# JackStrawPlot(seurat_spatialObj,dims = 1:20)
# ElbowPlot(seurat_spatialObj)

seurat_spatialObj <- RunUMAP(seurat_spatialObj, reduction = "pca", dims = 1:pc)
seurat_spatialObj <- FindNeighbors(seurat_spatialObj, reduction = "pca", dims = 1:pc)

## set different resolution
seurat_spatialObj <- FindClusters(seurat_spatialObj, verbose = FALSE, resolution = 0.4)
seurat_spatialObj <- FindClusters(seurat_spatialObj, verbose = FALSE, resolution = 0.8)
seurat_spatialObj <- FindClusters(seurat_spatialObj, verbose = FALSE, resolution = 1)
seurat_spatialObj <- FindClusters(seurat_spatialObj, verbose = FALSE, resolution = 1.2)
seurat_spatialObj <- FindClusters(seurat_spatialObj, verbose = FALSE, resolution = 0.6)
library(RColorBrewer)
col <- c(brewer.pal(9, "Set1")[c(1:5, 7, 6, 8, 9)], brewer.pal(9, "Pastel1")[c(7, 2, 3, 5, 1)], brewer.pal(12, "Set3")[c(-2, -4, -6, -9, -11)], brewer.pal(8, "Set2")[c(3, 4, 7)], brewer.pal(8, "Pastel2")[8], brewer.pal(12, "Paired")[c(5, 6, 9, 10, 12)])
if (length(levels(seurat_spatialObj)) > length(col)) {
    col <- colorRampPalette(col)(length(levels(seurat_spatialObj)))
} else {
    col <- col[1:length(levels(seurat_spatialObj))]
}

plot2 <- ElbowPlot(seurat_spatialObj, ndims = 50, reduction = "pca")

pt.s <- 0.2

pdf(paste0("figures/", pro, "_bin", bs, "_Spatial_res_UMAP.pdf"))
print(plot2)
DimPlot(seurat_spatialObj, reduction = "umap", label = TRUE, group.by = "Spatial_snn_res.0.4")
DimPlot(seurat_spatialObj, reduction = "umap", label = TRUE, group.by = "Spatial_snn_res.0.6")
DimPlot(seurat_spatialObj, reduction = "umap", label = TRUE, group.by = "Spatial_snn_res.0.8")
DimPlot(seurat_spatialObj, reduction = "umap", label = TRUE, group.by = "Spatial_snn_res.1")
DimPlot(seurat_spatialObj, reduction = "umap", label = TRUE, group.by = "Spatial_snn_res.1.2")
dev.off()

#Idents(seurat_spatialObj) <- seurat_spatialObj$Spatial_snn_res.1.2

pdf(paste0("figures/", pro, "_bin", bs, "_Spatial_res_spatial.pdf"))
SpatialDimPlot(seurat_spatialObj, label = FALSE, label.size = 3, stroke = 0, group.by = "Spatial_snn_res.0.4")
SpatialDimPlot(seurat_spatialObj, label = FALSE, label.size = 3, stroke = 0, group.by = "Spatial_snn_res.0.6")
SpatialDimPlot(seurat_spatialObj, label = FALSE, label.size = 3, stroke = 0, group.by = "Spatial_snn_res.0.8")
SpatialDimPlot(seurat_spatialObj, label = FALSE, label.size = 3, stroke = 0, group.by = "Spatial_snn_res.1")
SpatialDimPlot(seurat_spatialObj, label = FALSE, label.size = 3, stroke = 0, group.by = "Spatial_snn_res.1.2")
dev.off()

save(seurat_spatialObj, file = "seurat_spatialObj_umap.RData")

pdf(paste0("figures/", pro, "_bin", bs, "_Spatial_UMAP_split.pdf"), length(levels(seurat_spatialObj)) * 5, 6)
DimPlot(seurat_spatialObj, reduction = "umap", label = TRUE, split.by = "Spatial_snn_res.0.8", cols = col57)
dev.off()
save(seurat_spatialObj, file = "seurat_spatialObj_umap.RData")





if (organs == "BRCA") {
    genemarker <- BRCA_marker
} else if (organs == "PDAC") {
    genemarker <- PDAC_marker
}
genemarker <- unique(unlist(genemarker))
genemarker <- genemarker[genemarker %in% rownames(seurat_spatialObj)]
strommarker <- unique(SigGeneral[40:81])
strommarker <- strommarker[strommarker %in% rownames(seurat_spatialObj)]

immunemarker <- unique(SigGeneral[1:39])
immunemarker <- immunemarker[immunemarker %in% rownames(seurat_spatialObj)]
hei <- ceiling(length(strommarker) / 4) * 3
pdf(file = "figures/4_strommarker_featureplots.pdf", width = 14.5, height = hei)
print(FeaturePlot(seurat_spatialObj, features = strommarker, reduction = "umap", ncol = 4, raster = TRUE))
dev.off()
hei <- ceiling(length(immunemarker) / 4) * 3
pdf(file = "figures/4_immunemarker_featureplots.pdf", width = 14.5, height = hei)
print(FeaturePlot(seurat_spatialObj, features = immunemarker, reduction = "umap", ncol = 4, raster = TRUE))
dev.off()
hei <- ceiling(length(genemarker) / 4) * 3
pdf(file = "figures/4_genemarker_featureplots.pdf", width = 14.5, height = hei)
print(FeaturePlot(seurat_spatialObj, features = genemarker, reduction = "umap", ncol = 4, raster = TRUE))
dev.off()

# ## cluster find DEGs, classification and annotation
# normalgenes <- c(
#     "PTPRC", "CD163", "SLC11A1", "APOC1", "CD86", "CSF1R", "SLCO2B1", "CD68", "F13A1", "CD14", "AIF1", "CD80",
#     "FCER1G", "FCGR3A", "TYROBP", "PALLD", "ACTA2", "SULF1", "CTGF", "TAGLN", "TPM1", "GINS1", "THY1", "RBP1", "COL1A2",
#     "HSPG2", "LDB2", "GPR116", "PTPRB", "VWF", "DOCK9", "CDH5", "SELE", "KIF4A", "FANCI", "CHAF1A", "GTSE1", "ASPM",
#     "SPC25", "NCAPG2", "POLA2", "NCAPD3", "CD19", "CR2", "MS4A1", "CD79A", "CD79B", "BLNK", "GFAP", "BMPR1B", "CD44",
#     "SLC1A2", "AQP4", "S100B", "GJB7", "ALDH1L1", "ALDOC", "MLC1", "CD2", "CD3D", "CD3E", "CD3G", "CD4", "CD8A", "CD8B",
#     "GZMK", "MBP", "SOX10", "MOG", "CA2", "CNP", "RTN4", "PLP1", "PLP2", "OPALIN", "OMG", "OLIG1", "TNR", "ALCAM", "PLLP",
#     "ITGAM", "ITGAX", "CEACAM8", "FCGR3A", "NCAM1", "CD14", "CD33", "FCGR2A", "FCGR2B", "ITGB2", "SELL", "KIT", "ENPP3", "THBD", "CD1C",
#     "CDK4", "HIF1A", "MET", "PDGFRA", "PDGFRB", "VEGFA", "VEGFB", "VEGFC", "EGFR", "CHI3L1", "ELTD1", "MKI67", "DCX", "NCAM1", "NCAM2",
#     "NEUROD1", "NEUROD2", "ENO2", "RBFOX3", "RBFOX2", "RBFOX1", "MAP2", "TUBB3", "NEFM", "NEFH", "GAP43", "TH", "GAD1", "GAD2", "SLC17A7",
#     "SLC17A6", "FEV", "PET1", "SLC6A4", "DLG4", "SYP", "BSN"
# )
# normalgenes <- normalgenes[normalgenes %in% rownames(seurat_spatialObj)]
# # FeaturePlot(seurat_spatialObj_filter, features = c("PTPRC", "CD163", "SLC11A1", "APOC1", "CD86"), reduction = "tsne")

# seurat_spatialObj_markers <- FindAllMarkers(seurat_spatialObj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
# seurat_spatialObj_markers %>%
#     group_by(cluster) %>%
#     top_n(n = 10, wt = avg_log2FC) -> top10

# pdf(file = "./figures/seurat_spatialObj_markers_heatmap.pdf", width = 10, height = 24)
# DoHeatmap(seurat_spatialObj, features = top10$gene) + NoLegend()
# dev.off()

# featuregenes <- c(normalgenes, top10$gene)
# pdf(file = "./figures/normal_markers_featureplots.pdf", width = 12, height = 75) # width = 4*3, height = 2.5*(gene_num/4)
# FeaturePlot(seurat_spatialObj, features = normalgenes, reduction = "umap", ncol = 4,, raster = TRUE)
# dev.off()

# pdf(file = "./figures/markers_featureplots.pdf", width = 12, height = 120) # width = 4*3, height = 2.5*(gene_num/4)
# FeaturePlot(seurat_spatialObj, features = unique(top10$gene), reduction = "umap", ncol = 4, raster = TRUE)
# dev.off()


############  Cluster Ident
# if (ncol(seurat_spatialObj) > 10000) {
#     pt.s <- 1
# } else if (ncol(seurat_spatialObj) > 5000) {
#     pt.s <- 1.2
# } else if (ncol(seurat_spatialObj) > 2000) {
#     pt.s <- 1.4
# } else if (ncol(seurat_spatialObj) > 1000) {
#     pt.s <- 2
# } else {
#     pt.s <- 2.4
# }

# # if(ncol(seurat_spatialObj) > 3000){
# # 	pt.s <- 1.2}else {pt.s <- 1.3
# # }

# p <- SpatialDimPlot(seurat_spatialObj,
#     cells.highlight = CellsByIdentities(object = seurat_spatialObj), cols.highlight = c("#DE2D26", "grey90"),
#     facet.highlight = TRUE, stroke = 0, ncol = 5
# )
# # for(i in 1:length(levels(seurat_spatialObj))) {
# #   p[[i]] <- p[[i]] + scale_y_reverse()
# # }

# pdf(paste0("figures/", pro, "_bin", bs, "_Cluster_Ident.pdf"))
# # SpatialDimPlot(seurat_spatialObj, cells.highlight=CellsByIdentities(object=seurat_spatialObj),cols.highlight = c("#DE2D26", "grey90"), facet.highlight=TRUE, stroke=0, ncol=4) + scale_y_reverse()
# # cowplot::plot_grid(plotlist=p, ncol = 5)
# # print(p)
# for (i in 1:length(levels(seurat_spatialObj))) {
#     print(p[[i]])
# }
# dev.off()

# run NMF crossValidate to determine k
dir.create("figures/NMF", recursive = TRUE)
# remove cell with low number of genes
seurat_spatialObj_filter <- subset(seurat_spatialObj, subset = nFeature_Spatial > 100)
sce <- as.SingleCellExperiment(seurat_spatialObj_filter)
# seurat_spatialObj <- NormalizeData(seurat_spatialObj, assay = "Spatial", verbose = FALSE)
data <- assay(sce,"logcounts")
cv_mse <- crossValidate(data, method = "predict", k = 1:40, seed = 123)
cv_robust <- crossValidate(data, method = "robust", k = 1:40, seed = 123)
pdf("figures/NMF/NMF_CV.pdf", width = 10, height =7)
P1 <-plot(cv_mse)
P2 <-plot(cv_robust)
print(P1+P2)
dev.off()
