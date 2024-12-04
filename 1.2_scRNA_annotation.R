# # dendritic cell (DC)
# # mesenchymal stem cell (msc)
SigGeneral_all_sort <- list(
  thymic_epithelial_cell = c("Epcam","Foxn1","Krt14", "Aire","Psmb11",
                             "H2-Ab1","Il7","Ptn","Ctsl","Perp","Krt18",
                             "Mdk","Sfn","Krt5","Dlk2","Krt1","Prss16",
                             "Foxi1","Myod1","Chrna1","Meurod1","Chga"),
  adipocyte = c("Lpl", "Adipoq", "Fabp4", "Cebpa", "Cfd"),
  DC = c("Itgax", "Clec9a", "Xcr1"),
  endothelial_cells = c("Cdh5", "Pecam1", "Egfl7", "Esam", "Lyve1", "Vim", "Col4a1", "Sparc"),
  erythrocyte = c('Hba-a1', 'Hba-a2', 'Hbb-bt'),
  ETP = c("Notch1"),
  fibroblast = c("Pdgfra", "Colec11", "Aldh1a2", "Gdf10", "Fbn1", "Pi16", "Sema3d", "Col1a1", "Col1a2", "Col3a1", "Dcn"),
  innate_T = c("Zbtb16", "Klrb1", "Nkg7", "Klrd1"),
  #interesting = c("Apoe", "Apod", "Igfbp4", "Igfbp5"),
  macrophage = c("Cd68", "C1qa", "C1qb", "C1qc", "Lyz2"),
  monocyte = c("Cd14"),
  msc = c("Pdgfrb", "Pdgfra", "Dlk1", "Col3a1"),
  neutrophil = c("S100a8", "S100a9", "Ltf", 'Cebpe', 'Cebpa', 'Etv6', 'Foxp1', 'Gfi1', 'Spi1', 'Stat3', 'Elane', 'Gstm1', 'Lcn2', 'Lyz2', 'Mpo', 'Prtn3'),
  plasma_cells = c("Jchain","Ighg1", "Ighm"),
  smooth_muscle_cells = c("Acta2", "Myh11", "Pdgfrb", "Rgs5"),
  VDJ = c("Rag1", "Rag2"),
  B_cells = c("Ms4a1", "Cd19", "Cd79a", "Cd79b"),
  cycling = c("Mki67", "Top2a"),
  nucleated_hematopoietic_cells = c("Cd45"),
  NK_cell = c("Cd27"),
  Tgd = c("Cd24")
)



SigGeneral_mouse <- c("Cd45","Cd27","Cd24","Epcam","Foxn1","Krt14","Aire","Psmb11","H2-Ab1","Il7","Ptn","Ctsl","Perp","Krt18","Mdk","Sfn","Krt5","Dlk2","Krt1","Prss16","Foxi1","Myod1","Chrna1","Meurod1","Chga","Lpl","Adipoq","Fabp4","Cebpa","Cfd","Itgax","Clec9a","Xcr1","Cdh5","Pecam1","Egfl7","Esam","Lyve1","Vim","Col4a1","Sparc",'Hba-a1','Hba-a2','Hbb-bt',"Notch1","Pdgfra","Colec11","Aldh1a2","Gdf10","Fbn1","Pi16","Sema3d","Col1a1","Col1a2","Col3a1","Dcn","Zbtb16","Klrb1","Nkg7","Klrd1","Apoe","Apod","Igfbp4","Igfbp5","Cd68","C1qa","C1qb","C1qc","Lyz2","Cd14","Pdgfrb","Pdgfra","Dlk1","Col3a1","S100a8","S100a9","Ltf",'Cebpe','Cebpa','Etv6','Foxp1','Gfi1','Spi1','Stat3','Elane','Gstm1','Lcn2','Lyz2','Mpo','Prtn3',"Jchain","Ighg1","Ighm","Acta2","Myh11","Pdgfrb","Rgs5","Rag1","Rag2","Ms4a1","Cd19","Cd79a","Cd79b","Mki67","Top2a")


library('Seurat')
library('cowplot')
library('parallel')
library('ggplot2')
library('dplyr')
library("SingleR")
options(bitmapType='cairo')

if(!dir.exists(paste0(workdir, "Results/Figures/6_cluster_anno/1_resolution_selection"))){dir.create(paste0(workdir, "Results/Figures/6_cluster_anno/1_resolution_selection"), recursive = TRUE)}
if(!dir.exists(paste0(workdir, "Results/Figures/6_cluster_anno/2_cluster_mapping"))){dir.create(paste0(workdir, "Results/Figures/6_cluster_anno/2_cluster_mapping"), recursive = TRUE)}
if(!dir.exists(paste0(workdir, "Results/Figures/6_cluster_anno/3_heatmap_check"))){dir.create(paste0(workdir, "Results/Figures/6_cluster_anno/3_heatmap_check"), recursive = TRUE)}
if(!dir.exists(paste0(workdir, "Results/Figures/6_cluster_anno/4_anno_result"))){dir.create(paste0(workdir, "Results/Figures/6_cluster_anno/4_anno_result"), recursive = TRUE)}
# SigGeneral_all_sort <- SigGeneral_mouse_sort
# print(SigGeneral_all_sort)

# SigGeneral_all_sort <- list(
#     mycaf = c("ACTA2", "VIM", "CCN2", "COL1A1", "COL5A1", "COL6A1", "TNC", "TGFB1", "THY1", "TAGLN", "COL12A1", "PDGFRB"), # TGFβ/SMAD2/3
#     icaf = c("IL1A", "IL1B", "IL6", "IL11", "LIF", "CLEC3B", "COL14A1", "GSN", "LY6C1", "CXCL12", "CXCL14", "PDGFRA"), ## IL-1/JAK-STAT3
#     apcaf = c("SLPI", "SAA3", "CD74", "H2-Ab1", "NKAIN4", "IRF5"),
#     other = c('CD10','CD73','MMP11','CD34','CD248','CD146','PDGFRA','MKI67','MME','APC','IDO1','PDPN','CDH11')
# )

#*1.resolution选择
res4=DimPlot(seurat_comb_singlet, reduction = "umap", pt.size = 0.2, group.by = "RNA_snn_res.0.4", label = TRUE)
res5=DimPlot(seurat_comb_singlet, reduction = "umap", pt.size = 0.2, group.by = "RNA_snn_res.0.5", label = TRUE)
res6=DimPlot(seurat_comb_singlet, reduction = "umap", pt.size = 0.2, group.by = "RNA_snn_res.0.6", label = TRUE)
res8=DimPlot(seurat_comb_singlet, reduction = "umap", pt.size = 0.2, group.by = "RNA_snn_res.0.8", label = TRUE)
res10=DimPlot(seurat_comb_singlet, reduction = "umap", pt.size = 0.2, group.by = "RNA_snn_res.1", label = TRUE)
ress = plot_grid(res4,res5,res6,res8,res10, ncol = 2)

for (celltype in names(SigGeneral_all_sort)) {
  print(celltype)
  genemarker = SigGeneral_all_sort[[celltype]]
  genemarker = genemarker[genemarker %in% rownames(seurat_comb_singlet)]
  if(length(genemarker)<1){next}
  width <- 18+ceiling(length(genemarker)/3) * 6
  umap_celltype = FeaturePlot(seurat_comb_singlet, features = genemarker, reduction = "umap", ncol = ceiling(length(genemarker)/3))
  pdf(file = paste0(workdir, "Results/Figures/6_cluster_anno/1_resolution_selection/",celltype,".pdf"), width=width, height=18)
  print(plot_grid(ress, umap_celltype))
  dev.off()
  tiff(file = paste0(workdir, "Results/Figures/6_cluster_anno/1_resolution_selection/",celltype,".tiff"), res=300, width=width, height=18, compression="lzw", units="in")
  print(plot_grid(ress, umap_celltype, rel_widths = c(18,ceiling(length(genemarker)/3) * 6)))
  dev.off()
}

#*2.cluster注释
step2 <- function(res){
  if(!dir.exists(paste0(workdir, "Results/Figures/6_cluster_anno/2_cluster_mapping/res.", res))){dir.create(paste0(workdir, "Results/Figures/6_cluster_anno/2_cluster_mapping/res.", res), recursive = TRUE)}
  res_i = get(paste0('res', res*10))
  if (res==0.4) {
    Idents(seurat_comb_singlet) = seurat_comb_singlet$RNA_snn_res.0.4
  } else if (res==0.5) {
    Idents(seurat_comb_singlet) = seurat_comb_singlet$RNA_snn_res.0.5
  } else if (res==0.6) {
    Idents(seurat_comb_singlet) = seurat_comb_singlet$RNA_snn_res.0.6
  } else if (res==0.8) {
    Idents(seurat_comb_singlet) = seurat_comb_singlet$RNA_snn_res.0.8
  } else if (res==1) {
    Idents(seurat_comb_singlet) = seurat_comb_singlet$RNA_snn_res.1
  }
  #pdfs = list()
  for (celltype in names(SigGeneral_all_sort)) {
    print(c(res,celltype))
    genemarker = SigGeneral_all_sort[[celltype]]
    genemarker = genemarker[genemarker %in% rownames(seurat_comb_singlet)]
    if(length(genemarker)<1){next}
    width <- 18+ceiling(length(genemarker)/3) * 6
    umap_celltype = FeaturePlot(seurat_comb_singlet, features = genemarker, reduction = "umap", ncol = ceiling(length(genemarker)/3), label = T)
    pdf_i = plot_grid(res_i, umap_celltype)
    title = ggdraw() + draw_label(celltype, fontface = 'bold', size=26)
    pdf_i = plot_grid(title, pdf_i, nrow = 2, rel_heights = c(0.1,2))
    pdf(file = paste0(workdir, "Results/Figures/6_cluster_anno/2_cluster_mapping/res.",res,"/",celltype,".pdf"), width=width, height=18)
    print(pdf_i)
    dev.off()
    options(bitmapType='cairo')
    tiff(file = paste0(workdir, "Results/Figures/6_cluster_anno/2_cluster_mapping/res.",res,"/",celltype,".tiff"), res=300, width=width, height=18, compression="lzw", units="in")
    print(pdf_i)
    dev.off()
    #pdfs = c(pdfs, list(pdf_i))
    #print(length(pdf_i))
  }
  #cat(length(pdfs))
  #pdf(file = paste0(workdir, "Results/Figures/6_cluster_anno/2_cluster_mapping/res.",res,"/all.pdf"), width=width, height=18)
  #for (pdf_i2 in pdfs) {
  #    print(pdf_i2)
  #}
  #dev.off()
}
#clus <- makeCluster(5)
#clusterExport(clus,c('workdir', 'seurat_comb_singlet', 'SigGeneral_all_sort', 'res4', 'res5', 'res8', 'res10', 'res12'),envir = environment())
#clusterEvalQ(clus, library('Seurat'))
#clusterEvalQ(clus, library('cowplot'))
#parLapply(clus, c(0.4, 0.5, 0.8, 1, 1.2) ,fun = step2)
#stopCluster(clus)
for (res in c(1, 0.8, 0.4, 0.5, 1.2)) {
  step2(res)
}

#*3.heatmap校验
step3 <- function(res){
  if (res==0.4) {
    Idents(seurat_comb_singlet) = seurat_comb_singlet$RNA_snn_res.0.4
  } else if (res==0.5) {
    Idents(seurat_comb_singlet) = seurat_comb_singlet$RNA_snn_res.0.5
  } else if (res==0.8) {
    Idents(seurat_comb_singlet) = seurat_comb_singlet$RNA_snn_res.0.8
  } else if (res==1) {
    Idents(seurat_comb_singlet) = seurat_comb_singlet$RNA_snn_res.1
  } else if (res==1.2) {
    Idents(seurat_comb_singlet) = seurat_comb_singlet$RNA_snn_res.1.2
  }
  seurat_comb_singlet_markers <- FindAllMarkers(seurat_comb_singlet, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
  save(seurat_comb_singlet_markers, file = paste0(workdir, "saveData/seurat_merge_markers_",res,".RData"))
  seurat_comb_singlet_markers %>% group_by(cluster) %>% top_n(n=10, wt = avg_log2FC) -> top10
  
  height <- 0.8 + nrow(top10)*0.11
  pdf(file = paste0(workdir, "Results/Figures/6_cluster_anno/3_heatmap_check/res.",res,".pdf"), width = 10, height = height)
  print(DoHeatmap(seurat_comb_singlet, features = top10$gene) + NoLegend())
  dev.off()
  options(bitmapType='cairo')
  tiff(file = paste0(workdir, "Results/Figures/6_cluster_anno/3_heatmap_check/res.",res,".tiff"), res=300, width=10, height=height, compression="lzw", units="in")
  print(DoHeatmap(seurat_comb_singlet, features = top10$gene) + NoLegend())
  dev.off()
}
#clusterEvalQ(clus, library('ggplot2'))
#clusterEvalQ(clus, library('dplyr'))
#parLapply(clus, c(0.4, 0.5, 0.8, 1, 1.2) ,fun = step3)
#stopCluster(clus)
for (res in c(1, 0.8, 0.4, 0.5, 1.2)) {
  step3(res)
}

step3(0.4)

########singleR annotation
Idents(seurat_comb_singlet) = seurat_comb_singlet$RNA_snn_res.1
seurat_singlet_SingleR <- GetAssayData(seurat_comb_singlet, slot="data")
clusters=seurat_comb_singlet@meta.data$seurat_clusters
mouseImmu <- celldex::ImmGenData()
pred.mouseImmu <- SingleR(test = seurat_singlet_SingleR, ref = mouseImmu, labels = mouseImmu$label.main,
                          method = "cluster", clusters = clusters, 
                          assay.type.test = "logcounts", assay.type.ref = "logcounts")
mouseRNA <- celldex::MouseRNAseqData()
pred.mouseRNA <- SingleR(test = seurat_singlet_SingleR, ref = mouseRNA, labels = mouseRNA$label.fine ,
                         method = "cluster", clusters = clusters, 
                         assay.type.test = "logcounts", assay.type.ref = "logcounts")

cellType=data.frame(ClusterID=levels(seurat_comb_singlet@meta.data$seurat_clusters),
                    mouseImmu=pred.mouseImmu$labels,
                    mouseRNA=pred.mouseRNA$labels)

pred.mouseImmu.ids <- cellType$mouseImmu #把预测细胞类型与cluser进行匹配
names(pred.mouseImmu.ids) <-  cellType$ClusterID
seurat_comb_singlet_1 <- RenameIdents(seurat_comb_singlet, pred.mouseImmu.ids)#进行重新命名
pdf(file = paste0(workdir,"Results/Figures/6_cluster_anno/4_anno_result/umap_anno_singleR_mouseImmu.pdf"),width = 7, height = 6)
print(DimPlot(seurat_comb_singlet_1, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend())
dev.off()
rm(seurat_comb_singlet_1)
gc()
####
pred.mouseRNA.ids <- cellType$mouseRNA #把预测细胞类型与cluser进行匹配
names(pred.mouseRNA.ids) <-  cellType$ClusterID
seurat_comb_singlet_2 <- RenameIdents(seurat_comb_singlet, pred.mouseRNA.ids)#进行重新命名
pdf(file = paste0(workdir,"Results/Figures/6_cluster_anno/4_anno_result/umap_anno_singleR_mouseRNA.pdf"),width = 7, height = 6)
print(DimPlot(seurat_comb_singlet_2, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend())
dev.off()
rm(seurat_comb_singlet_2)
gc()


Idents(seurat_comb_singlet) = seurat_comb_singlet$RNA_snn_res.1
anno_mannual <- c(13:'CD44+MMP19+')
names(anno_mannual) <- c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15)
seurat_comb_singlet_3 <- RenameIdents(seurat_comb_singlet, anno_mannual)#进行重新命名
pdf(file = paste0(workdir,"Results/Figures/6_cluster_anno/4_anno_result/umap_anno_mannual.pdf"),width = 7, height = 6)
print(DimPlot(seurat_comb_singlet_3, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend())
dev.off()
rm(seurat_comb_singlet_3)
gc()
