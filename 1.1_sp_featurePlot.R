### Seurat clustering
library(Seurat)
library(dplyr)
library(cowplot)
library(patchwork)
#library(DoubletFinder)
library(ggplot2)
library(future)
library(Cairo)
options(bitmapType='cairo')
#library(psych)
#library(qgraph)
library(igraph)
library(Matrix)
#library(SeuratWrappers)
#library(pryr)
library(MuDataSeurat)

source("/jdfsbjcas1/ST_BJ/P21H28400N0232/wangdi/tumors_public/codes/scRNA_primary.R")
source("/jdfsbjcas1/ST_BJ/P21H28400N0232/wangdi/tumors_public/codes/colors.R")
source("/jdfsbjcas1/ST_BJ/P21H28400N0232/wangdi/tumors_public/codes/signatureGenes.R")
source("/jdfsbjcas1/ST_BJ/P21H28400N0232/wangdi/tumors_public/codes/0_initNEWanalysis.R")
source("/jdfsbjcas1/ST_BJ/P21H28400N0232/wangdi/tumors_public/codes/DEG_wilcox.R")

args <- commandArgs(T)
# topdir <- "/jdfsbjcas1/ST_BJ/P21H28400N0232/xieguixiang/BLCA/HBCP18/binbase/figures/"
workdir <- '/jdfsbjcas1/ST_BJ/P21H28400N0232/xieguixiang/BLCA/HBCP18/scfoundation_in/'

if (!dir.exists(workdir)) {
  dir.create(workdir)
}

seurat_spatialObj <- ReadH5AD('/jdfsbjcas1/ST_BJ/P21H28400N0232/xieguixiang/BLCA/HBCP18/scfoundation_in/adata_Epi_scFoundation.h5ad')
setwd(workdir)
# ### mannual annotation ###
# SigGeneral_all_sort <- list(
#   Basal_tumor_cell = c("KRT5","KRT14","KRT6A"),
#   muscle_cell = c("MYH11","NR2F2","CRYAB","LMOD1","TPPP3"),
#   lymphoid = c("PTPRC","CD247"),
#   myeloid = c("CSF2RA","CSTB"),
#   macrophage = c("ADGRE1","CD163", "SLC11A1", "APOC1", "CD86", "CSF1R", "SLCO2B1", "CD68", "F13A1", "CD14", "AIF1", "CD80","FCER1G", "FCGR3A", "TYROBP","LYZ","MS4A7"),
#   fibroblast = c("CALD1","LUM","COL1A1","COL1A2","ACTA2", "SULF1", "CTGF", "TAGLN","TPM1", "GINS1" ,"THY1", "RBP1", "COL1A2", "COL1A1", "C1R", "IGFBP7", "SFRP2", "MGP","C1S", "DCN", "CXCL14", "COL3A1","FAP","COL6A1","PDPN"),
#   Tcell = c("CCL5","CD52","KLRB1","CD2", "CD3D", "CD3E", "CD3G", "CD8A", "CD8B", "GZMK", "CD4","TNFRSF4","IL7R"),
#   Bcell = c("KIF4A","FANCI","CHAF1A", "GTSE1", "ASPM", "SPC25","NCAPG2","POLA2","NCAPD3", "CD19", "CD21", "MS4A1", "CD79A", "CD79B", "BLNK","MZB1"),
#   DCs = c("HLA-DRB1","ITGAX","CD83", "HLA-DMA","HLA-DQB1", "HLA-DPA1","HLA-DPB1","AIF1","LST1","FTL","HLA-C","CD209"),
#   Epithelial = c("HSPA6","S100A2","KRT17","KRT5","SLPI","CXCL17","C1orf56","S100A9","S100A8","IFI44L","EPCAM", "KRT3", "KRT14", "MUC1", "TP63","CDH1"),
#   Myeloid = c("SPP1","APOE","C1QB"),
#   Urothelial_cell =c("KRT13","UPK2","UPK1B","UPK3A"),
#   Epithe_KRT13_17 = c("KRT13","KRT17"),
#   Epithe_CDH12 = c("CDH18","CDH12"),
#   Epithe_cycling = c("TUBA1B","EEF2"),
#   Epithe_Lum = c("HIF1A","CEBPD","BTG1","KRT8","CD9","AQP3","MUC1","KRT18","SLPI","AGR2","LCN2","CLDN4","ANXA1","CD74"),
#   Epithe_Mam = c("PRLR","CLDN4","CSN3","CSN1S1","KRT19","KRT7","KRT8","KRT18"),
#   PVL = c("MCAM", "CD146", "ACTA2", "PDGFRB", "LYVE1","RGS5","COL4A1","CALD1","SPARCL1","SPARC","NID1"),
#   Pericyte = c("ACTA2","PDGFRB","MCAM","HIGD1B","ANGPT2","VIM","MFGE8","MYO1B","NOTCH3","COX4I2"),
#   endothelial = c("FLT1","PLVAP","SPARCL1","GNG11","PECAM1", "CD31", "CD34", "HSPG2", "LDB2", "GPR116", "PTPRB", "VWF", "DOCK9", "CDH5", "SELE","VCAM1","ENG"),
#   endocrine = c("CHGB", "CHGA", "TTR", "SCG5", "SLC30A8", "GCG", "CLU","CPE","SCG3","CRYBA2","TM4SF4","SCGN"),
#   monocyte = c("LYZ","MS4A7","CD14"),
#   tuftCells = c("AZGP1", "PLCG2", "HPGDS", "AVIL", "PAEP", "SH2D6", "BMX", "LRMP"),
#   neutrophils = c("A1BG","ALOX5","ASAH1","CD33","CD44","CD63","CTSG","DOCK2","HSPA1B","HSP90AA1"),
#   mastCells = c("RHOH","BTK","FER","GATA2","IL4R","KIT","LCP2","ENPP3","RAC2","LAT2","CD84","LAT","ADGRE2","UNC13D","NDEL1", 
#                 "CPA3","TPSB2","TPSAB1","MS4A2","SLC18A2","IL1RL1","PTGS1","HPGDS"),
#   NKcell = c("GNLY", "FCER1G", "KLRB1", "KLRC1", "AREG", "XCL1"),
#   NKT = c("NKG7", "GNLY", "GZMA", "GZMB", "FCGR3A", "KLRB1"),
#   ducat1Cell = c("AMBP", "FXYD2"),
#   plasmocyte = c("IGLC2","IGLC3","IGHG3","IGHG1","IGHG4","JCHAIN","IGLC1"),
#   DuctalCell = c("CAPS","TPPP3","MIA","RSPH1","PIFO","LCN2","GDF15","AGR3","CETN2"),
#   Cyclegene = c("MKI67","CDC20", "CENPF","PTTG1","TOP2A","CCNB1","PCNA"),
#   astrocyte = c("GFAP", "BMPR1B", "CD44", "SLC1A2", "AQP4", "S100B", "GJB7", "ALDH1L1", "ALDOC", "MLC1"),
#   oligodendrocyte = c("MBP", "SOX10", "MOG", "CA2", "CNP", "RTN4", "PLP1", "PLP2", "OPALIN", "OMG","OLIG1", "TNR", "ALCAM", "PLLP"),
#   MSC = c("CD44", "ITGA1","NT5E","THY1"),
#   mycaf = c("RGS5","ACTA2","VIM","CCN2", "COL1A1","COL5A1","COL6A1","TNC","TGFB1","THY1","TAGLN","COL12A1","PDGFRB"), #TGF尾/SMAD2/3 
#   icaf = c("PDGFRA","IL1A","IL1B","IL6","IL11","LIF","CLEC3B","COL14A1","GSN","LY6C1","CXCL12","CXCL14"), ##IL-1/JAK-STAT3
#   apcaf = c("SLPI","SAA3","CD74","H2-Ab1","NKAIN4", "IRF5"),
#   psc = c("DES", "GFAP", "CHRNA1"),
#   Luminal = c('KRT20','PPARG','FOXA1','GATA3','SNX31','UPK1A','UPK2','FGFR3'),
#   ECM_smooth_muscle = c('PGM5','DES','C7','SFRP4','COMP','SGCD'),
#   EMT_Claudin = c('ZEB1','ZEB2','SNAI1','TWIST1','CDH2','CLDN3','CLDN4','CLDN7'),
#   Basal = c('CD44','KRT6A','KRT5','KRT14','COL17A1'),
#   Squamous = c('DSC3','GSDMC','TGM1','PI3','TP63'),
#   Immune = c('CD274','PDCD1LG2','IDO1','CXCL11','L1CAM','SAA1'),
#   Neuronal_differentiation = c('MSI1','PLEKHG4B','GNG4','PEG10','RND2','APLP1','SOX2','TUBB2B'),
#   CIS_down = c('CRTAC1','CTSE','PADI3'),
#   CIS_up = c('MSN','NR3C1')
# )

# load('/jdfsbjcas1/ST_BJ/P21H28400N0232/xieguixiang/BLCA/SC/R28BC/fibroblast/intersectHVG_15/saveData/seurat_comb_markers.RData')
# seurat_comb_singlet_markers %>%
#      group_by(cluster) %>%
#      top_n(n = 30, wt = avg_log2FC) -> top30

# SigGeneral_all_sort <- list(
#   type8 = top30$gene[top30$cluster=='8'],
#   type11 = top30$gene[top30$cluster=='11'],
#   type4 = top30$gene[top30$cluster=='4'],
#   type2 = top30$gene[top30$cluster=='2'],
#   type15 = top30$gene[top30$cluster=='15']
# )

# library('parallel')
# options(bitmapType='cairo')

# #*2.cluster
# step2 <- function(res){
#   #res_i = get(paste0('Spatial_snn_resz', res))
#   if (res==0.4) {
#     Idents(seurat_spatialObj) <- seurat_spatialObj$Spatial_snn_res.0.4
#   } else if (res==0.5) {
#     Idents(seurat_spatialObj) <- seurat_spatialObj$Spatial_snn_res.0.5
#   } else if (res==0.8) {
#     Idents(seurat_spatialObj) <- seurat_spatialObj$Spatial_snn_res.0.8
#   } else if (res==1) {
#     Idents(seurat_spatialObj) <- seurat_spatialObj$Spatial_snn_res.1
#   } else if (res==1.2) {
#     Idents(seurat_spatialObj) <- seurat_spatialObj$Spatial_snn_res.1.2
#   }
#   #pdfs = list()
#   for (celltype in names(SigGeneral_all_sort)) {
#     print(c(res,celltype))
#     genemarker = SigGeneral_all_sort[[celltype]]
#     genemarker = genemarker[genemarker %in% rownames(seurat_spatialObj)]
#     if(length(genemarker)<1){next}
#     width <- 18+ceiling(length(genemarker)/3) * 6
#     umap_celltype = FeaturePlot(seurat_spatialObj, features = genemarker, reduction = "umap", ncol = ceiling(length(genemarker)/3), label = T)
#     pdf_i = plot_grid(umap_celltype)
#     title = ggdraw() + draw_label(celltype, fontface = 'bold', size=26)
#     pdf_i = plot_grid(title, pdf_i, nrow = 1, rel_heights = c(0.1,2))
    
#     pdf(file = paste0(workdir,"/",celltype,".pdf"), width=width, height=18)
#     print(pdf_i)
#     dev.off()
#     options(bitmapType='cairo')
#     #pdfs = c(pdfs, list(pdf_i))
#     #print(length(pdf_i))
#   }
#   #cat(length(pdfs))
#   #pdf(file = paste0(workdir, "Results/Figures/6_cluster_anno/2_cluster_mapping/res.",res,"/all.pdf"), width=width, height=18)
#   #for (pdf_i2 in pdfs) {
#   #    print(pdf_i2)
#   #}
#   #dev.off()
# }
#step2(0.6)

#Idents(seurat_spatialObj) <- seurat_spatialObj$Spatial_snn_res.0.4
#anno_mannual <- c("Myofibroblasts","Epithelials","Myofibroblasts","Fibroblasts","Epithelials","Epithelials","Endothelials","Plasmocytes","Epithelials","Endothelials","Fibroblasts","Fibroblasts","Epithelials")
#names(anno_mannual) <- c(0,1,2,3,4,5,6,7,8,9,10,11,12)
#seurat_spatialObj_2 <- RenameIdents(seurat_spatialObj, anno_mannual)#杩涜閲嶆柊鍛藉悕
#pdf(file = paste0(workdir,"Spatial_res_anno.pdf"),width = 7, height = 6)
#SpatialDimPlot(seurat_spatialObj_2, label = TRUE, label.size = 3, stroke = 0)
#LinkedDimPlot(seurat_spatialObj_2)
#SpatialFeaturePlot(seurat_spatialObj_2, features = c("SKA3","CCNB1","FOXM1","GSN","LAMC2"))
#SpatialFeaturePlot(seurat_spatialObj_2, features = c("ACTA2", "KRT19"))
#SpatialFeaturePlot(seurat_spatialObj_2, features = c("ACTA2", "TAGLN"))
#print(DimPlot(seurat_comb_singlet_2, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend())
#dev.off()
#rm(seurat_comb_singlet_2)
#gc()

#pdf(file = paste0(workdir,"cancer_marker.pdf"),width = 7, height = 6)
#SpatialFeaturePlot(seurat_spatialObj, features = c("SKA3")) +
#    theme(legend.position = "right")
#dev.off()


# for (celltype in names(SigGeneral_all_sort)) {
#   print(celltype)
#   genemarker = SigGeneral_all_sort[[celltype]]
#   genemarker = genemarker[genemarker %in% rownames(seurat_spatialObj)]
#   if(length(genemarker)<1){next}
#   width <- 18+ceiling(length(genemarker)/3) * 4
#   spatial_celltype = SpatialFeaturePlot(seurat_spatialObj, features = genemarker, crop = TRUE, stroke=0)
#   #spatial_celltype = SpatialFeaturePlot(seurat_spatialObj, features = genemarker, crop = TRUE, alpha = c(0.1, 1), pt.size.factor = 1.3)
#   pdf_i = plot_grid(spatial_celltype)
#   title = ggdraw() + draw_label(celltype, fontface = 'bold', size=26)
#   pdf_i = plot_grid(title, pdf_i, nrow = 2, rel_heights = c(0.1,2))
#   pdf(file = paste0(workdir,"/",celltype,".pdf"), width=120, height=100)
#   print(pdf_i)
#   dev.off()
#   options(bitmapType='cairo')
# }

Basal <- list(c("KRT5","KRT14","KRT6A","CD44","KRT6","KRT17","CDH3","MMP14"))
Luminal <- list(c('KRT20','PPARG','FOXA1','GATA3','SNX31','UPK1A','UPK2','FGFR3','ESR1','PGR','KRT8','KRT18'))
OXPHOS <- list(c("ABCB7","ACAA1","ACAA2","ACADM","ACADSB","ACADVL","ACAT1","ACO2","AFG3L2","AIFM1","ALAS1","ALDH6A1","ATP1B1","ATP5F1A","ATP5F1B","ATP5F1C","ATP5F1D","ATP5F1E","ATP5MC1","ATP5MC2","ATP5MC3","ATP5ME","ATP5MF","ATP5MG","ATP5PB","ATP5PD","ATP5PF","ATP5PO","ATP6AP1","ATP6V0B","ATP6V0C","ATP6V0E1","ATP6V1C1","ATP6V1D","ATP6V1E1","ATP6V1F","ATP6V1G1","ATP6V1H","BAX","BCKDHA","BDH2","CASP7","COX10","COX11","COX15","COX17","COX4I1","COX5A","COX5B","COX6A1","COX6B1","COX6C","COX7A2","COX7A2L","COX7B","COX7C","COX8A","CPT1A","CS","CYB5A","CYB5R3","CYC1","CYCS","DECR1","DLAT","DLD","DLST","ECH1","ECHS1","ECI1","ETFA","ETFB","ETFDH","FDX1","FH","FXN","GLUD1","GOT2","GPI","GPX4","GRPEL1","HADHA","HADHB","HCCS","HSD17B10","HSPA9","HTRA2","IDH1","IDH2","IDH3A","IDH3B","IDH3G","IMMT","ISCA1","ISCU","LDHA","LDHB","LRPPRC","MAOB","MDH1","MDH2","MFN2","MGST3","MPC1","MRPL11","MRPL15","MRPL34","MRPL35","MRPS11","MRPS12","MRPS15","MRPS22","MRPS30","MTRF1","MTRR","MTX2","NDUFA1","NDUFA2","NDUFA3","NDUFA4","NDUFA5","NDUFA6","NDUFA7","NDUFA8","NDUFA9","NDUFAB1","NDUFB1","NDUFB2","NDUFB3","NDUFB4","NDUFB5","NDUFB6","NDUFB7","NDUFB8","NDUFC1","NDUFC2","NDUFS1","NDUFS2","NDUFS3","NDUFS4","NDUFS6","NDUFS7","NDUFS8","NDUFV1","NDUFV2","NNT","NQO2","OAT","OGDH","OPA1","OXA1L","PDHA1","PDHB","PDHX","PDK4","PDP1","PHB2","PHYH","PMPCA","POLR2F","POR","PRDX3","RETSAT","RHOT1","RHOT2","SDHA","SDHB","SDHC","SDHD","SLC25A11","SLC25A12","SLC25A20","SLC25A3","SLC25A4","SLC25A5","SLC25A6","SUCLA2","SUCLG1","SUPV3L1","SURF1","TCIRG1","TIMM10","TIMM13","TIMM17A","TIMM50","TIMM8B","TIMM9","TOMM22","TOMM70","UQCR10","UQCR11","UQCRB","UQCRC1","UQCRC2","UQCRFS1","UQCRH","UQCRQ","VDAC1","VDAC2","VDAC3"))
EMT <- list(c("ABI3BP", "ACTA2", "ADAM12", "ANPEP", "APLP1", "AREG", "BASP1", "BDNF", "BGN", "BMP1", "CADM1", "CALD1", "CALU", "CAP2", "CAPG", "CCN1", "CCN2", "CD44", "CD59", "CDH11", "CDH2", "CDH6", "COL11A1", "COL12A1", "COL16A1", "COL1A1", "COL1A2", "COL3A1", "COL4A1", "COL4A2", "COL5A1", "COL5A2", "COL5A3", "COL6A2", "COL6A3", "COL7A1", "COL8A2", "COLGALT1", "COMP", "COPA", "CRLF1", "CTHRC1", "CXCL1", "CXCL12", "CXCL6", "CXCL8", "DAB2", "DCN", "DKK1", "DPYSL3", "DST", "ECM1", "ECM2", "EDIL3", "EFEMP2", "ELN", "EMP3", "ENO2", "FAP", "FAS", "FBLN1", "FBLN2", "FBLN5", "FBN1", "FBN2", "FERMT2", "FGF2", "FLNA", "FMOD", "FN1", "FOXC2", "FSTL1", "FSTL3", "FUCA1", "FZD8", "GADD45A", "GADD45B", "GAS1", "GEM", "GJA1", "GLIPR1", "GPC1", "GPX7", "GREM1", "HTRA1", "ID2", "IGFBP2", "IGFBP3", "IGFBP4", "IL15", "IL32", "IL6", "INHBA", "ITGA2", "ITGA5", "ITGAV", "ITGB1", "ITGB3", "ITGB5", "JUN", "LAMA1", "LAMA2", "LAMA3", "LAMC1", "LAMC2", "LGALS1", "LOX", "LOXL1", "LOXL2", "LRP1", "LRRC15", "LUM", "MAGEE1", "MATN2", "MATN3", "MCM7", "MEST", "MFAP5", "MGP", "MMP1", "MMP14", "MMP2", "MMP3", "MSX1", "MXRA5", "MYL9", "MYLK", "NID2", "NNMT", "NOTCH2", "NT5E", "NTM", "OXTR", "P3H1", "PCOLCE", "PCOLCE2", "PDGFRB", "PDLIM4", "PFN2", "PLAUR", "PLOD1", "PLOD2", "PLOD3", "PMEPA1", "PMP22", "POSTN", "PPIB", "PRRX1", "PRSS2", "PTHLH", "PTX3", "PVR", "QSOX1", "RGS4", "RHOB", "SAT1", "SCG2", "SDC1", "SDC4", "SERPINE1", "SERPINE2", "SERPINH1", "SFRP1", "SFRP4", "SGCB", "SGCD", "SGCG", "SLC6A8", "SLIT2", "SLIT3", "SNAI2", "SNTB1", "SPARC", "SPOCK1", "SPP1", "TAGLN", "TFPI2", "TGFB1", "TGFBI", "TGFBR3", "TGM2", "THBS1", "THBS2", "THY1", "TIMP1", "TIMP3", "TNC", "TNFAIP3", "TNFRSF11B", "TNFRSF12A", "TPM1", "TPM2", "TPM4", "VCAM1", "VCAN", "VEGFA", "VEGFC", "VIM", "WIPF1", "WNT5A"))
Glycolysis <- list(c("ABCB6", "ADORA2B", "AGL", "AGRN", "AK3", "AK4", "AKR1A1", "ALDH7A1", "ALDH9A1", "ALDOA", "ALDOB", "ALG1", "ANG", "ANGPTL4", "ANKZF1", "ARPP19", "ARTN", "AURKA", "B3GALT6", "B3GAT1", "B3GAT3", "B3GNT3", "B4GALT1", "B4GALT2", "B4GALT4", "B4GALT7", "BIK", "BPNT1", "CACNA1H", "CAPN5", "CASP6", "CD44", "CDK1", "CENPA", "CHPF", "CHPF2", "CHST1", "CHST12", "CHST2", "CHST4", "CHST6", "CITED2", "CLDN3", "CLDN9", "CLN6", "COG2", "COL5A1", "COPB2", "CTH", "CXCR4", "CYB5A", "DCN", "DDIT4", "DEPDC1", "DLD", "DPYSL4", "DSC2", "ECD", "EFNA3", "EGFR", "EGLN3", "ELF3", "ENO1", "ENO2", "ERO1A", "EXT1", "EXT2", "FAM162A", "FBP2", "FKBP4", "FUT8", "G6PD", "GAL3ST1", "GALE", "GALK1", "GALK2", "GAPDHS", "GCLC", "GFPT1", "GFUS", "GLCE", "GLRX", "GMPPA", "GMPPB", "GNE", "GNPDA1", "GOT1", "GOT2", "GPC1", "GPC3", "GPC4", "GPR87", "GUSB", "GYS1", "GYS2", "HAX1", "HDLBP", "HK2", "HMMR", "HOMER1", "HS2ST1", "HS6ST2", "HSPA5", "IDH1", "IDUA", "IER3", "IGFBP3", "IL13RA1", "IRS2", "ISG20", "KDELR3", "KIF20A", "KIF2A", "LCT", "LDHA", "LDHC", "LHPP", "LHX9", "MDH1", "MDH2", "ME1", "ME2", "MED24", "MERTK", "MET", "MIF", "MIOX", "MPI", "MXI1", "NANP", "NASP", "NDST3", "NDUFV3", "NOL3", "NSDHL", "NT5E", "P4HA1", "P4HA2", "PAM", "PAXIP1", "PC", "PDK3", "PFKFB1", "PFKP", "PGAM1", "PGAM2", "PGK1", "PGLS", "PGM2", "PHKA2", "PKM", "PKP2", "PLOD1", "PLOD2", "PMM2", "POLR3K", "PPFIA4", "PPIA", "PPP2CB", "PRPS1", "PSMC4", "PYGB", "PYGL", "QSOX1", "RARS1", "RBCK1", "RPE", "RRAGD", "SAP30", "SDC1", "SDC2", "SDC3", "SDHC", "SLC16A3", "SLC25A10", "SLC25A13", "SLC35A3", "SLC37A4", "SOD1", "SOX9", "SPAG4", "SRD5A3", "STC1", "STC2", "STMN1", "TALDO1", "TFF3", "TGFA", "TGFBI", "TKTL1", "TPBG", "TPI1", "TPST1", "TXN", "UGP2", "VCAN", "VEGFA", "VLDLR", "XYLT2", "ZNF292"))
Hypoxia <- list(c("ACKR3","ADM","ADORA2B","AK4","AKAP12","ALDOA","ALDOB","ALDOC","AMPD3","ANGPTL4","ANKZF1","ANXA2","ATF3","ATP7A","B3GALT6","B4GALNT2","BCAN","BCL2","BGN","BHLHE40","BNIP3L","BRS3","BTG1","CA12","CASP6","CAV1","CAVIN1","CAVIN3","CCN1","CCN2","CCN5","CCNG2","CDKN1A","CDKN1B","CDKN1C","CHST2","CHST3","CITED2","COL5A1","CP","CSRP2","CXCR4","DCN","DDIT3","DDIT4","DPYSL4","DTNA","DUSP1","EDN2","EFNA1","EFNA3","EGFR","ENO1","ENO2","ENO3","ERO1A","ERRFI1","ETS1","EXT1","F3","FAM162A","FBP1","FOS","FOSL2","FOXO3","GAA","GALK1","GAPDH","GAPDHS","GBE1","GCK","GCNT2","GLRX","GPC1","GPC3","GPC4","GPI","GRHPR","GYS1","HAS1","HDLBP","HEXA","HK1","HK2","HMOX1","HOXB9","HS3ST1","HSPA5","IDS","IER3","IGFBP1","IGFBP3","IL6","ILVBL","INHA","IRS2","ISG20","JMJD6","JUN","KDELR3","KDM3A","KIF5A","KLF6","KLF7","KLHL24","LALBA","LARGE1","LDHA","LDHC","LOX","LXN","MAFF","MAP3K1","MIF","MT1E","MT2A","MXI1","MYH9","NAGK","NCAN","NDRG1","NDST1","NDST2","NEDD4L","NFIL3","NOCT","NR3C1","P4HA1","P4HA2","PAM","PCK1","PDGFB","PDK1","PDK3","PFKFB3","PFKL","PFKP","PGAM2","PGF","PGK1","PGM1","PGM2","PHKG1","PIM1","PKLR","PKP1","PLAC8","PLAUR","PLIN2","PNRC1","PPARGC1A","PPFIA4","PPP1R15A","PPP1R3C","PRDX5","PRKCA","PYGM","RBPJ","RORA","RRAGD","S100A4","SAP30","SCARB1","SDC2","SDC3","SDC4","SELENBP1","SERPINE1","SIAH2","SLC25A1","SLC2A1","SLC2A3","SLC2A5","SLC37A4","SLC6A6","SRPX","STBD1","STC1","STC2","SULT2B1","TES","TGFB3","TGFBI","TGM2","TIPARP","TKTL1","TMEM45A","TNFAIP3","TPBG","TPD52","TPI1","TPST2","UGP2","VEGFA","VHL","VLDLR","WSB1","XPNPEP1","ZFP36","ZNF292"))



seurat_spatialObj <- AddModuleScore(seurat_spatialObj,
                          features = Basal,
                          ctrl = 100,
                          name = "Basal_score")
seurat_spatialObj <- AddModuleScore(seurat_spatialObj,
                          features = Luminal,
                          ctrl = 100,
                          name = "Luminal_score")
seurat_spatialObj <- AddModuleScore(seurat_spatialObj,
                          features = OXPHOS,
                          ctrl = 100,
                          name = "OXPHOS_score")
seurat_spatialObj <- AddModuleScore(seurat_spatialObj,
                          features = EMT,
                          ctrl = 100,
                          name = "EMT_score")
seurat_spatialObj <- AddModuleScore(seurat_spatialObj,
                          features = Glycolysis,
                          ctrl = 100,
                          name = "Glycolysis_score")
seurat_spatialObj <- AddModuleScore(seurat_spatialObj,
                          features = Hypoxia,
                          ctrl = 100,
                          name = "Hypoxia_score")
score_names <- c("Basal_score1", "Luminal_score1", "OXPHOS_score1", "EMT_score1", "Glycolysis_score1", "Hypoxia_score1")


# 循环处理每个得分列
for (score_name in score_names) {
  # 计算当前得分列的 10% 最大值阈值
  threshold <- quantile(seurat_spatialObj@meta.data[[score_name]], 0.9)
  # 突出最大值的 10%
  seurat_spatialObj@meta.data[[score_name]] <- ifelse(seurat_spatialObj@meta.data[[score_name]] >= threshold, seurat_spatialObj@meta.data[[score_name]] * 2 , seurat_spatialObj@meta.data[[score_name]])
}
head(seurat_spatialObj@meta.data)

# 重命名列名
colnames(seurat_spatialObj@meta.data)[which(colnames(seurat_spatialObj@meta.data) == "Basal_score1")] <- "Basal"
colnames(seurat_spatialObj@meta.data)[which(colnames(seurat_spatialObj@meta.data) == "Luminal_score1")] <- "Luminal"
colnames(seurat_spatialObj@meta.data)[which(colnames(seurat_spatialObj@meta.data) == "OXPHOS_score1")] <- "OXPHOS"
colnames(seurat_spatialObj@meta.data)[which(colnames(seurat_spatialObj@meta.data) == "EMT_score1")] <- "EMT"
colnames(seurat_spatialObj@meta.data)[which(colnames(seurat_spatialObj@meta.data) == "Glycolysis_score1")] <- "Glycolysis"
colnames(seurat_spatialObj@meta.data)[which(colnames(seurat_spatialObj@meta.data) == "Hypoxia_score1")] <- "Hypoxia"


pdf(file = paste0(workdir,"Basal_score.pdf"))
SpatialFeaturePlot(seurat_spatialObj, features = 'Basal', crop = TRUE, stroke=0)
dev.off()
pdf(file = paste0(workdir,"Luminal_score.pdf"))
SpatialFeaturePlot(seurat_spatialObj, features = 'Luminal', crop = TRUE, stroke=0)
dev.off()
pdf(file = paste0(workdir,"OXPHOS_score.pdf"))
SpatialFeaturePlot(seurat_spatialObj, features = 'OXPHOS', crop = TRUE, stroke=0)
dev.off()
pdf(file = paste0(workdir,"EMT_score.pdf"))
SpatialFeaturePlot(seurat_spatialObj, features = 'EMT', crop = TRUE, stroke=0)
dev.off()
pdf(file = paste0(workdir,"Glycolysis_score.pdf"))
SpatialFeaturePlot(seurat_spatialObj, features = 'Glycolysis', crop = TRUE, stroke=0)
dev.off()
pdf(file = paste0(workdir,"Hypoxia_score.pdf"))
SpatialFeaturePlot(seurat_spatialObj, features = 'Hypoxia', crop = TRUE, stroke=0)
dev.off()


MIBC_features <- c('COX6B1','ISG15','BST2','MAL2','CRABP2','CRIP1','IFI44','KRT5','TYMP','ST14','DHRS3','PMEPA1','LY6K','MYEOV')
tiff(file = paste0('MIBC_genes.tif'), width = 10, height = 10, units = 'in', res = 300)
SpatialFeaturePlot(seurat_spatialObj, features = MIBC_features, crop = TRUE, stroke=0)
dev.off()

metastatic_features <- c('GDPD2','TNNI2','ITM2C','PKP1','FAM3B','CYP1A1','CYP1B1','KRT16','FABP4','PYGL','DMKN','KRT17','TRIM29','TYMP','SLC7A5','BHLHE40','HAS3')
tiff(file = paste0('metastatic_genes.tif'), width = 10, height = 10, units = 'in', res = 300)
SpatialFeaturePlot(seurat_spatialObj, features = metastatic_features, crop = TRUE, stroke=0)
dev.off()

NMIBC_features <- c('ID1','ID4','ELF3','AQP3','GDF15','CCND1','FXYD5','NFKBIA','SNCG','MIEN1','CYP24A1','CRISP3','IFITM2','TOB1','DUSP2','SPINK1')
tiff(file = paste0('NMIBC_genes.tif'), width = 10, height = 10, units = 'in', res = 300)
SpatialFeaturePlot(seurat_spatialObj, features = NMIBC_features, crop = TRUE, stroke=0)
dev.off()

test <- c('PTPRR'
)
tiff(file = paste0('test.tif'), width = 10, height = 10, units = 'in', res = 300)
SpatialFeaturePlot(seurat_spatialObj, features = test, crop = TRUE, stroke=0)
dev.off()