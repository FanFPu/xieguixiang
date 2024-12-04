library('Seurat')
library('infercnv')
library('dendextend')
library('phylogram')
library('miscTools')
library('ggthemes')
library('RColorBrewer')
library('umap')
library('ggplot2')
library('car')
library('limma')
library('tibble')
library('dplyr')
library('stringr')
library('MuDataSeurat')
options(bitmapType = "cairo")

rm(list = ls())
gc()


# args = commandArgs(T)
sample = "P18"  # args[1]

#工作路径
topdir = '/jdfsbjcas1/ST_BJ/P21H28400N0232/xieguixiang/BLCA/Region_analysis/tumor'
indir = paste0(topdir, '/infercnv/', sample)
workdir = paste0(topdir, '/infercnv/', sample, '/analysis/')
# if(!dir.exists(indir)){dir.create(indir, recursive = TRUE)}
if(!dir.exists(workdir)){dir.create(workdir, recursive = TRUE)}
setwd(workdir)

# 画CNV聚类图
k1 = as.numeric(8)
tree = read.dendrogram(paste0(indir, '/infercnv.observations_dendrogram.txt'))
labels_hclust2 = cutree(tree, k=k1)   
labels_hclust2 = data.frame(row.names = names(labels_hclust2), 'subclone'=paste0('Subclo_', labels_hclust2))  #行名细胞，列名subclone

## 取肿瘤细胞ID
load(paste0(indir, "/infercnv_obj2.RData"))  #infercnv_obj
infercnv_obj2 <-  infercnv_obj
obs_annotations_names <- names(infercnv_obj@observation_grouped_cell_indices)
labels_hclust2 <- data.frame()
for (i in seq_along(obs_annotations_names)) {
    order_indices <- c(unlist(infercnv_obj@observation_grouped_cell_indices[obs_annotations_names[i]]))
    # labels <- infercnv_obj@tumor_subclusters$hc[['all_observations']]$labels
    labels <- colnames(infercnv_obj@expr.data)
    ordered_labels <- labels[order_indices]
    labels_hclust <-  data.frame(row.names = ordered_labels, 'subclone' = c(rep(obs_annotations_names[i], length(ordered_labels))))
    labels_hclust2 <- rbind(labels_hclust2, labels_hclust)
}

#将对照细胞提取出来
ref_ids = setdiff(colnames(infercnv_obj2@expr.data), rownames(labels_hclust2))   # setdiff(A, B) 返回在集合 A 中出现但不在集合 B 中出现的元素。换句话说，它返回属于 A 但不属于 B 的元素
#将其细胞名进行合并，不漏掉一个细胞
infercnv_obj2@expr.data = infercnv_obj2@expr.data[,union(rownames(labels_hclust2), ref_ids)]   #union()计算两个或多个向量的并集,且不包含重复的元素

name = list() #建立空列表
index = list()
index_length = list()
hc_list = list()


for (i in unique(labels_hclust2[,1])) {
    name[[i]] <- rownames(labels_hclust2)[labels_hclust2$subclone == i]  #将为subclone【K】的细胞名取出来
    index[[i]] <- which(colnames(infercnv_obj2@expr.data) %in% name[[i]])  #which建立满足括号里的索引
    index_length[[i]] <- seq(from=1,to=length(name[[i]]))   #建立等差数列，总共是subclone【K】的数量
    hc_list[[i]] = list(order=seq(from=1,to=length(name[[i]])),labels=name[[i]])  # 创建两个列表，一个是index_length的列表，一个是将每个亚型的细胞贴上标签，内容是细胞名，
}

infercnv_obj2@observation_grouped_cell_indices <- index_length
infercnv_obj2@tumor_subclusters$hc <- hc_list
infercnv_obj2@reference_grouped_cell_indices = list('ref' = which(colnames(infercnv_obj2@expr.data) %in% ref_ids))  #对照的细胞名，建立索引，然后单独拿出来,内容是索引号
length(table(infercnv_obj2@reference_grouped_cell_indices)) #对照细胞类型的细胞数量
saveRDS(infercnv_obj2, file=paste0(workdir, '/infercnv_obj2.rds'))

source('/jdfsbjcas1/ST_BJ/P21H28400N0232/fanfengpu/My_code/Single_cell/SC_1.4.10_infercnv_heatmap_plot.R')  #里面有plot_cnv函数
plot_cnv(infercnv_obj2,
        out_dir=workdir,
        title="inferCNV",
        obs_title="Observations (Cells)",
        ref_title="References (Cells)",
        cluster_by_groups=T,
        cluster_references=FALSE,
        plot_chr_scale=FALSE,
        chr_lengths=NULL,
        k_obs_groups = 1,
        contig_cex=1,
        x.center=mean(infercnv_obj2@expr.data),
        x.range="auto", #NA,
        hclust_method='ward.D',
        custom_color_pal=NULL,
        color_safe_pal=FALSE,
        output_filename="infercnv",
        output_format="png" , #pdf, png, NA
        png_res=300,
        dynamic_resize=0,
        ref_contig = NULL,
        write_expr_matrix=F,
        useRaster=TRUE)

#*根据cnv树状图区分亚型
#load(paste0(workdir, '/infercnv_obj2.RData'))
#sample = "HBCP10"
tree = read.dendrogram(paste0(indir, '/infercnv.observations_dendrogram.txt'))  #从文件或其他数据源中读取聚类树，其实并不直接读文件，而是对特定R对象创建dendrogram（聚类数）对象
labels = cutree(tree, k=k1, order_clusters_as_data = F)  #k具有所需组数的整数标量或向量
### 提取正常细胞的ID
cell_names <- names(labels[labels %in% c(6)])

#*画umap图--根据ID贴上亚型标签
x <- load('/jdfsbjcas1/ST_BJ/P21H28400N0232/xieguixiang/BLCA/SC/Mouse/yangzongzheng/infercnv/saveData/seurat_merge_annot.RData')
seurat_merge <- get(x)
#将其所有上皮命名为Unselected，然后将单个样本根据cnv种类的不同将单个样本中分成正常和肿瘤等等，其他样本还是Unselected
seurat_merge@meta.data <- data.frame(seurat_merge@meta.data, cellTypes = c("Unselected"))
seurat_merge@meta.data[rownames(seurat_merge@meta.data[rownames(seurat_merge@meta.data) %in% cell_names,]), "cellTypes"] = "subclo6"


#画所选细胞在整个上皮的表达图
pdf(file = paste0(workdir,"/subclo6.pdf"),width = 7, height = 6)
DimPlot(seurat_merge, reduction = "umap", pt.size = 0.05, group.by = "cellTypes", label = FALSE, cells.highlight = list("subclo6"=cell_names), cols.highlight = list("subclo6"="red"))
# print(DimPlot(seurat_merge, reduction = "umap", pt.size = 0.05, group.by = "cellTypes", label = FALSE, cols= c('subclo8&9'= "red",'Unselected'= "grey")))
# print(DimPlot(seurat_merge4, reduction = "umap", pt.size = 0.2, group.by = "cellTypes", label = TRUE,cells.highlight = list("HBCP10-Normal"=rownames(seurat_Epithe@meta.data[seurat_Epithe@meta.data$cellSubtype_infercnv %in% c("CNV-2") , ]) )), cols.highlight = list("HBCP10-Normal" = "red")) #cols= c("blue","grey","red")    cols= c('HBCP4B-tumor'= "blue",'Unselected'= "grey")
# print(DimPlot(seurat_merge4, reduction = "umap", pt.size = 0.2, group.by = "cellTypes", label = TRUE,cols= c("blue","red","grey"))) #cols= c("blue","grey","red")    cols= c('HBCP4B-tumor'= "blue",'Unselected'= "grey")
dev.off()



# #对整个上皮细胞进行注释
# load('/jdfsbjcas1/ST_BJ/P21H28400N0232/fanfengpu/BC/R28BC_cnv_1/Results/1.4.2_analyze_inferCNV/ref_T_B_macrophage_endothelial/HBCP19/1.4.2_seurat_Epithe_scaled_HBCP19.RData')
# load('/jdfsbjcas1/ST_BJ/P21H28400N0232/fanfengpu/BC/R28BC_cnv_1/Results/1.4.2_analyze_inferCNV/ref_T_B_macrophage_endothelial/HBCP22/1.4.2_seurat_Epithe_scaled_HBCP22.RData')
# seurat_merge4 <- seurat_merge3
# seurat_merge4@meta.data <- data.frame(seurat_merge4@meta.data, cellTypes = c("Tumor"))
# length(seurat_merge4@meta.data[rownames(seurat_merge3@meta.data[,seurat_merge3@meta.data$RNA_snn_res.1.3 %in% c("27")])])
# a <- rownames(seurat_merge3@meta.data[seurat_merge3@meta.data$RNA_snn_res.1.3 == "27",])
# b <- rownames(seurat_merge4@meta.data[seurat_merge4@meta.data$cellTypes == "Normal",])
# c <- intersect(a,b)
# seurat_merge4@meta.data[rownames(seurat_Epithe@meta.data[seurat_Epithe@meta.data$cellSubtype_infercnv %in% c("CNV-2","CNV-4"),]) , "cellTypes"] = "Undefined"
# seurat_merge4@meta.data[c , "cellTypes"] = "Normal"
# seurat_merge4$cellTypes <- Idents(seurat_merge4)
# table(seurat_merge4$cellTypes)
# seurat_merge4$cellTypes <- as.factor(seurat_merge4$cellTypes)

# pdf(file = "/jdfsbjcas1/ST_BJ/P21H28400N0232/fanfengpu/BC/R28BC/Results/Futher_Legend/umap_cluster3.1.pdf",width = 7, height = 6)
# print(DimPlot(seurat_merge4, reduction = "umap",pt.size = 0.2, group.by = "cellTypes",cols = c("red","#75bbfd","grey")))#,cols = c("grey","red")
# print(DimPlot(seurat_merge4, reduction = "umap", pt.size = 0.2, group.by = "cellTypes",cells.highlight = list("Normal"=rownames(seurat_merge4@meta.data[seurat_merge4$cellTypes %in% c("Normal") , ]) )), cols.highlight = list("Normal" = "red"))
# dev.off()
# save(seurat_merge4,file = '/jdfsbjcas1/ST_BJ/P21H28400N0232/fanfengpu/BC/R28BC/seurat_merge_Epithelial_annot.RData')
# load('/jdfsbjcas1/ST_BJ/P21H28400N0232/fanfengpu/BC/R28BC/seurat_merge_Epithelial_annot.RData')

pdf(file = paste0(workdir,"/test.pdf"),width = 7, height = 6) 
print(DimPlot(seurat_merge, reduction = "umap", pt.size = 0.05, group.by = "group", label = TRUE,cells.highlight = list("24w"=rownames(seurat_merge@meta.data[seurat_merge@meta.data$group == '24w',]) )), cols.highlight = list("24w" = "red"))
dev.off()