#### inferCNV copy number for each tumor samples, including all cells or only tumor cells####
#包不能在2-4，2-5，login-2上面安装，因为连不上网,只能在softwware上面安装

###整条染色体或大片段染色体的增加或丢失(gain or deletions)。
###工作原理是:以一组"正常"细胞作为参考，分析肿瘤基因组上各个位置的基因表达量强度变化. 通过热图的形式展示每条染色体上的基因相对表达量，相对于正常细胞，肿瘤基因组总会过表达或者低表达。
library(rjags)
library(infercnv)
library(AnnoProbe)
library(tidyverse)
library(Seurat)
library(parallel)
library(Cairo)
library(limma)
options(bitmapType = "cairo")
library('MuDataSeurat')

# args = commandArgs(T)
# sample_i = args[1] #P13

#工作路径
topdir = "/jdfsbjcas1/ST_BJ/P21H28400N0232/xieguixiang/BLCA/Region_analysis/tumor/infercnv"
workdir = "/jdfsbjcas1/ST_BJ/P21H28400N0232/xieguixiang/BLCA/Region_analysis/tumor/infercnv"
if(!dir.exists(workdir)){dir.create(workdir, recursive = TRUE)}
setwd(workdir)
#写入数据
P3A_tumor <- ReadH5AD('/jdfsbjcas1/ST_BJ/P21H28400N0232/xieguixiang/BLCA/HBCP3A/Spoint/image_extract/Metastatic_tumor_VS_NMIBC/bigcell_merged.h5ad')
P3A_tumor$group <- 'P3A'
P3A_infer <- ReadH5AD('/jdfsbjcas1/ST_BJ/P21H28400N0232/xieguixiang/BLCA/Union_allsample/P3A_bigcell_CNV_reference.h5ad')
P3A_infer$group <- 'P3A'

P3B_tumor <- ReadH5AD('/jdfsbjcas1/ST_BJ/P21H28400N0232/xieguixiang/BLCA/HBCP3B/Spoint/image_extract/Metastatic_tumor_VS_NMIBC/bigcell_merged.h5ad')
P3B_tumor$group <- 'P3B'
P3B_infer <- ReadH5AD('/jdfsbjcas1/ST_BJ/P21H28400N0232/xieguixiang/BLCA/Union_allsample/P3B_bigcell_CNV_reference.h5ad')
P3B_infer$group <- 'P3B'

P3C_tumor <- ReadH5AD('/jdfsbjcas1/ST_BJ/P21H28400N0232/xieguixiang/BLCA/HBCP3C/Spoint/image_extract/Metastatic_tumor_VS_NMIBC/bigcell_merged.h5ad')
P3C_tumor$group <- 'P3C'
P3C_infer <- ReadH5AD('/jdfsbjcas1/ST_BJ/P21H28400N0232/xieguixiang/BLCA/Union_allsample/P3C_bigcell_CNV_reference.h5ad')
P3C_infer$group <- 'P3C'

P11_tumor <- ReadH5AD('/jdfsbjcas1/ST_BJ/P21H28400N0232/xieguixiang/BLCA/HBCP11/Spoint/image_extract/Metastatic_tumor_VS_NMIBC/bigcell_merged.h5ad')
P11_tumor$group <- 'P11'
P11_infer <- ReadH5AD('/jdfsbjcas1/ST_BJ/P21H28400N0232/xieguixiang/BLCA/Union_allsample/P11_bigcell_CNV_reference.h5ad')
P11_infer$group <- 'P11'

P18_tumor <- ReadH5AD('/jdfsbjcas1/ST_BJ/P21H28400N0232/xieguixiang/BLCA/HBCP18/Spoint/image_extract/Metastatic_tumor_VS_NMIBC/bigcell_merged.h5ad')
P18_tumor$group <- 'P18'
P18_infer <- ReadH5AD('/jdfsbjcas1/ST_BJ/P21H28400N0232/xieguixiang/BLCA/Union_allsample/P18_bigcell_CNV_reference.h5ad')
P18_infer$group <- 'P18'

seurat_merge <- merge(P3A_tumor, y=c(P3A_infer,P3B_infer,P3B_tumor,P3C_infer,P3C_tumor,P11_infer,P11_tumor,P18_infer,P18_tumor))
# x <- load("/jdfsbjcas1/ST_BJ/P21H28400N0232/xieguixiang/BLCA/SC/Mouse/yangzongzheng/BLCA_scRNA_QC/Results/Figures/6_cluster_anno/4_anno_result/anno_result.RData")
# seurat_merge <- get(x)
# seurat_merge$cellTypes_new <- Idents(seurat_merge)


# table(seurat_merge$cellTypes_new)
##normalize,标准化
seurat_merge <- NormalizeData(seurat_merge)
# save(seurat_merge, file = paste0(topdir, "/saveData/seurat_merge_annot.RData"))
#读入数据，全部样本的注释信息
#load(paste0(topdir, "/saveData/seurat_merge_annot.RData"))
#保存文件路径
seurat_path = paste0(topdir, "/saveData/1.4.1_seurat_merge_infercnv.RData")  #按照样本进行分开，里面具有每个样本的细胞类型

#将数据筛选，其次按照样本进行分开，按照样本进行分析
if (!file.exists(seurat_path)) {
    #load(paste0(topdir, "/saveData/1.2.1_seurat_merge_annot.RData"))
    seurat_merge_infercnv <- seurat_merge
    seurat_lists = SplitObject(seurat_merge_infercnv, split.by = 'group') #按照split.by这一列进行拆分
    save(seurat_lists, file = seurat_path)
}
# table(seurat_merge_infercnv$cellTypes_new) #查看merge后的信息是否subset正确

load(seurat_path)

for (sample in names(seurat_lists)) {
    if(!dir.exists(paste0(workdir, '/', sample))){dir.create(paste0(workdir, '/', sample), recursive = TRUE)}
}

##run inferCNV
run_infercnv = function(sample) {
  cat(sample)
  options(bitmapType = "cairo")  #
  seurat_i = seurat_lists[[sample]]
  save(seurat_i, file=paste0(workdir, '/', sample, '/seurat_Epithe.RData'))
  print(table(seurat_i$merged_cluster))
  # rm(seurat_lists)
  gc() # gc函数用来执行垃圾收集，自动处理内存管理
  mat_comb <- as.matrix(GetAssayData(seurat_i, assay="RNA",slot = "counts"))

  geneinfo <- annoGene(rownames(mat_comb), "SYMBOL", "human")
  geneinfo <- geneinfo[with(geneinfo,order(chr,start)),c(1,4:6)] #with(对象，命令)
  geneinfo <- geneinfo[!duplicated(geneinfo[,1]),]
  rownames(geneinfo) <- geneinfo$SYMBOL#意思应该是将geneinfo中SYMBOL的信息作为行名，SYMBOL里是基因名
  
  ## for inferCNV the geneorder file just have 3 colums: chr, star, end, there is no name colum
  geneinfo$SYMBOL <- NULL
  geneinfo = geneinfo %>% mutate(chr_n=strsplit2(geneinfo$chr, split='chr')[,2]) %>%  #mutate函数是新增一列
    mutate(chr_n = replace(chr_n, chr_n=='X', 23)) %>%#
    mutate(chr_n = replace(chr_n, chr_n=='Y', 24)) %>%
    mutate(chr_n = replace(chr_n, chr_n=='M', 25)) %>% 
    mutate(chr_n = as.numeric(chr_n)) %>%
    arrange(chr_n, start, end) %>% 
    select(-chr_n) 
  #head(geneinfo)
  
  seurat_i$cluster <- Idents(seurat_i) 
  metada <- seurat_i@meta.data
  #
  metadata <- metada[,c("merged_cluster", "merged_cluster")]
  colnames(metadata) <- c("cellType", "group")
  #head(metadata)
  
  
  write.table(metadata, file=paste0(workdir,'/',sample,"/cellAnnotations.txt"), sep="\t", col.names = FALSE,quote = FALSE)
  write.table(geneinfo, file=paste0(workdir,'/',sample,"/gene_ordering_file.txt"), sep = "\t", col.names = FALSE, quote = FALSE)
  
  count_mat <- mat_comb[rownames(geneinfo),]
  rm(mat_comb)
  rm(seurat_i)
  rm(metada)
  gc()
  
  infercnv_obj <- CreateInfercnvObject(raw_counts_matrix = count_mat, 
                                       gene_order_file = paste0(workdir,'/',sample,"/gene_ordering_file.txt"), 
                                       delim = "\t", min_max_counts_per_cell = c(10, +Inf),
                                       ref_group_names = c('Bcell','Endothelial','Macrophage','Tcell'),
                                       annotations_file = paste0(workdir,'/',sample,"/cellAnnotations.txt"))
  
  rm(count_mat)
  gc()
  save(infercnv_obj, file = paste0(workdir,'/',sample,"/infercnv_obj1.RData"))
 
  infercnv_obj = infercnv::run(infercnv_obj,
                               cutoff=0.1,  # use 1 for smart-seq, 0.1 for 10x-genomics,过滤低表达基因
                               out_dir=sample,  # dir is auto-created for storing outputs
                               cluster_by_groups=F,   # cluster
                               scale_data=FALSE,
                               denoise=T, #denoise：默认FALSE，对CNV矩阵进行降噪
                               sd_amplifier=1.5,  
                               HMM=T, #
                               output_format = "png",  
                               num_threads=50 #线程数
  )
  save(infercnv_obj, file = paste0(workdir,'/',sample,"/infercnv_obj2.RData"))

}

#### 运行函数
run_infercnv("P11")
run_infercnv("P18")
run_infercnv("P3A")
run_infercnv("P3B")
run_infercnv("P3C")





