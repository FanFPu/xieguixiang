import os,sys
# del sys.path[4]
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib import rcParams


# sc.settings.verbosity = 3
# sc.logging.print_header()
# sc.settings.set_figure_params(dpi=100,dpi_save=300,facecolor='white',frameon=False,fontsize=10,vector_friendly=True,figsize=(3,3))

os.chdir('/jdfsbjcas1/ST_BJ/P21H28400N0232/xieguixiang/BLCA/SC/R28BC/all/40pc/figures/')

adata = sc.read_h5ad('/jdfsbjcas1/ST_BJ/P21H28400N0232/xieguixiang/BLCA/SC/R28BC/all/40pc/humanSC_allgenes.h5ad')
# print(adata.var.highly_variable.index)
# adata2 = sc.read_h5ad('/jdfsbjcas1/ST_BJ/P21H28400N0232/fanfengpu/BC/SingleCell/R28BC/Data/seurat_merge/seurat_merge_new_haromony_subtype.h5ad')
# adata.var_names = adata2.var_names[adata.var.highly_variable.index.to_numpy().astype(int)]
print(adata.obs.columns)  # 打印数据中可用的列名
print(adata.obs['Patients'].unique())  # 打印 'Patients' 列的唯一值
# sns_palette = sns.color_palette("Set1", n_colors=len(adata.obs['cellTypes_new'].unique()))
# adata.obs['cellTypes_new'] = adata.obs['cellTypes_new'].astype('str')
# adata = adata[adata.obs['cellTypes_new'] != 'Undefined']
# custom_palette = ["#B87A3D", "#8A2BE2", "#FEC643", "#FF99FF","#679966", "#663300", "#C71585", "#00FF00", "#FF0000", "#CCCC00","#43D9FE", "#437BFE", "#333399", "#43FE69", "#E5C494"]
# sc.pl.umap(adata, color='cellTypes_new', palette=custom_palette, save='umap.pdf')
# sc.pl.umap(adata, color='KRT17', save='KRT17.pdf')
# sc.pl.umap(adata, color='KRT5', save='KRT5.pdf')
# 设置默认图像大小
# 设置画布大小


# 定义颜色和标签
color_dict = {"HBCP1": "#ed1299", "HBCP2": "#09f9f5", "HBCP3A": "#246b93", "HBCP3B": "#cc8e12", "HBCP3C": "#d561dd", 
    "HBCP4A": "#c93f00", "HBCP4B": "#ddd53e", "HBCP5A": "#4aef7b", "HBCP5B": "#e86502", "HBCP6": "#9ed84e",
    "HBCP7": "#39ba30", "HBCP8": "#6ad157", "HBCP9": "#8249aa", "HBCP10": "#99db27", "HBCP11": "#e07233",
    "HBCP12": "#ff523f", "HBCP13": "#ce2523", "HBCP14": "#f7aa5d", "HBCP15": "#cebb10", "HBCP16": "#03827f",
    "HBCP17": "#931635", "HBCP18": "#373bbf", "HBCP19": "#a1ce4c", "HBCP20": "#ef3bb6", "HBCP21": "#d66551", 
    "HBCP22": "#1a918f"}

# 提取adata.obs['Patients']中所有唯一标签，并将其映射到color_dict中
patients = adata.obs['Patients'].unique().tolist()

# 根据患者标签顺序生成颜色列表
color_list = [color_dict[patient] for patient in patients]

# 设置图像大小
plt.rcParams['figure.figsize'] = [8, 8]

# 绘制UMAP图并应用颜色列表
sc.pl.umap(adata, color='Patients', use_raw=False, palette=color_dict, save='patients.pdf')

# 保存图像
plt.savefig('patients.svg', dpi=300)