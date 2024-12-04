import os,sys
del sys.path[4]
import numpy as np
import pandas as pd
import scanpy as sc

sc.settings.verbosity = 3
sc.logging.print_header()
sc.settings.set_figure_params(dpi=100,dpi_save=300,facecolor='white',frameon=False,fontsize=10,vector_friendly=True,figsize=(3,3))

def removeBiasGenes(adata):
    IGgenes = adata.var_names.str.startswith('IG')
    COLgenes = adata.var_names.str.startswith('COL')
    malat1 = adata.var_names.str.startswith('MALAT1')
    MTgenes = adata.var_names.str.startswith('MT')
    hb_genes = adata.var_names.str.contains('^HB[^(P)]')
    RPgenes = adata.var_names.str.startswith('RP') & adata.var_names.str.contains('-')
    RPgenes2 = adata.var_names.str.contains('^RP[SL]')
    CTCgenes = adata.var_names.str.startswith('CTC') & adata.var_names.str.contains('-')
    MIRgenes = adata.var_names.str.startswith('MIR')
    ACgenes = adata.var_names.str.contains('^AC[0-9]') & adata.var_names.str.contains('.')
    CTgenes = adata.var_names.str.startswith('CT') & adata.var_names.str.contains('-')
    LINCgenes = adata.var_names.str.contains('^LINC[0-9]')
    ALgenes = adata.var_names.str.contains('^AL') & adata.var_names.str.contains('.')
    # KRTgenes = adata.var_names.str.startswith('KRT')

    MESgenes = adata.var_names.isin(['ACTA2','ACTG2','TAGLN','DES','ACTC1','LUM','BGN','DCN','POSTN','MRC1',"CALD1","TPM1","MGP","C1S","PDPN",
                                      'FLT1','PLVAP','SPARCL1','PECAM1','CD31','HSPG2','VWF','CDH5','SELE','VCAM1','ENG','PDGFRB'])
    Tcellgenes = adata.var_names.isin(['TRAC','SPN','TAGAP','IL7R','PTPRC','IL10RA','LTB','RGS1','TRBC2','LCP2','FCMR','IL2RB','NLRC3','IL21R','IL18R1','SLAMF1','LAT2','CD2','CD52','KLRB1',
                                        "CD3D", "CD3E", "CD3G","CD8A", "CD8B", "GZMK", "CD4","TNFRSF4"])
    Plasmocytegenes = adata.var_names.isin(['FAM30A','MZB1','FCRL5','JCHAIN','LAX1','THEMIS2','SLAMF7','LY9','JSRP1','NCKAP1L'])
    Macrogenes = adata.var_names.isin(['LILRB5','MPEG1','MS4A7','FCGR2A','LYVE1','SIGLEC1','SIRPB2','RGS1','THEMIS2','ITGAX','C1QA','KCNE1','TYROBP','CSF3R','C1QB','CNR2','ADGRE1',
                                        'CD163','SLC11A1','APOC1','FCER1G','FCGR3A'])
    Bcellgenes = adata.var_names.isin(['MS4A1','BLK','CD79B','TAGAP','TNFRSF13C','FCMR','P2RX5','TLR10','FCRL1','CYBB','SCIMP','IRF8','LILRB1','THEMIS2','CR2','FCRL2','FCRL5','CD19','CD21','CD79A','CD79B','BLNK'])
    remove_genes = malat1 | MTgenes | hb_genes | RPgenes | RPgenes2 | CTCgenes | MIRgenes | ACgenes | CTgenes | LINCgenes | ALgenes | IGgenes | COLgenes | MESgenes | Tcellgenes | Plasmocytegenes | Macrogenes | Bcellgenes
    keep = np.invert(remove_genes)
    res = adata[:,keep]
    return res

os.chdir('/jdfsbjcas1/ST_BJ/P21H28400N0232/xieguixiang/BLCA/Union_allsample/epi_bigcell_15nei')
samplelist = [
    "HBCP1", "HBCP3A", "HBCP3B", "HBCP3C","HBCP5A","HBCP5B","HBCP6","HBCP7","HBCP8","HBCP9","HBCP10","HBCP11","HBCP12","HBCP13","HBCP14","HBCP15","HBCP17","HBCP18"]
st_ad_list = []
for sample in samplelist:
    adata = sc.read_h5ad(f'/jdfsbjcas1/ST_BJ/P21H28400N0232/xieguixiang/BLCA/{sample}/Spoint/image_extract/image_anno_rawdata.h5ad')
    # adata = adata[adata.obs['merged_cluster'].isin(['MIBC','Metastatic_tumor','NMIBC'])]
    st_ad_list.append(adata)
adata = sc.AnnData.concatenate(*st_ad_list,
                                join='outer',
                                batch_categories=samplelist,
                                batch_key='batch')

adata = removeBiasGenes(adata)

sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
print(adata.shape)

# # 找到每个细胞中是否存在NaN
# nan_cells = np.isnan(adata.X).any(axis=1)
# # 过滤掉包含NaN的细胞
# adata = adata[~nan_cells].copy()
# print(adata.shape)

sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
adata = adata[:, adata.var.highly_variable]
# adata.var_names = adata2.var_names[adata.var.highly_variable.index.to_numpy().astype(int)]
# 回归每个细胞的总计数和表达的线粒体基因的百分比的影响.
# sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt'])
# 将每个基因缩放到单位方差,阈值超过标准偏差10.
sc.pp.scale(adata, max_value=10)


sc.tl.pca(adata, svd_solver='arpack')
sc.pp.neighbors(adata, n_neighbors=15, n_pcs=40)
sc.tl.leiden(adata, resolution=0.6)
sc.tl.umap(adata)

adata.write_h5ad('epi_bigcell_umap.h5ad')
sc.pl.umap(adata, color='leiden', use_raw=False, save='epi_bigcell_umap.pdf')

# sc.pl.umap(adata, color='KRT5', use_raw=False, save='KRT5.pdf')
# sc.pl.umap(adata, color='KRT18', use_raw=False, save='KRT18.pdf')
sc.pl.umap(adata, color='KRT17', use_raw=False, save='KRT17.pdf')  ###选择基因画表达图

