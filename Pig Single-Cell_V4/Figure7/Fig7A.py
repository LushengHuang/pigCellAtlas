import numpy as np
import pandas as pd
import scanpy as sc
import os
import sys
from matplotlib.pyplot import rc_context
from math import exp 
import diopy


mydata = sc.read('HCL_A_selOrg.h5ad')
data = pd.read_csv(r'human_gene_retain.csv',sep=',')
mydata = mydata[:,data.human]

mydataP = diopy.input.read_h5(file = "PCA_A_selOrg.h5")
data = pd.read_csv(r'pig_to_human_gene.csv',sep=',')
mydataP = mydataP[:,data.piggene]
mydataP.var_names=data.humangene

mydataM = sc.read('MCA_A_selOrg.h5ad')
data = pd.read_csv(r'mouse_to_human_gene.csv',sep=',')
mydataM = mydataM[:,data.gene]
mydataM.var_names=data.humangene

dat_merge = mydataP.concatenate(mydata,mydataM)

cellty = 'HMP_A_merge_homo'
mydata = dat_merge
mydata.obs['celltype']=mydata.obs['orig.ident']
mydata.obs['celltype']=mydata.obs['celltype'].str.replace('A-','')
mydata.obs['celltype']=mydata.obs['celltype'].str.replace('P-','')
mydata.obs['celltype']=mydata.obs['celltype'].str.replace('H-','')
mydata.obs['celltype']=mydata.obs['celltype'].astype('category')
mydata.var_names_make_unique()
sc.pp.filter_genes(mydata,min_cells=3)
sc.pp.normalize_total(mydata,target_sum = 1e4) # normalization
sc.pp.log1p(mydata) # log transform
sc.pp.highly_variable_genes(mydata,n_top_genes=2000)
sc.tl.pca(mydata,svd_solver='arpack')

# Integrate harmony and bbknn
sc.external.pp.harmony_integrate(mydata,'species')
ratio = 1
sc.external.pp.bbknn(mydata, batch_key='tissue',use_rep='X_pca_harmony',set_op_mix_ratio=ratio)
sc.tl.umap(mydata)
res = 1
sc.tl.leiden(mydata,resolution=res)
sc.pl.umap(mydata,color=['celltype'],legend_loc="on data",frameon=False,legend_fontsize=15,legend_fontoutline=4,size =20.,save=cellty+'_harmony_bbknn'+'_resolution_'+str(res)+'_leiden'+'.png')
sc.pl.umap(mydata,color=['species'],legend_loc="right margin",frameon=False,legend_fontsize=5,legend_fontoutline=2,size =50.,save=cellty+'_harmony_bbknn_'+'species'+'.png')
sc.tl.louvain(mydata,resolution=res)
mydata.write('harmony_bbknn'+cellty+'.h5ad')

#########markgene
sc.tl.rank_genes_groups(mydata,'leiden',method='wilcoxon')
result = mydata.uns['rank_genes_groups']
groups = result['names'].dtype.names
markerGeneRank = pd.DataFrame({group + '_' + key[:1]: result[key][group] for group in groups for key in ['names','pvals']}).head(1000)
markerGeneRank.to_csv(cellty+'_leiden_markerGeneInfo_1000.csv')
