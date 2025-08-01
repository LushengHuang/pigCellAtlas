###Subclustering for type II myofibers

import numpy as np
import pandas as pd
import scanpy as sc
import os
import sys
from matplotlib.pyplot import rc_context
from math import exp

wd1='./03.harmony_4v4QC/'
tissue = '8LD_TypeIImyofiber'
os.chdir(wd1)

sc.settings.verbosity = 3
sc.logging.print_header()
sc.settings.set_figure_params(dpi=600,facecolor='white')

data = sc.read("Normalized_celltype1-2_celltype11_8LD.h5ad")
mydata=data[data.obs['celltype1'].isin(['Type II myofibers'])] 

#Preprocessing
sc.pp.highly_variable_genes(mydata,n_top_genes=2000)
sc.tl.pca(mydata,svd_solver='arpack')

#harmony
sc.external.pp.harmony_integrate(mydata, 'SampleID')

sc.pp.neighbors(mydata,use_rep='X_pca_harmony')
sc.tl.umap(mydata)
res = 0.5
sc.tl.louvain(mydata, resolution = res)

#visualization
sc.pl.umap(mydata,color=['louvain'],legend_loc="on data",legend_fontsize = "xx-small",save='_res_'+str(res)+'_harmony_'+tissue+'_ondata'+'.pdf')
pd.DataFrame(mydata.obs).to_csv('Normalized_8LD_TypeIImyofiber_louvain_cellInfo.csv')
mydata.write('Normalized_8LD_TypeIImyofiber_louvain.h5ad')

#Find marker genes
sc.tl.rank_genes_groups(mydata,'louvain',method='wilcoxon')
result = mydata.uns['rank_genes_groups']
groups = result['names'].dtype.names
groups
markerGeneRank500 = pd.DataFrame({group + '_' + key[:1]: result[key][group] for group in groups for key in ['names','pvals','logfoldchanges']}).head(500)
markerGeneRank500.to_csv('Normalized_8LD_TypeIImyofiber_louvain_markerGeneInfo_500.csv')

gene = {'markes':['ATP2A1','MYLPF','TNNT3','TNNC2','TNNI2'],
        'DEG':['CTNNA3','LRRTM3','PDE4B','KCNQ5','BMPR1B',
               'EGF','NRAP','B3GALT1','PEBP4','FAM13C']}  
#dotplot
sc.pl.dotplot(mydata, gene,'louvain',dendrogram=True, save=tissue+'_louvain_markers.pdf')