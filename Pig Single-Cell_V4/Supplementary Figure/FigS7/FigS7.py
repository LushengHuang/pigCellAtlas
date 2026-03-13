

outdir  = '~/17_withinOrgan_EC/res1.0'
datdir = '~/17_withinOrgan_EC/res1.0'

import numpy as np
import pandas as pd
import scanpy as sc
import os
import sys
from math import exp
import cosg as cosg
import importlib
importlib.reload(cosg)
from pathlib import Path

wd = '~/17_withinOrgan_EC/res1.0/Lymph nodes_9772/'
os.chdir(wd)
os.listdir('.')
sc.settings.verbosity = 3
sc.logging.print_header()
sc.settings.set_figure_params(dpi=600,facecolor='white')

h5ad_files = list(Path(".").glob("*.h5ad"))
h5ad_files_list = list(map(str, h5ad_files)) 

infile = 'Lymphnode'

dat = sc.read(h5ad_files_list[0])


dat.obs['celltype7'] = dat.obs.louvain
cellInfoFile = open("subtypes.txt","r")
celltype = []
for row in cellInfoFile:
    cell = row.replace('\n','').split('\t')
    celltype.append(cell[0])

cellInfoFile.close()

dat.rename_categories('celltype7',celltype)
cellTypeList = [i.split('-')[1] for i in dat.obs.celltype7]
dat.obs['celltype8'] = cellTypeList

dat = dat[dat.obs.celltype8 != 'Unknown',:]


sc.pl.umap(dat,color=['celltype8'],legend_fontoutline=1,legend_loc="on data",frameon=False,legend_fontsize=8,save='_'+infile+'_in2.pdf')
sc.pl.umap(dat,color=['celltype8'],legend_fontoutline=2,frameon=False,legend_fontsize=8,save='_'+infile+'_out2.pdf')

ct8colors = pd.DataFrame(dat.uns['celltype8_colors'])
ct8colors.to_csv(infile+'_ct8color.csv')

# calculate the differetial expressed genes
sc.tl.rank_genes_groups(dat,'celltype8',method='wilcoxon')
result = dat.uns['rank_genes_groups']
groups = result['names'].dtype.names
markerGeneRank50 = pd.DataFrame({group + '_' + key[:1]: result[key][group] for group in groups for key in ['names','pvals']}).head(500)
markerGeneRank50.to_excel(infile+'_markerGeneInfo_ct8_500.xlsx')
pd.DataFrame(dat.obs).to_csv(infile+'_cellInfo.csv')

# calculate the marker gene (rank1/rank2 > 4, pct > 0.05)
# Find cluster specfic genes (Tau and TSI values)
adata = dat

dat = adata
clusters = dat.obs['celltype8'].cat.categories
gene_ids = dat.var.index.values
dat = dat[:,gene_ids].X.toarray()
dat = pd.DataFrame(dat,columns=gene_ids,index=adata.obs['celltype8'])

dat = np.exp(dat) - 1
average_dat = dat.groupby(level=0).mean()
bool_dat = dat.astype(bool)
fraction_dat = bool_dat.groupby(level=0).sum()/bool_dat.groupby(level=0).count()
out = open(infile+'_specific_genes.txt','w')
out.write('genename'+'\t'+'average'+'\t'+'celltype'+'\t'+'TSI'+'\t'+'Tau'+'\t'+'fraction'+'\n')

n = len(average_dat)
for i in gene_ids:
    if average_dat[i].sum()==0:
        continue
    genename = i
    clust = average_dat[i].idxmax()
    average = str(average_dat[i][clust])
    fraction = str(fraction_dat[i][clust])
    TSI = str(float(average)/float(average_dat[i].sum()))
    x_hat = average_dat[i]/average_dat[i].max()
    Tau = str((n-x_hat.sum())/(n-1))
    out.write(genename+'\t'+average+'\t'+clust+'\t'+TSI+'\t'+Tau+'\t'+fraction+'\n')

out.close()

adata.write(infile+'.update.h5ad')

targetgenes = [gene for gene in adata.var_names if 'STAB2' in gene.upper()]
targetgenes

marker_genes_dict = {
    'General EC':['PECAM1','CDH5','CLDN5'], 
    'Artery':['TMEM100','GJA5','HEY1'],
    'Capillary':['CA4','KDR','RGCC','CLDN6'],
    'Lymphatic':['CCL21','MMRN1','PROX1','RELN'],
    'Vein':['ACKR1','ADGRG6','IL6','VWF'],
    'HEV':['CHST4','FUT7'],
    'Lym CLDN11+':['CLDN11'],
    'Lym STAB2+':['ENSSSCG00000000854-STAB2'],
    'Vein FCER1G+':['FCER1G'],
    'Vein IL6+':['IL6']
    }

sc.pl.dotplot(adata,
              var_names=marker_genes_dict,
              groupby="celltype8",
              standard_scale='var',
              dendrogram=False,
              swap_axes=False,
              save="_"+infile+"_markers.pdf")