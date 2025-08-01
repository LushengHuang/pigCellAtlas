import numpy as np
import pandas as pd
import scanpy as sc
import os
import sys
from matplotlib.pyplot import rc_context
from math import exp 
import cosg as cosg
import importlib
importlib.reload(cosg)







#------------------ Set parameters ------------------#
dat = dat[dat.obs.doublet_score<0.2,:]
dat = dat[dat.obs.pct_counts_mt<15,:]
res=2.5
sc.settings.verbosity = 3
sc.logging.print_header()
sc.settings.set_figure_params(dpi=1200,facecolor='white')



#------------------ Load data ------------------#

os.chdir('###')
dat=sc.read('10x_17samples_124635_cellinfoadded.h5ad')




#------------------ Preprocessing ------------------#
sc.pp.filter_genes(dat,min_cells=3)
sc.pp.highly_variable_genes(dat,n_top_genes=2000)
sc.tl.pca(dat,svd_solver='arpack')
rowdata1=dat


#------------------ Batch correction using BBKNN ------------------#
dat=rowdata1
sc.external.pp.bbknn(dat, batch_key='Individual_2')
sc.tl.umap(dat)
rowdata2=dat
sc.tl.leiden(dat, resolution = res)
platform='sn_10x_17sam_114400_batchI'


#------------------ UMAP visualization ------------------#
sc.pl.umap(dat,color="leiden",legend_loc="on data",legend_fontsize = 8,frameon=False,save='_'+'bbknn'+'_resolution_'+str(res)+'_'+platform+'_leiden_1.pdf')

sc.pl.umap(dat,color="Platform",save='_'+'bbknn'+'_resolution_'+str(res)+'_'+platform+'_Platform.png')
sc.pl.umap(dat,color="SampleID",save='_'+'bbknn'+'_resolution_'+str(res)+'_'+platform+'_SampleID.pdf')
sc.pl.umap(dat,color="TissueName",save='_'+'bbknn'+'_resolution_'+str(res)+'_'+platform+'_TissueName.pdf')

sc.pl.umap(dat,color="Group",save='_'+'bbknn'+'_resolution_'+str(res)+'_'+platform+'_Group.pdf')
sc.pl.umap(dat,color="S_I2",save='_'+'bbknn'+'_resolution_'+str(res)+'_'+platform+'_S_I2.pdf')


sc.pl.umap(dat,color="doublet_score",save='_'+'bbknn'+'_resolution_'+str(res)+'_'+platform+'_doublet_score.png')
sc.pl.umap(dat,color=['total_counts'],save='_'+'bbknn'+'_resolution_'+str(res)+'_'+platform+'_total_counts.png',frameon=False)
sc.pl.umap(dat,color=['pct_counts_mt'],save='_'+'bbknn'+'_resolution_'+str(res)+'_'+platform+'_pct_counts_mt.png',frameon=False)

#------------------Known marker gene visualization ------------------#
marker_genes_dict = {
 'B cells':['ENSSSCG00000038719-LV319','ENSSSCG00000036224-KVD28','CD79A'],
 'Cardiomyocytes': ['RYR2','DMD','PDE4D','PDE3A','SLC8A1'],
 'Endothelial cells':['PECAM1','VWF','CDH5','ID1'],
 'Fibroblasts':['DCN','GSN','COL1A2','COL3A1'],
 'lymEC':['CCL21','MMRN1','LYVE1'],
 'DC':['FLT3','IRF8','BCL11A','ITGAE'],
 'M1_MP':['SLA-DQB1','CD86','GPR183','CD80','FN1'],
 'M2_MP':['CD163','C1QB','C1QC','C1QA','MRC1','MERTK','LYVE1'],
 'Monocytes':['TBXAS1','ENSSSCG00000036618-FCG3A','CXCL2'],
 'NKT cells': ['GNLY','KLRK1','NKG7','GZMH'],
 'Pericytes': ['RGS5','MYO1B','EGFLAM','CPM','GUCY1A2'],
 'Schwann cells': ['SCN7A','CHL1','GRIK3'],
 'Smooth muscle cells': ['ACTA2','ITIH4','MYH11','MYLK'],
 'T cells': ['CD3D','CD3E','CD2','CCL5','PTPRC'],
 'Adipocytes': ['PLIN1','PPARG','ADIG','ADIPOQ']
 } 
 
sc.pl.dotplot(dat,marker_genes_dict,groupby='leiden',swap_axes=False,save=platform+'_resolution_'+str(res)+'_allmarker.pdf')
 


#------------------ Save processed data ------------------#
dat.write('bbknn'+'_resolution_'+str(res)+'_'+platform+'.h5ad')
pd.DataFrame(dat.obs).to_csv('bbknn'+'_resolution_'+str(res)+'_'+platform+'_cellInfo.csv')



#------------------ Marker gene detection ------------------#
sc.tl.rank_genes_groups(dat,'leiden',method='wilcoxon')
result = dat.uns['rank_genes_groups']
groups = result['names'].dtype.names
markerGeneRank200 = pd.DataFrame({group + '_' + key[:1]: result[key][group] for group in groups for key in ['names','pvals','logfoldchanges']}).head(200)
markerGeneRank200.to_csv('bbknn'+'_resolution_'+str(res)+'_'+platform+'_markerGeneInfo_200.csv')



#------------------ top marker gene detection ------------------#
sc.tl.rank_genes_groups(dat,'leiden',method='wilcoxon')
sc.pl.rank_genes_groups_dotplot(dat,n_genes=10,save=platform+'_resolution_'+str(res)+'_top10.png')

sc.pl.rank_genes_groups_dotplot(
    dat,
    n_genes=5,
    values_to_plot="logfoldchanges",
    cmap='bwr',
    vmin=-7,
    vmax=7,
    min_logfoldchange=1.5,
    colorbar_title='log fold change',
    save=platform+'_resolution_'+str(res)+'_top5_lgf1.pdf')
    
    
#------------------ top marker gene detection using COSG------------------#

cosg.cosg(dat,
    key_added='cosg',
    mu=1,
    n_genes_user=50,
    groupby='leiden')

sc.pl.rank_genes_groups_dotplot(dat,
    groupby='leiden',
    cmap='Spectral_r',
    standard_scale='var',
    n_genes=5,
    key='cosg',
    save=platform+'_resolution_'+str(res)+"_cosg.pdf")
	
