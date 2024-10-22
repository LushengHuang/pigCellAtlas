#Inegrate Chorioallantoic Membrane data (5 NW vs 5 GR fetal pig )
import numpy as np
import pandas as pd
import scanpy as sc
import os
import sys
from matplotlib.pyplot import rc_context
from math import exp 

data=sc.read("chomem_harmony_5v5.h5ad")
tissue='10cho_mem'

ratio=1
res = 0.5
colorList1 = ['louvain']
sc.external.pp.harmony_integrate(data, key='batch_new')
sc.pp.neighbors(data,use_rep='X_pca_harmony')
sc.tl.umap(data)
sc.tl.louvain(data,resolution=res)
with rc_context({
    'figure.figsize': (5, 5)}):
    sc.pl.umap(data,color="louvain",legend_loc="on data",legend_fontsize='small',save='_harmony'+'_resolution_'+str(res)+'_'+tissue+'_louvain.png')

with rc_context({
    'figure.figsize': (5, 5)}):
    sc.pl.umap(data,color=['celltype'],legend_loc="on data",legend_fontoutline=1,frameon=False,size=4,legend_fontsize=8,title=['Maternal-fetal-interface'],save='_'+tissue+'_celltype.png')

with rc_context({
    'figure.figsize': (5, 5)}):
    sc.pl.umap(data,color=['celltype'],legend_fontoutline=1,frameon=False,size=4,legend_fontsize=8,title=['Maternal-fetal-interface'],save='_'+tissue+'_celltype_legend.png')

mydata=data
marker_genes_dict = {
'Proliferating plasma cells':['JCHAIN','TOP2A','HMGB2'],
'B cells':['CD79A','CD79B','CCR4'],
'NKT cells':['ENSSSCG00000000640-NKG2A','NKG7','GNLY','CCL5','PRF1'],
'T cells':['PTPRC','CD3D','CD8A','CD8B','CD4'],
'Neutrophils':['CD177','LTF','MMP8','MMP9','S100A12','S100A8','ENSSSCG00000006588-S10A9','CXCL8','RETN','NCF1','FCGR1A'],
'Monocytes':['ENSSSCG00000036618-FCG3A','TYROBP','CD14'],
'pDCs':['FLT3','IRF8','BLNK'],
'DCs':['SLA-DRA','SLA-DQB1','XCR1','ENSSSCG00000021576-CD83','RGS1','APOE'],
'Fibroblasts':['COL1A1','COL3A1','TBX3','PEG10','DCN','TNN'],
'Myfibro':['ACTA2','MYLK','MYH11','MYL9','RGS5','COL1A1','COL3A1'],
'Macrophages':['CD209','C1QB','C1QC','CD163','CD68','MRC1','MARCO','CD86'],
'Trophoblasts':['PERP','TACSTD2','KRT8','AQP3','UPK2','ADAM28','CYP3A29','SLC15A1','LRP2','CDH16','UPK1A','AGR2'],
'mast cells':['MS4A2','FCER1G','CD300C'],
'VCTs':['PAGE4','CYP19A1','CYP3A29'],
'EVTs':['MFAP5','DEFB1','SLC15A1'],
'VECs':['PECAM1','ACKR1','PLVAP','CCL21','CD63'],
'LECs':['LYVE1'],
'Epi':['WFDC2','EPCAM','CDH1','TFF2','ENSSSCG00000033382-MUC5A','C5','C6','ELOVL2','ELOVL6','MS4A2'],
}
sc.pl.dotplot(mydata, marker_genes_dict,'louvain',dendrogram=True,save=tissue+'_louvain_markergene.pdf')


sc.tl.rank_genes_groups(mydata,'louvain',method='wilcoxon',pts=True)
sc.pl.rank_genes_groups(mydata,n_genes=25,sharey=False,save='_'+tissue+'.png')
result = mydata.uns['rank_genes_groups']
groups = result['names'].dtype.names
markerGeneRank200 = pd.DataFrame({group + '_' + key[:1]: result[key][group] for group in groups for key in ['names','pvals','logfoldchanges']}).head(200)
markerGeneRank200.to_csv('/harmony_'+'_resolution_'+str(res)+'_'+tissue+'_markerGeneInfo_200.csv')
########################################################
dat=mydata
gene_ids = dat.var.index.values
clusters = dat.obs['louvain'].cat.categories
dat0 = dat
dat = dat[:,gene_ids].X.toarray()
dat = pd.DataFrame(dat,columns=gene_ids,index=dat0.obs.louvain)

def exp1(x):
    return exp(x)-1

dat = dat.applymap(exp1)
average_dat = dat.groupby(level=0).mean()
bool_dat = dat.astype(bool)
fraction_dat = bool_dat.groupby(level=0).sum()/bool_dat.groupby(level=0).count()
out = open(tissue+'_markergenes.txt','w')
out.write('tissueSample'+'\t'+'GeneName'+'\t'+'averageExprLevel'+'\t'+'celltype'+'\t')
out.write('expFraction'+'\t'+'R1R2ratio'+'\t'+'TSI'+'\t'+'Tau'+'\n')

n = len(average_dat)
for i in gene_ids:
    if average_dat[i].sum()==0:
        continue
    X = average_dat[i]
    PCT = fraction_dat[i]
    r1 = max(X)
    Y = list(X)
    maxIdx = Y.index(r1)
    pct = PCT[maxIdx]
    if pct < 0.05:
        continue
    Y.remove(r1)
    r2 = max(Y)
    TSI = str(r1/X.sum())
    x_hat = X/r1
    Tau = str((n-x_hat.sum())/(n-1)) 
    if r2 == 0:
        out.write('\t'+i+'\t'+str(r1)+'\t'+clusters[maxIdx]+'\t')
        out.write(str(pct)+'\t'+'99'+'\t'+TSI+'\t'+Tau+'\n')
        continue
    ratio = r1/r2
    out.write('\t'+i+'\t'+str(r1)+'\t'+clusters[maxIdx]+'\t')
    out.write(str(pct)+'\t'+str(ratio)+'\t'+TSI+'\t'+Tau+'\n')

out.close()


mydata.obs['celltype1']=mydata.obs.louvain
cellInfoFile = open("./celltype1.txt","r")
celltype1 = []
for row in cellInfoFile:
    cell = row.replace('\n','').split('\t')
    celltype1.append(cell[0])

cellInfoFile.close()
mydata.rename_categories('celltype1',celltype1)
cellTypeList = [i.split('_')[1] for i in mydata.obs.celltype1]
mydata.obs['celltype1'] = cellTypeList
with rc_context({
    'figure.figsize': (5, 5)}):
    sc.pl.umap(mydata,color=['celltype1'],legend_fontsize='small',legend_loc='on data',save='_harmony'+'_resolution_'+str(res)+'_'+tissue+'_celltype1.pdf')



#Add Infomation
mydata=sc.read("harmony_resolution_0.5_10cho_mem_celltype1_group.h5ad")
mydata.obs['Family'] = mydata.obs['SampleID']
mydata.obs['Family'] = mydata.obs['Family'].str.replace('LA230709066','335129')
mydata.obs['Family'] = mydata.obs['Family'].str.replace('LA230709079','335129')
mydata.obs['Family'] = mydata.obs['Family'].str.replace('LA230903019','333588')
mydata.obs['Family'] = mydata.obs['Family'].str.replace('LA230903021','333588')
mydata.obs['Family'] = mydata.obs['Family'].str.replace('LA230903029','331941')
mydata.obs['Family'] = mydata.obs['Family'].str.replace('LA230709063','335129')
mydata.obs['Family'] = mydata.obs['Family'].str.replace('LA230709069','335129')
mydata.obs['Family'] = mydata.obs['Family'].str.replace('LA230903015','333588')
mydata.obs['Family'] = mydata.obs['Family'].str.replace('LA230903017','333588')
mydata.obs['Family'] = mydata.obs['Family'].str.replace('LA230903027','331941')


mydata.write('/harmony'+'_resolution_'+str(res)+'_'+tissue+'_celltype1_group.h5ad')
pd.DataFrame(mydata.obs).to_csv('/harmony'+'_resolution_'+str(res)+'_'+tissue+'_celltype1_group_cellInfo.csv')

mydata.obs['broadcelltype1']=mydata.obs.louvain
cellInfoFile = open("./broadcelltype1.txt","r")
broadcelltype1 = []
for row in cellInfoFile:
    cell = row.replace('\n','').split('\t')
    broadcelltype1.append(cell[0])

cellInfoFile.close()
mydata.rename_categories('broadcelltype1',broadcelltype1)
cellTypeList = [i.split('_')[1] for i in mydata.obs.broadcelltype1]
mydata.obs['broadcelltype1'] = cellTypeList

with rc_context({
    'figure.figsize': (5, 5)}):
    sc.pl.umap(mydata,color=['broadcelltype1'],legend_fontsize='small',legend_loc="on data",save='_harmony'+'_resolution_'+str(res)+'_'+tissue+'_broadcelltype1.pdf')

mydata.write('/harmony'+'_resolution_'+str(res)+'_'+tissue+'_broadcelltype1_group.h5ad')
pd.DataFrame(mydata.obs).to_csv('/harmony'+'_resolution_'+str(res)+'_'+tissue+'_broadcelltype1_group_cellInfo.csv')

import diopy
dat=sc.read('/harmony'+'_resolution_'+str(res)+'_'+tissue+'_broadcelltype1_group.h5ad')
diopy.output.write_h5(dat, file ='/harmony'+'_resolution_'+str(res)+'_'+tissue+'_broadcelltype1_group.h5')
