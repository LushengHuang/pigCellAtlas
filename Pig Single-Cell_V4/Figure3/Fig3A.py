

import numpy as np
import pandas as pd
import scanpy as sc
import os
import sys
from math import exp
import cosg as cosg
import importlib
importlib.reload(cosg)

wd = 'path'
os.chdir(wd)
sc.settings.verbosity = 3
sc.logging.print_header()
sc.settings.set_figure_params(dpi=600,facecolor='white')

# --------------- Extract endothelial cells from whole data ------------------------------------
dat = sc.read('/path/merge_405Samples_ct6bt6_20220408.h5ad')
celltypes = ['LEC','VEC','Endothelial cells']
idx = [(True if i in celltypes else False) for i in dat.obs.celltype6]
dat = dat[idx,:]

dat = dat[dat.obs.doublet_score<0.1,:]
dat = dat[dat.obs.n_genes<3000,:]

uniInd = set(dat.obs.Individual)
uniOrgan = set(dat.obs.Organ)

Individual = list(dat.obs.Individual)
Organ = list(dat.obs.Organ)

# --------- prepare marker genes for identifying the EC subtypes -----------------------
cap = ['CA4',"KDR",'RGCC','CD36','SPARC','SPARCL1','GPIHBP1','PPFIBP1','EPAS1','EDNRB']
art = ['TMEM100','SEMA3G','GJA5','GJA4','HEY1','BMX','CLU','EFNB2']
vein = ['ACKR1','ADGRG6','IL6','IL1R1','CCL26','VWF','ICAM1']
larv = ['VWF','ICAM1','VCAM1']
genv = ['FLT1','SOX17','GATA2','BMPR2']
prol = ['TOP2A','HMGB2','TPX2','ENSSSCG00000026302-KI67']
hem = ['HBB','AHSP','ALAS2']
general = ['CLDN5','PECAM1','CALCRL']
lym = ['CCL21','MMRN1','PROX1','RELN']
capArt = ['TGFB2','GLUL']
capVein = ['TFRC']
aero = ['EDNRB','TBX2']
other = ['ENSSSCG00000010992-AQP7','FABP5','RSPO3','WNT2','WNT9B','ENSSSCG00000000854-STAB2','CHST4','FUT7','ST6GAL1','NFATC1','ALCAM','NPR3','TBX20']
LSEC = ['ENSSSCG00000000854-STAB2']
PortalVein = ['RSPO3','WNT2','WNT9B']
CentralVein = ['VWF','GJA5']
HEV = ['CHST4','FUT7']

genes = {'General':general,'Lymphatic':lym,'Art':art,'Cap':cap,'Vein':vein,'HEV':HEV}


# ----------- expBatch: b1 and b2 are scRNA-seq in two batches, b3 is the snRNA-seq ----------------
# ----------- expBatch is used to correct for batch effect in the EC subclustering analysis --------
expBatch=['b2']*len(dat.obs)
Individual2 = list(set(dat.obs.Individual_2))
Individual2.remove("Sow2")
Individual2.remove("Sow2_E1")
Individual2.remove("Boar")

for i in range(len(dat.obs)):
    ind = dat.obs.Individual_2[i]
    platform = dat.obs.Platform[i]
    if ind in Individual2:
        expBatch[i] = 'b1'
    if platform == 'sn':
        expBatch[i] = 'b3'

dat.obs['expBatch'] = expBatch

dat.obs['ind_batch'] = [f"{ind}_{batch}" for ind, batch in zip(dat.obs['Individual'], dat.obs['expBatch'])]

# ------------  Perform clustering analysis by organs -------------
uniOrgan

for j in uniOrgan:
    v1 = [(True if Organ[k] == j else False) for k in range(len(Individual))]
    cellnum = v1.count(True)
    if cellnum > 100:
        tmp = dat[v1,:]
        filename = j+"_"+str(cellnum)
        if os.path.exists(wd+filename):
            continue 
        os.mkdir(wd+filename)
        os.chdir(wd+filename)
        #sc.pp.normalize_total(tmp, target_sum=1e4)
        #sc.pp.log1p(tmp)
        tmp.layers['scaled'] = sc.pp.scale(tmp,copy=True).X
        sc.pp.highly_variable_genes(tmp,n_top_genes=2000)
        sc.tl.pca(tmp,svd_solver='arpack')
        sc.pp.neighbors(tmp,n_neighbors=10)
        sc.external.pp.harmony_integrate(tmp,key='ind_batch')
        sc.pp.neighbors(tmp,n_neighbors=10,use_rep='X_pca_harmony')
        sc.tl.umap(tmp)
        sc.tl.louvain(tmp,resolution=1)
        uniqclt = list(set(tmp.obs.louvain))
        lst = list(tmp.obs.louvain)
        cltFrq = [lst.count(k) for k in uniqclt] 
        clt2keep = [uniqclt[k] for k in range(len(cltFrq)) if cltFrq[k] > 9]
        idx = [(True if k in clt2keep else False) for k in tmp.obs['louvain']]
        tmp = tmp[idx,:]
        sc.pl.dotplot(tmp,genes,groupby='louvain',swap_axes=False,standard_scale='var',dendrogram=False,save=filename+'.pdf')
        sc.pl.umap(tmp,color=['louvain'],ncols=1,save=filename+'_louvain.pdf',legend_loc="on data",legend_fontoutline=2,frameon=False)
        sc.pl.umap(tmp,color=['Platform'],  ncols=1,save=filename+'_platform.pdf',legend_fontoutline=2,frameon=False)
        sc.pl.umap(tmp,color=['TissueName'],ncols=1,save=filename+'_tissue.pdf',legend_fontoutline=2,frameon=False)
        sc.pl.umap(tmp,color=['Individual'],ncols=1,save=filename+'_Individual.pdf',legend_fontoutline=2,frameon=False)
        cltnum = len(set(tmp.obs.louvain))
        if cltnum == 1:
            continue
        sc.tl.rank_genes_groups(tmp,'louvain',method='wilcoxon')
        result = tmp.uns['rank_genes_groups']
        groups = result['names'].dtype.names
        markerGeneRank50 = pd.DataFrame({group + '_' + key[:1]: result[key][group] for group in groups for key in ['names','pvals']}).head(500)
        markerGeneRank50.to_excel(filename+'_markerGeneInfo_500.xlsx')
        pd.DataFrame(tmp.obs).to_csv(filename+'_cellInfo.csv')
        cosg.cosg(tmp,key_added='cosg',mu=1,n_genes_user=50,groupby='louvain')
        cosg_names = tmp.uns['cosg']['names']
        top5_genes_per_cluster = {
            cluster: list(cosg_names[cluster][:5])
            for cluster in cosg_names.dtype.names
        }
        sc.pl.dotplot(tmp,var_names=top5_genes_per_cluster,groupby="louvain",standard_scale='var',save=filename+"_celltype.cosg.png")
        sc.pl.dotplot(tmp,var_names=top5_genes_per_cluster,groupby="louvain",standard_scale='var',save=filename+"_celltype.cosg.pdf")
        tmp.write(filename+'.h5ad')




