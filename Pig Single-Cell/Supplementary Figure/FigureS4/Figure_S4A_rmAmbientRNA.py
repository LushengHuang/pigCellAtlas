import numpy as np
import pandas as pd
import scanpy as sc
import os
import sys
import dropkick as dk

wd = '/home/SCell/yangbin/03_cellCrossTissue/'
sc.settings.verbosity = 3
sc.logging.print_header()
sc.settings.set_figure_params(dpi=600,facecolor='white')
datfile = '/home/SCell/SC_Singleron/05.ATLAS2021/04.Merge405_0815/CellTypeAdded0828/merge_405samples_cellInfo_added_0908.h5ad'
dat = sc.read(datfile)

cellType = 'endothelial'
outdir = wd+cellType
#os.mkdir(outdir)
os.chdir(outdir)
cell2use = ['Endothelial cells']
idx = [(True if i in cell2use else False) for i in dat.obs['broadcelltype1']]
dat2 = dat[idx,:]
del dat

# 1 Find sample with > 20 ECs
uniqSampleID = list(set(dat2.obs.SampleID))
lst = list(dat2.obs.SampleID)
sampleFrq = []
for i in uniqSampleID:
    sampleFrq.append(lst.count(i))

sample2keep = []
for i in range(len(sampleFrq)):
    if sampleFrq[i] > 20:
        sample2keep.append(uniqSampleID[i])

del dat2

pd.DataFrame(sample2keep).to_csv("EC_sample.csv")

# 2 Find the Ambient RNA list 
dat = sc.read("/home/SCell/yangbin/rawData/merge_405Samples_raw1001.h5ad")
nkeep = 50
ambientList = []
for i in sample2keep:
    datSample = dat[dat.obs.SampleID==i,:]
#    datSample_dk = dk.dropkick(datSample,n_jobs=5)
#    datSample_empty = datSample[datSample.obs.dropkick_score<0.5,:]
#    datSample_empty = dk.recipe_dropkick(datSample_empty,n_hvgs=None, X_final="raw_counts")
    datSample = dk.recipe_dropkick(datSample,n_hvgs=None, X_final="raw_counts")
    tmplist = list(datSample.var_names[range(nkeep)])
    ambientList = ambientList+tmplist

pd.DataFrame(ambientList).to_csv('ambientRNA_all.csv')
#ambientListUniq = list(np.unique(ambientList))
#pd.DataFrame(ambientListUniq).to_csv("sn_ambientRNA.csv")

del dat
del datSample

# 3 Fine the top EC marker genes in each sample
dat = sc.read("/home/SCell/SC_Singleron/05.ATLAS2021/04.Merge405_0815/CellTypeAdded0828/merge_405samples_cellInfo_added_0908.h5ad")
mkrList = []
for i in sample2keep:
    datTmp = dat[dat.obs.SampleID==i,:]
    sc.tl.rank_genes_groups(datTmp,'celltype1',groups=['Endothelial cells'],method='wilcoxon')
    tmplist = pd.DataFrame(datTmp.uns['rank_genes_groups']['names']).head(nkeep)
    mkrList = mkrList+list(tmplist['Endothelial cells'])

pd.DataFrame(mkrList).to_csv('mkrGenes_all.csv')
#mkrListUniq = list(np.unique(mkrList))
#pd.DataFrame(mkrListUniq).to_csv("sn_mkrGenes.csv")
del dat




