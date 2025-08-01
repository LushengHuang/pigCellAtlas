import numpy as np
import pandas as pd
import scanpy as sc
import os
import sys
from matplotlib.pyplot import rc_context
from math import exp 

# --------------------- Load dataaset ---------------------#
data=sc.read("input_files.h5ad")
dat0=data
dat0.obs['Family_SampleID_group_broadcelltype1_index'] = dat0.obs['Family'].astype(str) +'-'+dat0.obs['SampleID'].astype(str) +'-'+ dat0.obs['group'].astype(str) +'-'+dat0.obs['broadcelltype1'].astype(str) +'-'+dat0.obs.index
dat1 = dat0

#--------------------- Calculate gene expression ---------------------#
gene_ids =  ['SLC1A5','SLC38A2']
dat1 = dat1[:,gene_ids].X.toarray()
dat1 = pd.DataFrame(dat1,columns=gene_ids,index=dat0.obs.Family_SampleID_group_broadcelltype1_index)
def exp1(x):
    return exp(x)-1

dat1 = dat1.applymap(exp1)

average_dat = dat1.groupby(level=0).mean()

#--------------------- Save the result file ---------------------#
pd.DataFrame(average_dat).to_csv('output_file.csv')