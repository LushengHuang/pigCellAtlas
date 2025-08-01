
import numpy as np
import pandas as pd
import scanpy as sc
import os
import sys
from matplotlib.pyplot import rc_context
from math import exp 

#------------------ Load data ------------------#
mydata=sc.read("input_file.h5ad")

#------------------ UMAP visualization ------------------#
tissue="chorioallantoic membrane"
with rc_context({
    'figure.figsize': (5, 5)}):
    sc.pl.umap(mydata,color=['broadcelltype1'],legend_fontsize='small',legend_loc="on data",save='_harmony'+'_resolution_res0.5_'+tissue+'_broadcelltype1.pdf')