import numpy as np
import pandas as pd
import scanpy as sc

mydata=sc.read("Vagina_sampleNum_2_cellNum_16290.h5ad")
mydata.layers['scaled'] = sc.pp.scale(mydata, copy=True).X

markers = {'Endothelial cells LYVE1+':['PECAM1','CCL21','LYVE1'], 'Endothelial cells VWF+':'VWF',
       'Epithelial cells':['EPCAM'], 'Fibroblasts':['COL1A1','DCN'], 'Macrophages CD163+':'CD163',
       'Macrophages GPR183+':['GPR183'], 'Plasma cells':['MZB1','CD79A'], 'Proliferating Fibroblasts':['TOP2A','TPX2'],
       'Proliferating epithelial cells':'CDH1', 'Schwann cells':['SCN7A'],
       'Smooth muscle cells':['ACTA2','MYLK'], 'Smooth muscle cells RGS5+':'RGS5'}

sc.pl.umap(mydata,color=['celltype3'],frameon=False,legend_loc="on data",save="_celltype3.pdf")
sc.pl.dotplot(mydata, markers, groupby='celltype3', dendrogram=True,layer="scaled",save="_markers.pdf")
