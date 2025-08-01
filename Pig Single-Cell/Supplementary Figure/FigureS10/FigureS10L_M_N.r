# ingest--------------------------

import numpy as np
import pandas as pd
import scanpy as sc
import scanpy.external as sce
import os
import sys
import anndata

outdir = 'res0.3'
os.chdir(outdir)
sc.settings.verbosity = 3
sc.logging.print_header()
sc.settings.set_figure_params(dpi=600,facecolor='white')

#Pig type II myofiber subclusters
adata_ref = sc.read("Normalized_8LD_TypeIImyofiber_louvain.h5ad")
adata_ref.obs.rename(columns={"louvain": "celltype"}, inplace=True)
adata_ref.obs['Project']='Pig'

#Human type II myofibers
dat1=sc.read("harmony_sampleID_HumanMuscleCell_res0.3_celltype.h5ad")
# FigureS10L
sc.pl.umap(mydata,color=['celltype'],legend_loc="on data",legend_fontsize = "xx-small",save=tissue+'_celltype_ondata.pdf')
adata=dat1[dat1.obs.celltype.isin(['Type II myofiber'])]
adata.obs['Project']='Human'

var_names = adata_ref.var_names.intersection(adata.var_names)
adata_ref = adata_ref[:, var_names].copy()
adata = adata[:, var_names].copy()
adata.X = adata.X.toarray()
sc.tl.ingest(adata, adata_ref, obs="celltype")

adata_concat = anndata.concat([adata_ref, adata], label="Project", keys=["Pig", "Human"]) 
adata_concat.obs["celltype"] = (
    adata_concat.obs["celltype"]
    .astype("category")
    .cat.reorder_categories(adata_ref.obs["celltype"].cat.categories)
)

# FigureS10M-N
sc.pl.umap(adata_concat, color=["Project", "celltype"],save="Project_TypeIIsubtype_merge.pdf")