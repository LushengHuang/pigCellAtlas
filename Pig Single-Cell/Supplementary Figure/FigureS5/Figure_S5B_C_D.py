import numpy as np
import pandas as pd
import scanpy as sc
import scanpy.external as sce
import os
import sys
import anndata
import matplotlib.pyplot as plt
import seaborn as sns 

# Fig. S5B and C
adata_ref = sc.read("discovery_sow_adipose_EC_15667.h5ad")
adata = sc.read("validaton_sow_adipose_EC_11303.h5ad")

adata_ref.obs['Cohort']='discovery'
adata.obs['Cohort']='validaton'


var_names = adata_ref.var_names.intersection(adata.var_names)
adata_ref = adata_ref[:, var_names].copy()
adata = adata[:, var_names].copy()
adata.X = adata.X.toarray()
sc.tl.ingest(adata, adata_ref, obs="celltype8")
adata_concat = anndata.concat([adata_ref, adata], label="Cohort", keys=["discovery", "validaton"]) 

adata_concat.obs["celltype8"] = (
    adata_concat.obs["celltype8"]
    .astype("category")
    .cat.reorder_categories(adata_ref.obs["celltype8"].cat.categories)
)
del adata_concat.obs['n_genes_by_counts']
adata_concat.write("adata_concat_2.h5ad")

sc.pl.umap(adata_ref,color=['celltype8'],save='_celltype8.pdf',legend_loc="right margin",legend_fontoutline=2,frameon=False)
sc.pl.umap(adata_concat,color=['Cohort'],save='_Cohort.pdf',legend_loc="right margin",legend_fontoutline=2,frameon=False)



# Fig. S5D
df = adata_concat.obs[["celltype8",'Cohort']]
pivot_table = df.pivot_table(index='celltype8', columns='Cohort', aggfunc='size', fill_value=0)
pivot_table = pivot_table.div(pivot_table.sum(axis=1), axis=0)
plt.figure(figsize=(10, 6))
pivot_table.plot(kind='bar', stacked=True, colormap='tab20', edgecolor='black')
plt.title('Proportion of Cohort in Each celltype8', fontsize=14)
plt.xlabel('celltype8', fontsize=12)
plt.ylabel('Proportion', fontsize=12)
plt.legend(title='Cohort', bbox_to_anchor=(1.05, 1), loc='upper left')
plt.tight_layout()

plt.savefig('celltype8_project_distribution.pdf', format='pdf', dpi=300, bbox_inches='tight')
plt.close()

