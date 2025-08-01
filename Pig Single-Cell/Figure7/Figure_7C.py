import os
import pandas as pd
import numpy as np
from scipy.io import mmread
import scipy.sparse as sp
import matplotlib.pyplot as plt
from IPython.display import Image
import scanpy as sc
from cnmf import cNMF

np.random.seed(6)
# ------------------------- load data and QC --------------------------
adata = sc.read("harmony_bbknn_772062cell_HMP_A_merge_homo.h5ad")
adata.var_names_make_unique()
sc.pp.filter_cells(adata, min_genes=300)
sc.pp.filter_genes(adata, min_cells=10)
counts = adata.obs['species'].value_counts()


# --------------------------- start the cNMF analysis -----------------
numiter = 200
numhvgenes = 2000

output_directory = 'cNMF'
if not os.path.exists(output_directory):
    os.mkdir(output_directory)
    
adata.write("threeSpecies.h5ad")
run_name = 'CrossSpecial_cNMF'

# Specify the Ks to use as a space separated list in this case '5 6 7 8 9 10'
seed = 14
countfn = 'threeSpecies.h5ad' # path to the filtered counts dataset we output previously 
cnmf_obj = cNMF(output_dir=output_directory, name=run_name)
cnmf_obj.prepare(counts_fn=countfn, components=np.arange(8,20), n_iter=50, seed=14, num_highvar_genes=2000)
cnmf_obj.factorize(worker_i=0, total_workers=20)
cnmf_obj.combine()
cnmf_obj.k_selection_plot(close_fig=False)
print("This saves the corresponding figure to the following file: %s" % cnmf_obj.paths['k_selection_plot'])

selected_K = 17 # K value = 17 maximized the stability score 
density_threshold = 0.10
cnmf_obj.consensus(k=selected_K, density_threshold=density_threshold, show_clustering=True, close_clustergram_fig=False)

# We are most interested in the usage and gene_spectra_score files for the density threshold of 0.1
adata = sc.read(countfn)
hvgs = open('/cNMF/CrossSpecial_cNMF/CrossSpecial_cNMF.overdispersed_genes.txt').read().split('\n')
sc.pp.normalize_per_cell(adata,counts_per_cell_after=10**4)

# Set log-normalized data to the raw attribute of the AnnData Oject to make it easy to plot expression levels of individuals
# This does not log normalize the actual AnnData data matrix
adata.raw = sc.pp.log1p(adata.copy(), copy=True)
adata = adata[:,hvgs]
sc.pp.scale(adata) # This step will substantially extend the size of the adata
sc.pp.pca(adata)
sc.pl.pca_variance_ratio(adata, log=True,save='pca.pdf')
sc.pp.neighbors(adata, n_neighbors=50, n_pcs=30)
sc.tl.umap(adata)


# -------------- load results from cNMF analysis for downstream analysis ------------------ 
usage_norm, gep_scores, gep_tpm, topgenes = cnmf_obj.load_results(K=selected_K, density_threshold=density_threshold)
usage_norm.columns = ['Usage_%d' % i for i in usage_norm.columns]
usage_file = cnmf_obj.paths['consensus_usages__txt'] % (selected_K,'0_1')
gene_scores_file = cnmf_obj.paths['gene_spectra_score__txt'] % (selected_K,'0_1')
gene_tpm_file = cnmf_obj.paths['gene_spectra_tpm__txt'] % (selected_K,'0_1')

usageDf = pd.read_csv("CrossSpecial_cNMF/CrossSpecial_cNMF.usages.k_17.dt_0_1.consensus.txt",sep=r"\s+", index_col=0,header=0)
usageDf.columns = list(range(1, len(usageDf.columns) + 1))
usageDf["Max_Column"] = usageDf.idxmax(axis=1).astype(str)

adata.obs = pd.merge(left=adata.obs, right=usageDf, how='left', left_index=True, right_index=True)
sc.pl.umap(adata,color=['Max_Column'],frameon=False, legend_fontsize=15, save="_usage.png")

# UMAP visulization of a particular GEP 
adata.obs.columns = adata.obs.columns.map(str)
data_p = adata[adata.obs.species=="P",:]
data_h = adata[adata.obs.species=="H",:]
data_m = adata[adata.obs.species=="M",:]
sc.pl.umap(data_p,color=['5'],frameon=False, legend_fontsize=15, save="_pig_GEP5.png")
sc.pl.umap(data_h,color=['5'],frameon=False, legend_fontsize=15, save="_human_GEP5.png")
sc.pl.umap(data_m,color=['5'],frameon=False, legend_fontsize=15, save="_mice_GEP5.png")



# ---------------- get Top200 gene associated with each gene expression program ----------------（
df = pd.read_csv("CrossSpecial_cNMF/CrossSpecial_cNMF.spectra.k_17.dt_0_1.consensus.txt",sep="\t", index_col=0)
top_genes = []

for gep_name, row in df.iterrows():
    top200_genes = row.sort_values(ascending=False).head(200).index.tolist()
    top_genes.extend([(gep_name, gene) for gene in top200_genes])

output_df = pd.DataFrame(top_genes, columns=["GEP", "Gene"])
output_df.to_csv("top200_genes_per_GEP.tsv", sep="\t", index=False)


# load gene scores and identify the genes that are most associated with each program to see if we can interpret
top_genes = []
ngenes = 200
for gep in gep_scores.columns:
    top_genes.append(list(gep_scores.sort_values(by=gep,ascending=False).index[:ngenes]))

top_genes = pd.DataFrame(top_genes, index=gep_scores.columns).T
sc.pl.umap(adata, color=['CDC6','MKI67'], use_raw=True, ncols=4)
sc.pl.umap(adata, color=['MMD','CLEC1B'], use_raw=True, ncols=4)







































