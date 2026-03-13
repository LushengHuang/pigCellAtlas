############################################
#1.Extract 500 cells from each sample
import scanpy as sc
import numpy as np
import anndata as ad


adata = sc.read("merge_405Samples_bbknn_batch_platform_res3_0807_all_correctV5_20251023.h5ad")
samples = adata.obs["SampleID"].unique()
sample_sizes = [500]

for n in sample_sizes:
    idx_list = []
    for s in samples:
        s_idx = np.where(adata.obs["SampleID"] == s)[0]
        
        if len(s_idx) <= n:
            idx_list.extend(s_idx)
        else:
            idx_list.extend(np.random.choice(s_idx, n, replace=False))
    
    sampled_adata = adata[idx_list, :].copy()
    sampled_adata
    sampled_adata.obs["SampleID"].value_counts()
    sampled_adata.write(f"merge_405Samples_Sample_New_{n}.h5ad")


#############################################
#2. Integrate using Harmony
import scanpy as sc
import numpy as np
import pandas as pd

SamSize = 500
bat = "Pt"
Met = "harmony"

dat = sc.read("merge_405Samples_Sample_New_500.h5ad")

sc.pp.highly_variable_genes(dat, n_top_genes = 2000) #this function expects logarithmized data
sc.tl.pca(dat,svd_solver='arpack')

sc.external.pp.harmony_integrate(dat, 'Platform')
sc.pp.neighbors(dat,use_rep='X_pca_harmony')
sc.tl.umap(dat)

sc.settings.verbosity = 3
sc.logging.print_header()
sc.settings.set_figure_params(dpi=3600,facecolor='white')

sc.pl.umap(dat,color=['CT76'],legend_loc="right margin",frameon=False,legend_fontsize=15,palette = color_181,size = 2,title = " ",save='_CT76_4_20251027_Dpi3600.pdf')
sc.pl.umap(dat,color=['CT76'],legend_loc="on data",frameon=False,legend_fontsize=3,palette = color_181,size = 2,title = " ",save='_CT76_4_ondata_20251027_Dpi3600.pdf')


sc.tl.rank_genes_groups(dat,'CT76',method='wilcoxon')
result = dat.uns['rank_genes_groups']
groups = result['names'].dtype.names

markerGeneRank_all = pd.DataFrame({group + '_' + key[:1]: result[key][group] for group in groups for key in ['names', 'pvals']})
markerGeneRank_all.to_csv('AtlasSampling500_markerGeneInfo_all1024.csv')

markerGeneRank500 = pd.DataFrame({group + '_' + key[:1]: result[key][group] for group in groups for key in ['names','pvals']}).head(500)
#markerGeneRank500.to_excel(sample+'_markerGeneInfo_500.xls')
markerGeneRank500.to_csv('markerGeneInfo_500_1024.csv')

dat.write('merge_405Samples_Sample_New_500_harmonyPlatForm2000_diff_marker_5.h5ad',compression='gzip')



color_181 = [
"#FFCCCCFF",
"#8A4198FF",
"#91D1C2FF",
"#0099B4FF",
"#CD534CFF",
"#7AA6DCFF",
"#003399FF",
"#357EBDFF",
"#FFD147FF",
"#B24745FF",
"#802268FF",
"#42B540FF",
"#E18727FF",
"#155F83FF",
"#CDDEB7FF",
"#58593FFF",
"#17BECFFF",
"#6A6599FF",
"#7E6148FF",
"#99CCFFFF",
"#9900CCFF",
"#D62728FF",
"#8C564BFF",
"#0073C2FF",
"#A20056FF",
"#666600FF",
"#79CC3DFF",
"#46B8DAFF",
"#E377C2FF",
"#FFA319FF",
"#4A6990FF",
"#99991EFF",
"#BB0021FF",
"#5DB1DDFF",
"#2CA02CFF",
"#9632B8FF",
"#8A9045FF",
"#FF1463FF",
"#CC0C00FF",
"#FF410DFF",
"#80796BFF",
"#F39B7FFF",
"#5A655EFF",
"#CCFFFFFF",
"#749B58FF",
"#E64B35FF",
"#D58F5CFF",
"#D595A7FF",
"#358000FF",
"#CE3D32FF",
"#FFDC91FF",
"#1F77B4FF",
"#BC3C29FF",
"#3B4992FF",
"#008EA0FF",
"#00FFFFFF",
"#990033FF",
"#008B45FF",
"#996600FF",
"#AD002AFF",
"#6EE2FFFF",
"#99CC00FF",
"#95CC5EFF",
"#C16622FF",
"#9467BDFF",
"#003C67FF",
"#D0DFE6FF",
"#8F3931FF",
"#350E20FF",
"#6699FFFF",
"#00CC99FF",
"#339900FF",
"#374E55FF",
"#0099CCFF",
"#FDAF91FF",
"#B09C85FF"
]
