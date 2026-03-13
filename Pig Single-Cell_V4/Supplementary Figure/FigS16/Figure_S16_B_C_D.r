library(Seurat)
library(tidyverse)
library(cowplot)
library(patchwork)
# co-expression network analysis packages:
library(WGCNA)
library(hdWGCNA)
# using the cowplot theme for ggplot
theme_set(theme_cowplot())
# set random seed for reproducibility
set.seed(12345)

#optionally enable multithreading
#enableWGCNAThreads(nThreads = 10)


#------------------ Load data ------------------#
setwd('###')
seurat_obj <- readRDS('17sam_113072_batchI_bbknn_ct2_v1.rds')
seurat_obj <-seurat_obj[,seurat_obj@meta.data$ct1=='CM']

#------------------ Set up Seurat object for WGCNA ------------------#
seurat_obj <- SetupForWGCNA(
  seurat_obj,
  gene_select = "fraction", # the gene selection approach
  fraction = 0.05, # fraction of cells that a gene needs to be expressed in order to be included
  wgcna_name = "CM" # the name of the hdWGCNA experiment
)


#------------------  construct metacells  in each group ------------------#
seurat_obj <- MetacellsByGroups(
  seurat_obj = seurat_obj,
  group.by = c("ct1","S_I2","Individual_2","ShortName","SampleID",'Group'),
  reduction = 'pca', 
  k = 25, # nearest-neighbors parameter
  max_shared = 10, # maximum number of shared cells between two metacells
  ident.group = 'S_I2'
)
#------------------  normalize metacell expression matrix ------------------#
seurat_obj <- NormalizeMetacells(seurat_obj)


#------------------  Extract cell types and construct the expression matrix ------------------#
seurat_obj <- SetDatExpr(
  seurat_obj,# the name of the group of interest in the group.by column
  group_name='CM',
  group.by='ct1', # the metadata column containing the cell type info. This same column should have also been used in MetacellsByGroups
  assay = 'RNA', # using RNA assay
  slot = 'data' # using normalized data
)
#------------------  Select soft-thresholding power ------------------#
seurat_obj <- TestSoftPowers(
  seurat_obj,
  networkType = 'signed' # you can also use "unsigned" or "signed hybrid"
)
plot_list <- PlotSoftPowers(seurat_obj)

filename<-'CM'
#------------------   assemble with patchwork ------------------#
png(paste0(filename,"_softpower.png"),width=2800,height=2800,res=400)
wrap_plots(plot_list, ncol=2)
dev.off()
power_table <- GetPowerTable(seurat_obj)
head(power_table)

#------------------  Construct co-expression network ------------------#
seurat_obj <- ConstructNetwork(
  seurat_obj, soft_power=5,
  setDatExpr=FALSE,
  tom_name = 'CM' # name of the topoligical overlap matrix written to disk
)

png(paste0(filename,"_PlotDendrogram.png"),width=2800,height=1500,res=600)
PlotDendrogram(seurat_obj, main='CM hdWGCNA Dendrogram')
dev.off()
#------------------  Optional: Check the Topological Overlap Matrix (TOM) ------------------#
TOM <- GetTOM(seurat_obj)

#------------------  Calculate consensus module eigengenes ------------------#
seurat_obj <- ScaleData(seurat_obj, features=VariableFeatures(seurat_obj))
seurat_obj <- ModuleEigengenes(
 seurat_obj,
 group.by.vars="S_I2"
)

hMEs <- GetMEs(seurat_obj)

#------------------   module eigengenes  ------------------#
MEs <- GetMEs(seurat_obj, harmonized=FALSE)

#------------------  Calculate module connectivity ------------------#

seurat_obj <- ModuleConnectivity(
  seurat_obj,
  group.by = "S_I2", group_name = c('P1','P2','A1','R1','R2')
 )

seurat_obj <- ResetModuleNames(
  seurat_obj,
  new_name = c('CM-M')
)

#------------------  Obtain the module assignment table ------------------#
modules <- GetModules(seurat_obj)
write.csv(modules,paste0(filename,"_all_modules_gene.csv"))
hub_df <- GetHubGenes(seurat_obj, n_hubs = 10)



#------------------  Calculate the eigengene score ------------------#

seurat_obj <- ModuleExprScore(
  seurat_obj,
  n_genes = 25,
  method='Seurat'
)




MEs <- GetMEs(seurat_obj, harmonized=TRUE)
mods <- colnames(MEs); mods <- mods[mods != 'grey']
seurat_obj@meta.data <- cbind(seurat_obj@meta.data, MEs)

#------------------Module Eigengene Score visualization ------------------#
p <- VlnPlot(
    seurat_obj,
    features = mods,
    group.by = 'ShortName',
    split.by = 'I2',
	split.plot=TRUE,
	pt.size = 0
   )
     
pdf(paste0(filename, '_hME_vln-site_I2.pdf'), width=20, height=10)
p
dev.off()

#------------------ Save processed data ------------------#
saveRDS(seurat_obj, file=paste0(filename,'_hdWGCNA_object_I2.rds'))