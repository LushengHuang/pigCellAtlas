
rm(list = ls())
wd = "path"
setwd(wd)
library(patchwork)
library(ggplot2)
library(CellID)
library(tidyverse)
library(Seurat)
library(SeuratDisk)
library(cli)
library(crayon)
library(hdf5r)
library(Matrix)
library(R6)
library(rlang)
library(Seurat)
library(SeuratObject)
library(stringi)
library(withr)
library(tools)
library(scHCL)
library(scMCA)
library(data.table)

wd = "path"
setwd(wd)
##h5ad转seurat
filename <- 'endothelial'
Convert(paste0(filename,".h5ad"), dest = "h5seurat", overwrite = TRUE)
dat <- LoadH5Seurat(paste0(filename,".h5seurat")) 
dat.list <- SplitObject(dat, split.by = "expBatch")
dat.list <- lapply(X = dat.list, FUN = function(x){
#    x <- NormalizeData(x, verbose = FALSE)
    x <- FindVariableFeatures(x, verbose = FALSE)
})

##---------------3.2.1 Integreated Data -----------------------
options(future.globals.maxSize=1048576000000)
features <- SelectIntegrationFeatures(object.list = dat.list)
dat.anchors <- FindIntegrationAnchors(object.list = dat.list, anchor.features = features)
dat.combined <- IntegrateData(anchorset = dat.anchors)

#----------- specify that we will perform downstream analysis on the corrected data note that the
#------------ original unmodified data still resides in the 'RNA' assay
DefaultAssay(dat.combined) <- "integrated"

# -----------Run the standard workflow for visualization and clustering
dat.combined <- ScaleData(dat.combined, verbose = FALSE)
dat.combined <- RunPCA(dat.combined, npcs = 30, verbose = FALSE)
dat.combined <- RunUMAP(dat.combined, reduction = "pca", dims = 1:30)
dat.combined <- FindNeighbors(dat.combined, reduction = "pca", dims = 1:30)
dat0 <- dat

res <- c(0.6,0.8,1,1.2,1.4,1.6)
for(i in res){
    dat <- dat0
    dat.combined <- FindClusters(dat.combined,resolution=i)
    p1 <- DimPlot(dat.combined, group.by = "Individual_2")
    p3 <- DimPlot(dat.combined, reduction = "umap",label = TRUE)
    p5 <- DimPlot(dat.combined,group.by = "Organ")
    p7 <- DimPlot(dat.combined,group.by = "expBatch")

    ggsave(paste0(filename,i,"_Individual_2.png"), plot = p1, width = 10, height = 7)
    ggsave(paste0(filename,i,"_umap.png"), plot = p3,width = 10, height = 7)
    ggsave(paste0(filename,i,"_Organ.png"), plot = p5,width = 10, height = 7)
    ggsave(paste0(filename,i,"_expBatch.png"), plot = p7,width = 10, height = 7)

    saveRDS(dat.combined,file = paste0(filename,"_",i,"_mnn.rds"))
    write.csv(dat.combined@meta.data,file = paste0(filename,i,"_mnn_cellInfo.csv"))

    ###finding differentially expressed features
    dat.markers <- FindAllMarkers(dat.combined,only.pos = TRUE,min.pct = 0.25,logfc.threshold = 0.25)
    write.csv(dat.markers,file = paste0(filename,i,"_mnn_makergene.csv"))
    top10 <- dat.markers %>% group_by(cluster) %>% slice_max(n = 10, order_by = avg_log2FC)
    write.csv(top10,file = paste0(filename,i,"_mnn_makergene_top10.csv"))	
    all.difG <- dat.markers %>% select(gene,everything()) %>% subset(p_val < 0.05)
    write.csv(all.difG,paste0(filename,i,".all.difG.csv"),row.names = F,col.names = T)
			 
    ##top markergene可视化				  
    top5 <- dat.markers %>% group_by(cluster) %>% slice_max(n = 5, order_by = avg_log2FC)
    pHtmap <- DoHeatmap(dat.combined, features = top5$gene,group.by = "seurat_clusters", group.bar = T, size = 3,angle = 45)	
    ggsave(paste0(filename,i,"top5_htmap.pdf"),plot = pHtmap,width = 15,height = 10)
}



	



