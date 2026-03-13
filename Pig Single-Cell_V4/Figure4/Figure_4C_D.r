rm(list=ls());options(stringsAsFactors=FALSE)
library(tidyverse)
library(ggsci)
library(ggplot2) 
library(ggpattern)
library(ComplexHeatmap)
library(gplots)
library(reshape2)
library(gg.gap)

#------------------ Set parameters ------------------#
setwd('Pathway')
Plotdata = read.table(file = 'Celltype_tissue_specificTF_average_expression.csv', sep = ',',header = T)

#------------------ Preprocessing ------------------#
mat <- Plotdata[,-2]
rownames(mat)=Plotdata$org_celltype
normalize <- function(x) {
  return((x - min(x)) / (max(x) - min(x)))
}
normalized_mat <- apply(mat, 2, normalize)
mat=t(normalized_mat)
cellinfo <- Plotdata[,1:2]
cellinfo$celltype=factor(cellinfo$celltype,levels = unique(cellinfo$celltype))
colbroad<-colorRampPalette((pal_npg("nrc")(9)))(23)
names(colbroad) <- levels(cellinfo$celltype)
top_anno = HeatmapAnnotation(
  celltype = cellinfo$celltype,
  col = list(celltype = colbroad), show_legend = FALSE)
col_fun  <- colorpanel(100,low = "#F2F2F2", mid = "#fee1d4",high = "#e73827")

pdf(paste0('Tissue_celltype_specific_TF.pdf'), width=8, height=10)
Heatmap(mat,
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        show_column_names= F,
        show_row_names = T,
        row_names_gp = gpar(fontsize = 12),
        column_split = cellinfo$celltype,
        column_gap = unit(0.5, "mm"),
        column_title_gp = gpar(fontsize = 12),
        column_title_side = 'bottom',
        column_title_rot = 90,
        bottom_annotation = top_anno,
        heatmap_legend_param = list(
          title = "Scale Expression",
          title_position = "leftcenter-rot"),
        col = col_fun)
dev.off()

