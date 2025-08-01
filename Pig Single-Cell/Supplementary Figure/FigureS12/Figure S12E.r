rm(list=ls())
library(stringr)
library(reshape2)
library(ggplot2) 
library(ggsci)
library(gridExtra)
library(cowplot)
library(tidyverse)
library(circlize)
library(grid)
library(ComplexHeatmap)
library(dplyr)
library(colorspace)

#------------------ Set parameters ------------------#
setwd('Pathway')
Plotdata = read.table(file = 'Cel-specific transcription factors conserved across species.csv', sep = ',',header = T)

#------------------ Preprocessing ------------------#
TF=unique(Plotdata$TF)
data1=read.table('Mouse_A_celltype_TF_expression.csv', header = T, stringsAsFactors = F, sep = ',')
mat <- data1
mat=mat[,-1]
mat=scale(mat)
data1[,-1]=mat
names(data1) <- gsub("\\.", "-", names(data1))
data1=data1[,TF]
data1$species="M"
data2=read.table('Human_A_celltype_TF_expression.csv', header = T, stringsAsFactors = F, sep = ',')
mat <- data2
mat=mat[,-1]
mat=scale(mat)
data2[,-1]=mat
names(data2) <- gsub("\\.", "-", names(data2))
data2=data2[,TF]
data2$species="H"
data3=read.table('Pig_A_celltype_TF_expression.csv', header = T, stringsAsFactors = F, sep = ',')
mat <- data3
mat=mat[,-1]
mat=scale(mat)
data3[,-1]=mat
names(data3) <- gsub("\\.", "-", names(data3))
data3=data3[,TF]
data3$species="P"
data=rbind(data1,data2,data3)

celltype=unique(Plotdata$celltype)

data$celltype=factor(data$celltype, levels = celltype, ordered = T) 
data$species=factor(data$species, levels = c('H','M','P'), ordered = T)
data=data[order(data$celltype),]

col_fun <- colorRamp2(
  c(-2, 2, 4, 6), 
  c("#F7F7F7", "#F1F1F1" ,"#fb694a", "#67000d")
)

colspe<-pal_npg("nrc")(3)
names(colspe) <- levels(data$species)

colmod <- qualitative_hcl(28, palette = "Dynamic")
names(colmod) <- levels(data$celltype)

numeric_cols <- sapply(data[,2:200], is.numeric)
data_numeric <- data[,2:200][, numeric_cols]
data_matrix <- as.matrix(data_numeric)

pdf("Adult_Cell_type_specific_TFs_across_species_heatmap.pdf", width=28, height=15)

Heatmap(data$celltype, name = "Celltype", col = colmod, width = unit(5, "mm"), 
        row_split = data$celltype, row_gap = unit(1.2, "mm"), 
        row_title_rot = 0, row_title_gp = gpar(fontsize = 14),
        heatmap_legend_param = list(labels_gp = gpar(fontsize = 10), 
                                    title_gp = gpar(fontsize = 12, fontface = "bold"))) +
  Heatmap(data_matrix, 
          cluster_rows = FALSE, 
          cluster_columns = FALSE, 
          column_names_gp = gpar(fontsize = 8), 
          show_row_names = FALSE, 
          col = col_fun, 
          heatmap_legend_param = list(labels_gp = gpar(fontsize = 10), 
                                      title_gp = gpar(fontsize = 12, fontface = "bold"), 
                                      title = "Scale Expression")) +
  Heatmap(data$species, name = "Species", col = colspe, width = unit(5, "mm"), 
          heatmap_legend_param = list(labels_gp = gpar(fontsize = 10), 
                                      title_gp = gpar(fontsize = 12, fontface = "bold")))
dev.off()


