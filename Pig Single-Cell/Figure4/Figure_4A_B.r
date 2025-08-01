rm(list=ls())
library(stringr)
library(reshape2)
library(ggplot2) 
library(ggsci)
library(gridExtra)
library(cowplot)
library(circlize)
library(grid)
library(ComplexHeatmap)
library(dplyr)

#------------------ Set parameters ------------------#
setwd('Pathway')
Plotdata = read.table(file = 'TF_average_expression_in_fatus.csv', sep = ',',header = T)
specificTF <- read.csv("E_Celltype_specific_TF_tau.csv",sep = ",")

#------------------ Preprocessing ------------------#
expMat <- as.matrix(Plotdata)
names(specificTF)[1]='NO'
head(specificTF)
specificTF$TFName=gsub('-','.',specificTF$TFName)
mat <- expMat[,specificTF$TFName]
mat=t(scale(mat))
cellinfo <- as.data.frame(colnames(mat))
names(cellinfo) <- c('celltype')
cellinfo=left_join(cellinfo,celltype_info[,2:4],by='celltype')
cellinfo=cellinfo[order(cellinfo$broadcelltype_order,cellinfo$celltype),]
cellinfo$celltype=factor(cellinfo$celltype,levels = unique(cellinfo$celltype))
cellinfo$broadcelltype=factor(cellinfo$broadcelltype,levels = unique(cellinfo$broadcelltype))
rownames(cellinfo)=NULL
head(cellinfo)
colbroad<-c("#FF7F0EFF","#BCBD22FF","#1F77B4FF","#4E91C5","#7EACD6","#AEC7E8FF",
            "#17BECFFF","#9EDAE5FF","#D62728FF","#E34C4C",
            "#F17271","#FF9896FF","#E377C2FF","#F7B6D2FF","#DBDB8DFF","#9467BDFF",
            "#AC8BC9","#C5B0D5FF","#8C564BFF","#A8796F","#C49C94FF","#2CA02CFF",
            "#61BF5A","#98DF8AFF","#7F7F7FFF","#C7C7C7FF")
names(colbroad) <- levels(cellinfo$broadcelltype)

top_anno = HeatmapAnnotation(
        celltype = cellinfo$broadcelltype,
        col = list(celltype = colbroad))
mat=as.data.frame(mat)
newcol=cellinfo$celltype
mat = mat %>% select(newcol)

###########Select TF to be displayed
showTF=c("SOX17","ERG","SOX18","BCL6B","SOX7","GATA2","STAT4","FOXS1","TEAD3",
         "NR2F2","TBX18","SOX10","MEF2A","MAF","IRF5","IRF8","SPI1","TFEC","CEBPA",
         "BHLHE41","NR2F1","NFATC1","FOXA1","SOX9")

gene_pos <- match(showTF,rownames(mat))
row_anno <- rowAnnotation(gene = anno_mark(at = gene_pos, labels = gene))

#color
col = colorpanel(100,low = "#fcfbe6", mid = "#fffcc0",high = "#003973")
col_fun <- colorRamp2(
        c(-2, 0, 2), 
        c("#8c510a", "white", "#01665e")
)

#
pdf(paste0('Fatus_celltype_specific_TF_average_expression.pdf'), width=9, height=9)
Heatmap(mat,
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        show_column_names= FALSE,
        show_row_names = FALSE,
        column_split = cellinfo$broadcelltype,
        column_gap = unit(0.5, "mm"),
        top_annotation = top_anno,
        column_title = NULL,
        right_annotation = row_anno,
        heatmap_legend_param = list(
                title = "Scale Expression",
                title_position = "leftcenter-rot"),
        col = col_fun)
dev.off()

###############Adult_celltype_specific_TF_average_expression
Plotdata = read.table(file = 'TF_average_expression_in_adult.csv', sep = ',', row.names = 1,header = T)
expMat <- as.matrix(Plotdata)
specificTF <- read.csv("Adult_celltype_specific_TF_tau.csv",sep = ",")

names(specificTF)[1]='NO'
specificTF$TFName=gsub('-','.',specificTF$TFName)
mat <- expMat[,specificTF$TFName]
mat=t(scale(mat))
cellinfo <- as.data.frame(colnames(mat))
names(cellinfo) <- c('celltype')
cellinfo=left_join(cellinfo,celltype_info[,2:4],by='celltype')
cellinfo=cellinfo[order(cellinfo$broadcelltype_order,cellinfo$celltype),]
cellinfo$celltype=factor(cellinfo$celltype,levels = unique(cellinfo$celltype))
cellinfo$broadcelltype=factor(cellinfo$broadcelltype,levels = unique(cellinfo$broadcelltype))
rownames(cellinfo)=NULL
head(cellinfo)
colbroad<-c("#FF7F0EFF","#BCBD22FF","#1F77B4FF","#4E91C5","#7EACD6","#AEC7E8FF",
            "#17BECFFF","#9EDAE5FF","#D62728FF","#E34C4C",
            "#F17271","#FF9896FF","#E377C2FF","#F7B6D2FF","#DBDB8DFF","#9467BDFF",
            "#AC8BC9","#C5B0D5FF","#8C564BFF","#A8796F","#C49C94FF","#2CA02CFF",
            "#61BF5A","#98DF8AFF","#C7C7C7FF")
names(colbroad) <- levels(cellinfo$broadcelltype)

top_anno = HeatmapAnnotation(
        celltype = cellinfo$broadcelltype,
        col = list(celltype = colbroad))
mat=as.data.frame(mat)
newcol=cellinfo$celltype
mat = mat %>% select(newcol)

###########Select TF to be displayed
showTF=c("SOX17","ERG","SOX18","BCL6B","SOX7","EDN1","KLF2","POU2F3","RUNX3",
         "TCF7","TEAD3","FOXS1","JUNB","TFDP2","SOX10","MEF2A","IRF5","IRF8",
         "TFEC","SPIC","SMAD1","NFATC1","FOXA3","MYB","PAX5","RORB","LHX2","NR2E1","THRB")

gene_pos <- match(showTF,rownames(mat))
row_anno <- rowAnnotation(gene = anno_mark(at = gene_pos, labels = gene))

col_fun <- colorRamp2(
        c(-2, 0, 2), 
        c("#8c510a", "white", "#01665e")
)

pdf(paste0('Adult_celltype_specific_TF_average_expression0608.pdf'), width=9, height=9)
Heatmap(mat,
        cluster_rows = FALSE,
        #show_row_dend = F,
        cluster_columns = FALSE,
        show_column_names= FALSE,
        show_row_names = FALSE,
        #col = col_fun,
        column_split = cellinfo$broadcelltype,
        column_gap = unit(0.5, "mm"),
        top_annotation = top_anno,
        column_title = NULL,
        right_annotation = row_anno,
        heatmap_legend_param = list(
                title = "Scale Expression",
                title_position = "leftcenter-rot"),
        col = col_fun)
dev.off()
