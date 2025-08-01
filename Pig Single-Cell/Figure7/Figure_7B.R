rm(list=ls())
library(Seurat)
library(RColorBrewer)
library(tidyverse)

human<-readRDS("A_Hum_merge_Pseu50_subanno.rds") #Hum
mouse<-readRDS("A_Mou_merge_Pseu50_subanno.rds") #Mou
Pig<-readRDS("A_Pig_merge_Pseu50_subanno.rds") #Pig

human
mouse
Pig

human <- GetAssayData(human, slot = "data")
mouse <- GetAssayData(mouse, slot = "data")
Pig <- GetAssayData(Pig, slot = "data")

human<-as.data.frame(human)
mouse<-as.data.frame(mouse)
Pig<-as.data.frame(Pig)


HCL_cell = NULL
HCL_phe = NULL

HCL_cell$Pseu_cell = colnames(human)
HCL_cell =as.data.frame(HCL_cell)
HCL_Sep = separate(HCL_cell, col = Pseu_cell, into = c("cell", "Pseu_num"), sep = "\\|Cell") 
HCL_Sep$Species = "HCL"
HCL_Sep$Pseu_cell = HCL_cell$Pseu_cell
HCL_phe = subset(HCL_Sep, select = c("Pseu_cell","Species","cell"))


MCA_cell = NULL
MCA_phe = NULL

MCA_cell$Pseu_cell = colnames(mouse)
MCA_cell =as.data.frame(MCA_cell)
MCA_Sep = separate(MCA_cell, col = Pseu_cell, into = c("cell", "Pseu_num"), sep = "\\|Cell") 
MCA_Sep$Species = "MCA"
MCA_Sep$Pseu_cell = MCA_cell$Pseu_cell
MCA_phe = subset(MCA_Sep, select = c("Pseu_cell","Species","cell"))


Pig_cell = NULL
Pig_phe = NULL

Pig_cell$Pseu_cell = colnames(Pig)
Pig_cell =as.data.frame(Pig_cell)
Pig_Sep = separate(Pig_cell, col = Pseu_cell, into = c("cell", "Pseu_num"), sep = "\\|Cell") 
Pig_Sep$Species = "Pig"
Pig_Sep$Pseu_cell = Pig_cell$Pseu_cell
Pig_phe = subset(Pig_Sep, select = c("Pseu_cell","Species","cell"))


P1<-rbind(Pig_phe,HCL_phe)
P1<-rbind(P1, MCA_phe)
colnames(P1)<-c("Sample_ID","Study_ID","Celltype") 

data<-cbind(Pig,human)
data<-cbind(data, mouse)
data[is.na(data)]<-0
data1<-data[,as.character(P1$Sample_ID)]


source("2017-08-28-runMN-US.R")
library(gplots)
library(RColorBrewer)

celltypes1 <-unique(as.character(P1$Celltype))

var.genes1=get_variable_genes(data1,P1)
length(var.genes1)

write.table(var.genes1,"var.genes_75_PHM.out",sep="\t",quote=F)#####--------
celltype.NV=run_MetaNeighbor_US(var.genes1,data1,celltypes1,P1)
write.table(celltype.NV,file="celltype.NV_SRS_75_Insec_Adu_PHM.out",sep="\t",quote=F)

cols=rev(colorRampPalette(brewer.pal(11,"RdYlBu"))(100))
breaks=seq(0,1,length=101)
pdf("celltype.NV_SRS_75_Adu_Insec1119_PHM.pdf")  #########--------------------------
heatmap.2(celltype.NV,trace="none",density.info="none",col=cols,breaks=breaks,cexRow=0.3,cexCol=0.3)
dev.off()
top_hits=get_top_hits(celltype.NV,P1,threshold=0.9,filename="top_hits_SRS_75_Insec_PHM.out") 
top_hits=get_top_hits(celltype.NV,P1,threshold=0.8,filename="top_hits_SRS_0.8_75_Insec_PHM.out")
top_hits=get_top_hits(celltype.NV,P1,threshold=0.7,filename="top_hits_SRS_0.7_75_Insec_PHM.out")
top_hits=get_top_hits(celltype.NV,P1,threshold=0.6,filename="top_hits_SRS_0.6_75_Insec_PHM.out")


#Plot Heatmap
#################################################################
rm(list=ls())
library("data.table") 
library("R.utils")
library("ggsci")
library("ComplexHeatmap")
library("viridis")
library("circlize")
library("RColorBrewer")

mat <- fread("celltype.NV_SRS_75_Insec_Adu_PHM.out",data.table=F)
rownames(mat)=mat[,1]
mat=mat[,-1]
mat=mat[,-9]
mat=mat[-9,]
mat=mat[-41,]
mat=mat[,-41]
mat=mat[-71,]
mat=mat[,-71]
info <- read.csv("Adu_ct3.csv", head = F)
names(info) <- c('species','celltype','module','ctmod5','ctnum',"categ")
info$module=factor(info$module)
info$ctmod5=factor(info$ctmod5, levels = c("Other_Cells","Brain_Cells", "Stromal_Cells", "Epithelial_Cells", "Immune_Cells"))
info$ctnum=factor(info$ctnum, levels = c("7","6","5","4","3","2","1"))
info$categ=factor(info$categ, levels = c("Other_Cells",
                                         "Brain_Cells-Astrocytes",
                                         "Brain_Cells-Granule_Cells",
                                         "Brain_Cells-Oligodendrocytes",
                                         "Stromal_Cells-Endothelial_Cells",
                                         "Stromal_Cells-Fibroblasts",
                                         "Stromal_Cells-Smooth_Muscle_Cells",
                                         "Epithelial_Cells-Epithelial_Cells",
                                         "Immune_Cells-Neutrophils",
                                         "Immune_Cells-B_Cells",
                                         "Immune_Cells-T_Cells",
                                         "Immune_Cells-Macrophages",
                                         "Immune_Cells-Dendritic_Cells"
))
rownames(info) <- info$celltype
##############人########鼠########猪###
colspe<-c('#E64B35FF','#4DBBD5FF','#00A087FF')
colmod<-c(
"#8F7700FF",
"#4DBBD5FF",
"#00A087FF",
"#3C5488FF",
"#F39B7FFF",
"#8491B4FF",
"#91D1C2FF",
"#DC0000FF",
"#7E6148FF",
"#B09C85FF",
"#3B4992FF",
"#FF0000FF",
"#008B45FF",
"#631879FF",
"#008280FF",
"#A20056FF",
"#EEA236FF"
)

mod5 <- c(
"#8F7700FF",
"#4DBBD5FF",
"#00A087FF",
"#3C5488FF",
"#F39B7FFF"
)

ctnumber<-c(
  "#8F7700FF",
  "#4DBBD5FF",
  "#00A087FF",
  "#3C5488FF",
  "#F39B7FFF",
  "#8491B4FF",
  "#91D1C2FF"
)

names(colspe) <- levels(info$specie)
names(colmod) <- levels(info$module)
names(mod5) <- levels(info$ctmod5)
names(ctnumber) <- levels(info$ctnum)

col_fun <- colorRamp2(
  c(0, 0.5, 1), 
  c("#e6e6e6", "#e6e6e6" ,'#ff0016'))#转置

order = c(
  "P-Hepatocytes",
  "H-Hepatocytes",
  "M-Hepatocytes",
  "P-Adipocytes",
  "H-Adipocytes",
  "M-Pancreas_exocrine_cells",
  "P-Astrocytes",
  "H-Astrocytes",
  "M-Astrocytes",
  "P-Neutrophils",
  "H-Neutrophils",
  "M-Neutrophils",
  "P-Proliferating_B_cells",
  "H-Proliferating_B_cells",
  "M-Proliferating_B_cells",
  "P-B_cells",
  "H-B_cells",
  "M-B_cells",
  "P-Plasma_cells",
  "H-Plasma_cells",
  "M-Plasma_cells",
  "P-Skeletal_muscle_cells",
  "H-Skeletal_muscle_cells",
  "M-Fibroblasts",
  "P-Smooth_muscle_cells",
  "H-Smooth_muscle_cells",
  "P-Proliferating_T_cells",
  "H-Proliferating_T_cells",
  "M-Proliferating_T_cells",
  "P-T_cells",
  "H-T_cells",
  "M-Skeletal_muscle_cells",
  "P-Oligodendrocytes",
  "H-Oligodendrocytes",
  "M-Oligodendrocytes",
  "M-Smooth_muscle_cells",
  "M-T_cells",
  "P-Dendritic_cells",
  "H-Dendritic_cells",
  "M-Dendritic_cells",
  "P-VEC",
  "H-VEC",
  "M-VEC",
  "P-Endothelial_cells",
  "H-Endothelial_cells",
  "M-Endothelial_cells",
  "P-Adipose_fibroblasts",
  "H-Adipose_fibroblasts",
  "M-Adipose_fibroblasts",
  "P-Adrenal_Gland_epithelial_cells",
  "H-Adrenal_Gland_epithelial_cells",
  "M-Adrenal_Gland_epithelial_cells",
  "P-Alveolar_type_II_cells",
  "H-Alveolar_type_II_cells",
  "M-Alveolar_type_II_cells",
  "P-Pancreas_exocrine_cells",
  "H-Pancreas_exocrine_cells",
  "P-Proximal_tubule_cells",
  "H-Proximal_tubule_cells",
  "M-Proximal_tubule_cells",
  "P-Stem_cells",
  "H-Stem_cells",
  "M-Stem_cells",
  "P-Epithelial_cells",
  "H-Epithelial_cells",
  "M-Epithelial_cells",
  "P-Gastric_chief_cells",
  "H-Gastric_chief_cells",
  "M-Gastric_chief_cells",
  "P-Stomach_epithelial_cells",
  "H-Stomach_epithelial_cells",
  "M-Stomach_epithelial_cells",
  "P-Interstinal_epithelial_cells",
  "H-Interstinal_epithelial_cells",
  "M-Interstinal_epithelial_cells",
  "P-Enterocytes",
  "H-Enterocytes",
  "M-Enterocytes",
  "P-Fibroblasts",
  "H-Fibroblasts",
  "P-Granule_cells",
  "H-Granule_cells",
  "M-Granule_cells",
  "P-Adipose_macrophages",
  "H-Adipose_macrophages",
  "P-Lung_macrophages",
  "M-Lung_macrophages",
  "H-Lung_macrophages",
  "M-Liver_macrophages",
  "P-Liver_macrophages",
  "H-Liver_macrophages",
  "P-Macrophages",
  "H-Macrophages",
  "M-Macrophages")
  
mat = mat[, order]
mat = mat[order,]
head(mat,4)
tail(mat, 4)

pdf("Adu_PHM_heatmap20230101_13.pdf", width=22, height=20)

Heatmap(as.matrix(mat),
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        col = col_fun,
        #cluster_columns = F,
        show_column_names= TRUE,
        show_row_names = FALSE,
        show_heatmap_legend = FALSE,
        #top_annotation = top_anno,
        column_title = NULL,
        column_names_max_height = unit(14, "cm"),
        row_split = info$categ,
        column_split = info$categ,
        row_gap = unit(c(4,2,2,4,2,2,4,4,2,2,2,2,2), "mm"),
        column_gap=unit(c(4,2,2,4,2,2,4,4,2,2,2,2,2),"mm"),##隔断间隙
        heatmap_legend_param = list(title_gp = gpar(fontsize = 8, fontface = "bold"))
)+Heatmap(info$species, name = "species", col = colspe, width = unit(5, "mm"),
          heatmap_legend_param = list(labels_gp = gpar(fontsize = 14),title_gp = gpar(fontsize = 16, fontface = "bold")))+
  Heatmap(info$module, name = "Main Celltype", col = colmod, width = unit(5, "mm"),
          heatmap_legend_param = list(labels_gp = gpar(fontsize = 14),title_gp = gpar(fontsize = 16, fontface = "bold")))

dev.off()
