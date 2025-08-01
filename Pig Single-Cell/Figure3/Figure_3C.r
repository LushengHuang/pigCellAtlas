
wd <- "path"
setwd(wd)

rm(list=ls());options(stringsAsFactors=FALSE)
library(openxlsx)
library(tidyverse)
library(circlize)	
library(RColorBrewer)
library(ggsci)
library(ComplexHeatmap)
library(plyr)

godsnot_102 <- c(
  "#FFFF00","#1CE6FF","#FF34FF","#FF4A46","#008941","#006FA6","#A30059","#FFDBE5","#7A4900","#0000A6",
  "#63FFAC","#B79762","#004D43","#8FB0FF","#997D87","#5A0007","#809693","#6A3A4C","#1B4400","#4FC601",
  "#3B5DFF","#4A3B53","#FF2F80","#61615A","#BA0900","#6B7900","#00C2A0","#FFAA92","#FF90C9","#B903AA",
  "#D16100","#DDEFFF","#000035","#7B4F4B","#A1C299","#300018","#0AA6D8","#013349","#00846F","#372101",
  "#FFB500","#C2FFED","#A079BF","#CC0744","#C0B9B2","#C2FF99","#001E09","#00489C","#6F0062","#0CBD66",
  "#EEC3FF","#456D75","#B77B68","#7A87A1","#788D66","#885578","#FAD09F","#FF8A9A","#D157A0","#BEC459",
  "#456648","#0086ED","#886F4C","#34362D","#B4A8BD","#00A6AA","#452C2C","#636375","#A3C8C9","#FF913F",
  "#938A81","#575329","#00FECF","#B05B6F","#8CD0FF","#3B9700","#04F757","#C8A1A1","#1E6E00","#7900D7",
  "#A77500","#6367A9","#A05837","#6B002C","#772600","#D790FF","#9B9700","#549E79","#FFF69F","#201625",
  "#72418F","#BC23FF","#99ADC0","#3A2465","#922329","#5B4534","#FDE8DC","#404E55","#0089A3","#CB7E98",
  "#A4E804","#324E72")
 
# KEGG immune pathways 
kegg1 <- read.table("Im_Kegg_int_exp_rmsn_sel_b1_Kegg0831.txt",sep = "\t",header = T,check.names = F)  
kegg1 <- kegg1[,-2]

# scavenger
kegg2 <- read.table("Im_Kegg_int_exp_rmsn_sel_b1_add_pt0901.txt",sep = "\t",header = T,check.names = F)
kegg2 <- kegg2[,-c(2:3,5:14)]
pathway <- "Scavenger"
names(kegg2)[2] <- "pathway"

subset1 <- kegg1[grep("SLA-DRA|SLA-DQB1|SLA-DRB1|IL6|CXCL12",kegg1[,1]),]
subset1 <- subset1[!duplicated(subset1[,"gene"]),]
subset1 <- subset1[-c(3,4),]

data <- rbind(kegg2,subset1)


# ---------------- classify and order the EC subtypes -------------------------------
  col_info <- as.data.frame(names(data)[-c(1:2)])
  names(col_info) <- "id"
  col_info$stage = as.data.frame(str_split_fixed(col_info$id,"[|]",n=3))[,1]
  col_info$Organ = as.data.frame(str_split_fixed(col_info$id,"[|]",n=3))[,2]
  col_info$celltype = as.data.frame(str_split_fixed(col_info$id,"[|]",n=3))[,3]
  
  ct10up <- read.csv("celltype10_update.csv",header = T)
  ct10_int = ct10up[ct10up$celltype %in% col_info$celltype ,]
  dim(ct10_int)
  col_info2 = merge(col_info,ct10_int,by = "celltype") #按照col_info的信息进行merge合并。
  col_info2 <- col_info2[col_info2$category != "others" ,] #过滤掉others
  length(table(ct10_int$celltype))
  
  col_info2 <- col_info2 %>% select(c("celltype","id","stage","Organ","category","Big_cate", "color")) #选择需要的列
  col_info2$stage_clt <- paste0(col_info2$stage,"-",col_info2$category)

  tissue_code <- read.xlsx("TableS1.xlsx",2)
  tissue_code <-tissue_code[,c(4:5,7)]
  #tissue_code$Organ<-gsub(" ",".",tissue_code$Organ)
  tissue_code <- distinct(tissue_code,Organ,.keep_all = T)
  
  col_info2 <- left_join(col_info2,tissue_code,by ="Organ")
  col_info2 <- arrange(col_info2,Organ.code)
  length(table(col_info2$celltype))
  head(col_info2,5)

  #
  od1 <- col_info2[grep("Fetus-artery",col_info2$stage_clt),]
  od2 <- col_info2[grep("Fetus-capillary",col_info2$stage_clt),]
  od4 <- col_info2[grep("Fetus-vein",col_info2$stage_clt),]
  od3 <- col_info2[grep("Fetus-lymphatic",col_info2$stage_clt),]
  od5 <- col_info2[grep("Sow-artery",col_info2$stage_clt),]
  od6 <- col_info2[grep("Sow-capillary",col_info2$stage_clt),]
  od8 <- col_info2[grep("Sow-vein",col_info2$stage_clt),]
  od7 <- col_info2[grep("Sow-lymphatic",col_info2$stage_clt),]

  all_od <- rbind(od1,od2,od4,od3,od5,od6,od8,od7)
  order_col <- as.character(all_od$id)
  data2 <- data %>% select(1,order_col)
  
  mat <- data2
  rownames(mat) <- mat$gene
  mat <- mat[,-1]
  mat <- as.matrix(mat) #25 393
 
  all_od$Organ <- factor(all_od$Organ, levels = unique(all_od$Organ))
  all_od$stage <- factor(all_od$stage, levels = c("Fetus","Sow"))
  all_od$System_abbr <- factor(all_od$System_abbr, unique(all_od$System_abbr))
  all_od$category <- factor(all_od$category, levels=c("artery","capillary","vein","lymphatic"))
  all_od$Big_cate <-  factor(all_od$Big_cate, levels=c("VEC","LEC"))

  panel_Organ <- godsnot_102[1:length(unique(all_od$Organ))]
  names(panel_Organ) <- unique(all_od$Organ)

  panel_stage <- c("Fetus"="#ea5455", "Sow" = "#2d4059")  
  panel_sys <- c("#1F77B4FF","#FF7F0EFF","#2CA02CFF","#D62728FF","#9467BDFF",
                 "#8C564BFF","#E377C2FF","#7F7F7FFF","#BCBD22FF","#17BECFFF",
                 "#AEC7E8FF","#FFBB78FF")
  sys <- c("NS","ES","CS","IS","RS","DS","FRS","MFI","MRS","US","MS","Others")
  names(panel_sys) <- sys
  panel_celltype <- c("artery"="#F5C014", "capillary" = "#4DBBD5FF","vein"="#FF7F0EFF",
                      "lymphatic" = "#2CA02CFF")
  panel_Big_cate <- c("VEC" = "#D62728FF","LEC" = "#BCBD22FF")

  col_ha <- HeatmapAnnotation(Stage = all_od$stage,
                              System = all_od$System_abbr,
                              col = list(Stage = panel_stage,
                                         System = panel_sys))
  mat <- mat[!(rownames(mat) %in% c("ENSSSCG00000035651-STAB2","ENSSSCG00000008769-CD68")),]
  rownames(mat)[rownames(mat) == "ENSSSCG00000000854-STAB2"] <- "STAB2"  
  max(mat)
  min(mat)
  
  
  col_fun  <- colorRamp2(c(0,ceiling(max(mat))/10,ceiling(max(mat))),c("#FFFFFF","#5178F1","#E0002C"))
  
  identical(as.character(colnames(mat)),as.character(all_od$id))  
  pdf(paste0("./kegg_scavenger.pdf"),width=12,height=4) #Scavenger
  rowsplit <- c(rep("Immune genes",5),rep("Scavenger genes",23))
  p <-Heatmap(mat,
              cluster_rows = FALSE,cluster_columns = FALSE,
              show_column_names= FALSE,show_row_names = TRUE,
              column_split = all_od$category,
			        row_split = rowsplit,
              top_annotation = col_ha,
              column_title = "",
              row_names_max_width = unit(14, "cm"), #设置行名的最大宽度
              row_names_side = "right",
              column_names_side = 'top',
              col = col_fun,
              rect_gp = gpar(col="white"),
              border = TRUE,
              column_order = colnames(mat),
			  row_names_gp = gpar(fontsize = 14),
			  column_names_gp = gpar(fontsize=14),
			  row_title_gp = gpar(fontsize=14),
			  column_title_gp = gpar(fontsize=14))
  draw(p)
  dev.off()

  