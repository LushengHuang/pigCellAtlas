
rm(list=ls());options(stringsAsFactors=FALSE)

wd <- "~\\04_2022_scATLAS\\202409_immuneSubtypeHeatmap"
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
kegg1 <- read.table("Im_Kegg_int_exp_rmsn_sel_b1_Kegg0831.txt",sep = "\t",header=T,check.names=F)  
names(kegg1) <- gsub("Accessory genitalg land","Accessory genital gland",names(kegg1))
kegg1 <- kegg1[,-2]

idx <- grep("^IL",kegg1$gene)
ILgenes <- unique(kegg1[idx,"gene"])
ILgenes <- ILgenes[-grep("R|S",ILgenes)]

# scavenger
kegg2 <- read.table("Im_Kegg_int_exp_rmsn_sel_b1_add_pt0901.txt",sep = "\t",header = T,check.names = F)
names(kegg2) <- gsub("Accessory genitalg land","Accessory genital gland",names(kegg2))
kegg2 <- kegg2[,-c(2:3,5:14)]
pathway <- "Scavenger"
names(kegg2)[2] <- "pathway"

subset1 <- kegg1[grep("SLA-DRA|SLA-DQB1|SLA-DRB1|IL6|CXCL12",kegg1[,1]),]
subset1 <- subset1[!duplicated(subset1[,"gene"]),]
subset1 <- subset1[-c(3,4),]

subset2 <- kegg1[kegg1$gene %in% ILgenes,]
subset2 <- subset2[!duplicated(subset2[,"gene"]),]

data <- rbind(kegg2,subset2)


  #对列进行注释的准备
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

  #将选定的匹配到的子集行分别选出
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
  
  mat <- data2[data2$gene %in% ILgenes,]
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
  
#  col_ha = HeatmapAnnotation(BroadcellType = all_od$Big_cate,
#                             Celltype = all_od$category,
#                             Stage = all_od$stage,
#                             System = all_od$System_abbr,
#                             Organ = all_od$Organ,
#                             col = list(BroadcellType = panel_Big_cate,
#                                        Celltype= panel_celltype,
#                                        Stage = panel_stage,
#                                        System = panel_sys,
#                                        Organ = panel_Organ))


  col_ha <- HeatmapAnnotation(Stage = all_od$stage,
                              System = all_od$System_abbr,
                              col = list(Stage = panel_stage,
                                         System = panel_sys))
  mat <- mat[!(rownames(mat) %in% c("ENSSSCG00000035651-STAB2","ENSSSCG00000008769-CD68")),]
  rownames(mat)[rownames(mat) == "ENSSSCG00000000854-STAB2"] <- "STAB2"  
  max(mat)
  min(mat)
  
  #col_fun <- structure(c("white","red"), names = c("0", "1")) # black, red, green, blue
  col_fun  <- colorRamp2(c(0,ceiling(max(mat))/8,ceiling(max(mat))),c("#FFFFFF","#5178F1","#E0002C"))
  
  identical(as.character(colnames(mat)),as.character(all_od$id))  
  pdf(paste0("./ILgenes.pdf"),width=12,height=7) #Scavenger
  p <-Heatmap(mat,
              cluster_rows = FALSE,cluster_columns = FALSE,
              show_column_names= FALSE,show_row_names = TRUE,
              column_split = all_od$category,
              top_annotation = col_ha,
              column_title = "",
              row_names_max_width = unit(14, "cm"), #设置行名的最大宽度
              row_names_side = "right",
              column_names_side = 'top',
              col = col_fun,
              rect_gp = gpar(col="white"),
              border = TRUE,
              column_order = colnames(mat),
			  row_names_gp = gpar(fontsize = 15),
			  column_names_gp = gpar(fontsize=15),
			  row_title_gp = gpar(fontsize=15),
			  column_title_gp = gpar(fontsize=15))
  draw(p)
  dev.off()



# Draw the violin plot for different genes

plotdat <- data.frame(col_info2,t(mat2),group1 = "")
plotdat[plotdat$System_abbr=="NS","group1"] <- "Brain"
plotdat[plotdat$System_abbr!="NS","group1"] <- "Peripheral"

plotdat$group2 <- paste(plotdat$group1,plotdat$stage,sep="-")
plotdat$category <- as.factor(plotdat$category)

library(ggplot2)
library(cowplot)

pdf("IL7.pdf",width=16,height=3)
ggplot(plotdat, aes(x=category,y=IL7)) + 
  geom_violin(aes(fill=category)) + 
  geom_jitter(aes(fill=category),width=0.1) +
  theme_cowplot() +
  labs(y="Normalized expression-IL7") +
  labs(x="Vessel types") +
  guides(fill=FALSE) + facet_grid(.~group2)
dev.off()






















  mat2 <- mat[1:23,]
  col_fun  <- colorRamp2(c(0,ceiling(max(mat2))/20,ceiling(max(mat2))),c("#FFFFFF","#5178F1","#E0002C"))
  
  pdf(paste0("./kegg_scavenger_only.pdf"),width=20,height=12) #Scavenger
  geneOrder <- rev(c('CD36',"LY75",'SCARA3','SCARB2',"STAB2","MRC1","STAB1",'SCARB1',"OLR1",
                     "MARCO","CD209","COLEC12","SCARF1","SCARF2","CD14","CD163","CD68","SSC5D","MSR1","CLEC7A",
					 "AGER","SSC4D","CD207"))
  mat2 <- mat2[geneOrder,]

  p <-Heatmap(mat2,
              cluster_rows = FALSE,cluster_columns = FALSE,
              show_column_names= FALSE,show_row_names = TRUE,
              column_split = all_od$category,
              top_annotation = col_ha,
              column_title = "",
              row_names_max_width = unit(14, "cm"), #设置行名的最大宽度
              row_names_side = "right",
              column_names_side = 'top',
              col = col_fun,
              rect_gp = gpar(col="white"),
              border = TRUE,
              column_order = colnames(mat),
			  row_names_gp = gpar(fontsize = 16),
			  column_names_gp = gpar(fontsize=16),
			  row_title_gp = gpar(fontsize=16),
			  column_title_gp = gpar(fontsize=16))
  draw(p)
  dev.off()
    
  
  mat1 <- mat[c(1:5),]
  pdf(paste0("./kegg_selectedGenes.pdf"),width=20,height=4) #Scavenger
  geneOrder <- rev(c("SLA-DQB1","SLA-DRA","SLA-DRB1","IL6","CXCL12"))
  mat1 <- mat1[geneOrder,]
  p <-Heatmap(mat1,
              cluster_rows = FALSE,cluster_columns = FALSE,
              show_column_names= FALSE,show_row_names = TRUE,
              column_split = all_od$category,
              top_annotation = col_ha,
              column_title = "",
              row_names_max_width = unit(14, "cm"), #设置行名的最大宽度
              row_names_side = "right",
              column_names_side = 'top',
              col = col_fun,
              rect_gp = gpar(col="white"),
              border = TRUE,
              column_order = colnames(mat),
			  row_names_gp = gpar(fontsize = 16),
			  column_names_gp = gpar(fontsize=16),
			  row_title_gp = gpar(fontsize=16),
			  column_title_gp = gpar(fontsize=16))
  draw(p)
  dev.off()
     
  
  
  
  
  
  row_ord <- as.data.frame(row_order(p))
  idx <- row_ord$`row_order(p)`
  ord_mat <- mat[idx,]
  ord_mat <- as.data.frame(ord_mat)
  ord_mat$gene <- rownames(ord_mat)
  ord_mat$pathway <- pathway
  ord_mat <- ord_mat %>% select(394:395,1:393)
  write.xlsx(ord_mat,paste0("./",pathway,"-meanExp-rmsn-sel-b1-Kegg0831.xlsx"),rownames=FALSE,colnames = TRUE)


