rm(list=ls());options(stringsAsFactors=FALSE)
wd <- "~path/09_crossSpecies/cNMF"
setwd(wd)

dat <- read.csv(file="cellInfo.csv",header=TRUE)
rownames(dat) <- dat[,"X"]
dat <- dat[,-1]

usageMat <- read.table(file="CrossSpecial_cNMF/CrossSpecial_cNMF.usages.k_17.dt_0_1.consensus.txt",header=TRUE)

dat <- data.frame(dat,usageMat)

uniCelltypes <- unique(dat$subannotation)
uniSpecies <- unique(dat$species)
i <- uniCelltypes[1]
j <- uniSpecies[1]
res <- NULL

for(i in uniCelltypes){
    for(j in uniSpecies){
        tmp <- subset(dat,subannotation==i & species==j)
        if(nrow(tmp)>10){
            aveUsage <- c(celltype=i,species=j,colMeans(tmp[,paste0("X",1:17)]))
            res <- rbind(res,aveUsage)
        }
    }
}

res <- data.frame(res)
for(i in paste0("X",1:17)){
    res[,i] <- as.numeric(res[,i])
}

res1 <- subset(res,!(celltype %in% c("Adipose macrophages","Doublets")))
res1 <- res1[order(res1$celltype),]


uniqueCelltype <- unique(res1$celltype)

celltypeMat <- matrix(NA,length(uniqueCelltype),17)
rownames(celltypeMat) <- uniqueCelltype

for(i in uniqueCelltype){
    celltypeMat[i,] <- colMeans(subset(res1,celltype==i)[,3:19])
}

hclustSub <- hclust(as.dist(1 - cor(t(celltypeMat))))

celltypeOrd <- uniqueCelltype[hclustSub$order]

# order the res1
res2 <- NULL
for(i in celltypeOrd){
    res2 <- rbind(res2,subset(res1,celltype==i))
}

# order modules
celltypeMat <- celltypeMat[celltypeOrd,]
max_col_index <- apply(celltypeMat,1,which.max)
moduleOrd <- c("celltype","species",paste0("X",c(unique(max_col_index),17)))
df <- res2[,moduleOrd]
write.table(df,file="df.txt",sep="\t")


library(openxlsx)
library(tidyverse)
library(circlize)   
library(RColorBrewer)
library(ggsci)
library(ComplexHeatmap)
library(plyr)

df <- read.table(file="df.txt",header=T)
df2 <- df[,-c(1,2)]

allCelltypes = c('Epithelial cells', 'B cells', 'Doublets', 'Fibroblasts', 'Neutrophils', 
                'Proximal tubule cells', 'Smooth muscle cells', 'Adipose fibroblasts', 
                'Liver macrophages', 'Granule cells', 'Oligodendrocytes', 'Stem cells', 
                'Plasma cells', 'Enterocytes', 'Pancreas exocrine cell', 'Skeletal muscle cells', 
                'Dendritic cells', 'Alveolar type II cells', 'T cells', 'Interstinal epithelial cells', 
                'Endothelial cells', 'Stomach epithelial cells', 'Adipose macrophages', 'Adipocytes', 
                'VEC', 'Lung macrophages', 'Hepatocytes', 'Gastric chief cells', 
                'Adrenal Gland epithelial cells', 'Proliferating T cells', 'Macrophages', 
                'Proliferating B cells', 'Astrocytes')

celltypeColors <- c("#E64B35FF","#4DBBD5FF","#00A087FF","#3C5488FF","#F39B7FFF","#8491B4FF","#91D1C2FF","#DC0000FF","#7E6148FF","#B09C85FF", 
            "#3B4992FF","#EE0000FF","#008B45FF","#631879FF","#008280FF","#BB0021FF","#5F559BFF","#A20056FF","#808180FF","#BC3C29FF",
            "#0072B5FF","#E18727FF","#20854EFF","#7876B1FF","#6F99ADFF","#FFDC91FF","#EE4C97FF","#00468BFF","#ED0000FF","#42B540FF",
            "#0099B4FF","#925E9FFF","#FDAF91FF")

names(celltypeColors) <- c("Adipocytes","Adipose fibroblasts","Adipose macrophages","Adrenal Gland epithelial cells","Alveolar type II cells",
                           "Astrocytes","B cells","Dendritic cells","Doublets","Endothelial cells",
                           "Enterocytes","Epithelial cells","Fibroblasts","Gastric chief cells","Granule cells",
                           "Hepatocytes","Interstinal epithelial cells","Liver macrophages","Lung macrophages","Macrophages",
                           "Neutrophils","Oligodendrocytes","Pancreas exocrine cell","Plasma cells","Proliferating B cells",
                           "Proliferating T cells","Proximal tubule cells","Skeletal muscle cells","Smooth muscle cells","Stem cells",
                           "Stomach epithelial cells","T cells","VEC")


uniqueCelltypes <- unique(df$celltype)
panel_celltype <- celltypeColors[uniqueCelltypes]
panel_species <- c("P"="#00A087FF", "H"="#E64B35FF", "M"="#4DBBD5FF")

df$celltype <- factor(df$celltype,levels=uniqueCelltypes)
df$species <- factor(df$species, levels=c("P","H","M"))

col_ha <- HeatmapAnnotation(celltype = df$celltype,
                              species = df$species,
                              col = list(celltype = panel_celltype,
                                         species = panel_species))
  mat <- t(df2) 
  max(mat)
  min(mat)
  
  
  #col_fun <- structure(c("white","red"), names = c("0", "1")) # black, red, green, blue
  col_fun  <- colorRamp2(c(0,ceiling(max(mat)/2)),c("#FFFFFF","#E0002C"))
  
  pdf(paste0("crossSpecies_network2.pdf"),width=10,height=5.6) #Scavenger
  p <-Heatmap(mat,
              cluster_rows = FALSE,cluster_columns = FALSE,
              show_column_names= FALSE,show_row_names = TRUE,
              top_annotation = col_ha,
              #column_split = df$celltype,
              column_title = "",
              row_names_max_width = unit(14, "cm"), #设置行名的最大宽度
              row_names_side = "right",
              column_names_side = 'top',
              col = col_fun,
              rect_gp = gpar(col="grey90"),
              border = TRUE,
              column_order = colnames(mat),
              row_names_gp = gpar(fontsize = 14),
              column_names_gp = gpar(fontsize=14),
              row_title_gp = gpar(fontsize=14),
              column_title_gp = gpar(fontsize=14))
  draw(p)
  dev.off()



