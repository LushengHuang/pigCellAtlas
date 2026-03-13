#
rm(list=ls());options(stringsAsFactors=FALSE)
setwd("~/04_2022_scATLAS/202510_EC_immuneGene_HeatMap")

#if (!requireNamespace("BiocManager", quietly = TRUE))
#install.packages("BiocManager")
# The following initializes usage of Bioc devel
#BiocManager::install(version='devel')
#BiocManager::install("M3C")
#BiocManager::install("umap")

library("umap")
library("ggplot2")
#source("umap.r")
#load("pollen_test_data.RData")
#umap(pollen$data,colvec=c('skyblue'))
#umap(pollen$data,labels=as.factor(pollen$celltypes),controlscale=TRUE,scale=3)

load("ave.Rdata")
organ2 <- gsub("Accessory genitalg land","Accessory genital gland",res$organ)

brainOrgans <- c("Olfactory bulb","Cerebral cortex","Striatum",
                 "Corpus callosum","Thalamus","Hippocampus",
                 "Cingulate gyrus","Hypothalamus","Cerebellum",
                 "Brain stem","Spinal cord","Retina","Pituitary","Pineal body")

organ2[organ2 %in% brainOrgans] <- "Brain"
celltype8 <- as.character(res$celltype8)
vesseltype <- sapply(celltype8,function(x){strsplit(x," ")[[1]][1]})

vesseltype[(vesseltype %in% c("Art_Cap","Vein_Cap","Endocardial"))] <- "Other"
vesseltype[(vesseltype %in% c("Central","Portal","HEV"))] <- "Vein"
vesseltype[(vesseltype %in% c("LSEC","Aerocytes"))] <- "Cap"

rownames(res) <- paste0(res[,2],"-",res[,1])
dat <- res[,-c(1,2)]

SDs <- apply(dat,2,sd)
qthr <- quantile(SDs,1 - (2000/ncol(dat)))
gene2select <- which(SDs > qthr)
dat <- dat[,gene2select]

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

uniOrgans <- unique(organ2)
organColors <- godsnot_102[1:length(uniOrgans)]
names(organColors) <- uniOrgans

uniCelltypes <- unique(vesseltype)
vesselColors <- godsnot_102[1:length(uniCelltypes)]
names(vesselColors) <- uniCelltypes

cols = organColors[organ2]
pchs = as.numeric(nrow(dat))

pchs[vesseltype=="Art"] = 19
pchs[vesseltype=="Cap"] = 18
pchs[vesseltype=="Lym"] = 2
pchs[vesseltype=="Vein"] = 17
pchs[vesseltype=="Other"] = 8

umapRes = umap(as.matrix(dat))
scores = data.frame(umapRes$layout)
colnames(scores) = c("UMAP1","UMAP2")

organLegends = names(organColors)
organColors = organColors
organPch <- 15

vesselLegends <- names(vesselColors)
vesselColors <- "grey"
vesselPch <- c(19,18,2,17,8)

tiff("UMAP_412ECsubtypes_organ.tif",width=6000,height=3200,res=600,compression="lzw")
layout(matrix(1:3,nrow=1),widths=c(1,0.6,0.3),height=2)
#par(lab=c(5,10,0.1))
par(mai=c(0.7,0.7,0.1,0.1),omi=c(0,0,0,0),mgp=c(2,0.5,0),tcl=-0.3)
plot(scores[,1],scores[,2],xlab="UMAP1",ylab="UMAP2",
     col=cols,pch=pchs,cex=1.5,lwd=1.2,cex.lab=1.8,cex.axis=1.8)
par(mai=c(0.7,0,0.1,0.1))
plot(x=0,y=0,axes=F,type="n",col="white",xlab="",ylab="")
legend("topleft",legend=organLegends,col=organColors,pch=organPch,bty='n',cex=1.3,text.font=1,ncol=2)
plot(x=0,y=0,axes=F,type="n",col="white",xlab="",ylab="")
legend("topleft",legend=vesselLegends,col=vesselColors,pch=vesselPch,bty='n',cex=1.3,text.font=1)
dev.off()




























devStages = rep(c("adult","embryo"),length(colorInfo))
legends = paste0(legends,"-",devStages)
colorLegends = rep(colorInfo,each=2)
pchLegends = rep(c(19,17),length(colorInfo))

dat = read.csv(file="4CHIP_accumulate_len_expressed.csv",header=T)
colorInfo = c("dark blue","dark green","red","dark orange","Chocolate4")
names(colorInfo) = c("Brain","Heart","Liver","Muscle","SI")
rownames(dat) = dat[,"ID"]
rownames(dat) = sub(":","-",rownames(dat))
peaks10k = read.table(file="PhaseI_peak.txt",header=T)[,1]
peaks10k = peaks10k[peaks10k %in% rownames(dat)]

dat = dat[peaks10k,]
dat2plot = dat[,-c(1:7)]
dat2plot = log2(dat2plot+1)
#dat2plot = t(scale(t(dat2plot)))

tissues = sapply(colnames(dat2plot), function(x){strsplit(x,"_")[[1]][2]})
devs = substring(colnames(dat2plot),1,1)





# PCA analysis
pca = prcomp(dat2plot,scale=T)
tiff("PCA_CHIP_70samples.tif",width=4500,height=3200,res=600,compression="lzw")
layout(matrix(1:2,nrow=1),widths=c(1,0.35),height=2)
#par(lab=c(5,10,0.1))
par(mai=c(0.7,0.7,0.1,0.1),omi=c(0,0,0,0),mgp=c(2,0.5,0),tcl=-0.3)
plot(pca$rotation[,1],pca$rotation[,2],xlab="PC1",ylab="PC2",
     col=cols,pch=pchs,cex=1.5,lwd=1.2,cex.lab=1.2,cex.axis=1.2)
par(mai=c(0.7,0,0.1,0.1))
plot(x=0,y=0,axes=F,type="n",col="white",xlab="",ylab="")
legend("topleft",legend=legends,col=colorLegends,
       pch=pchLegends,bty='n',cex=1,text.font=1)
dev.off()

# For RNA-Seq data
rna = read.table(file="vst_400_threshold.txt",header=T,row.names=1)
colnames(rna) = sub("_wasp_sorted.bam","",colnames(rna))
dat2plot = rna

rnaID = read.table(file="RNA_ID.txt",header=F)[,1]
devIdx = substring(rnaID,1,1)
rnaID[devIdx=="a"] = sub("a","",rnaID[devIdx=="a"])
dat2plot = dat2plot[,rnaID]

tissues = sapply(colnames(dat2plot), function(x){strsplit(x,"_")[[1]][2]})
devs = substring(colnames(dat2plot),1,1)
devs[devs!="e"] = "a"
cols = colorInfo[tissues]
pchs = as.numeric(length(devs))

pchs[devs=="a"] = 19; pchs[devs=="e"] = 17
umapRes = umap(t(as.matrix(dat2plot)))
scores = data.frame(umapRes$layout)
colnames(scores) = c("UMAP1","UMAP2")

legends = rep(names(colorInfo),each=2)
devStages = rep(c("adult","fetal"),length(colorInfo))
legends = paste0(legends,"-",devStages)
colorLegends = rep(colorInfo,each=2)
pchLegends = rep(c(19,17),length(colorInfo))

tiff("UMAP_RNA_VST_60samples.tif",width=4500,height=3200,res=600,compression="lzw")
layout(matrix(1:2,nrow=1),widths=c(1,0.35),height=2)
#par(lab=c(5,10,0.1))
par(mai=c(0.7,0.7,0.1,0.1),omi=c(0,0,0,0),mgp=c(2,0.5,0),tcl=-0.3)
plot(scores[,1],scores[,2],xlab="UMAP1",ylab="UMAP2",
     col=cols,pch=pchs,cex=1.5,lwd=1.2,cex.lab=1.2,cex.axis=1.2)
par(mai=c(0.7,0,0.1,0.1))
plot(x=0,y=0,axes=F,type="n",col="white",xlab="",ylab="")
legend("topleft",legend=legends,col=colorLegends,
       pch=pchLegends,bty='n',cex=1.4,text.font=1)
dev.off()


# PCA analysis
pca = prcomp(dat2plot,scale=T)
tiff("PCA_RNA_65samples.tif",width=4500,height=3200,res=600,compression="lzw")
layout(matrix(1:2,nrow=1),widths=c(1,0.35),height=2)
#par(lab=c(5,10,0.1))
par(mai=c(0.7,0.7,0.1,0.1),omi=c(0,0,0,0),mgp=c(2,0.5,0),tcl=-0.3)
plot(pca$rotation[,1],pca$rotation[,2],xlab="PC1",ylab="PC2",
     col=cols,pch=pchs,cex=1.5,lwd=1.2,cex.lab=1.2,cex.axis=1.2)
par(mai=c(0.7,0,0.1,0.1))
plot(x=0,y=0,axes=F,type="n",col="white",xlab="",ylab="")
legend("topleft",legend=legends,col=colorLegends,
       pch=pchLegends,bty='n',cex=1,text.font=1)
dev.off()
















