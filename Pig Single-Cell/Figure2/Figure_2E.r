rm(list=ls());options(stringsAsFactors=FALSE)
wd <- 'path'
setwd(wd)

library(tidyverse)
library(dplyr)
library(ggplot2)
library(ggsci)
library(RColorBrewer)

# ---- endothelial1_mnn_cellInfo.csv is the adata.obs, where adata is Scanpy anndata object integrating all ECs
# Prepare the data
dat <- read.csv("endothelial1_mnn_cellInfo.csv",header=T)
dat[dat$Individual %in% c("Sow","Boar"),"Individual"] <- "Adult"
dat$celltype <- "VEC"
dat$celltype[dat$seurat_clusters %in% c(11,17,44,46,48,53,54)] <- "LEC"
test <- dat %>% group_by(System_abbr,Organ,Individual,celltype,seurat_clusters) %>% summarise(counts=n()) %>% 
        ungroup() %>% group_by(Organ,Individual,celltype) %>% 
		mutate(ratio = counts/sum(counts)) %>% mutate(allcounts=sum(counts))
test[test$Organ=="Adipose tissue","System_abbr"] <- "Others"


# Set up the organ orders and colors
organOrder <- c("Olfactory bulb","Cerebral cortex","Striatum","Corpus callosum","Thalamus",
                "Hippocampus","Cingulate gyrus","Hypothalamus","Cerebellum","Brain stem",
				"Spinal cord","Retina",# NS
				"Pituitary","Pineal body","Glands adrenal","Glands thyroid",# ES
				"Blood vessel","Heart", # CS
				"Spleen","Tonsil","Lymph nodes",# IS
				"Turbinate","Vocal cords","Trachea","Throat","Pharynx","Bronchi","Lung","Diaphragm",#RS
				"Salivary gland","Tongue","Esophagus","Pancreas","Liver","Gallbladder",
				"Stomach","Small intestine","Large intestine","Anus",#DS
				"Ovary","Oviduct","Uterus","Urethra","Vagina", #FRS
				"Placenta","Umbilical cord", # MFI
				"Testis","Epididymis","Vas deferens","Accessory genital gland","Corpus cavernosum",
				"Penis",# MRS
				"Kidney","Ureter","Bladder", # US
				"Muscle","Cartilage","Ligament", # MSS
				"Adipose tissue","Mammary glands" # Others
			    )		
test$Organ <- factor(test$Organ, levels=organOrder)

sysOrg <- data.frame(test[!duplicated(test[,1:2]),1:2])
sysOrg <- sysOrg[!duplicated(sysOrg$Organ),]
rownames(sysOrg) <- sysOrg[,"Organ"]
sysOrg <- sysOrg[organOrder,]
#tabSys <- table(sysOrg[,1])
#numSys <- length(tabSys)
#xTextCol <- rep(colorRampPalette(brewer.pal(8,'Accent'))(numSys),tabSys)
panel_sys <- c("#1F77B4FF","#FF7F0EFF","#2CA02CFF","#D62728FF","#9467BDFF",
               "#8C564BFF","#E377C2FF","#7F7F7FFF","#BCBD22FF","#17BECFFF",
               "#AEC7E8FF","#FFBB78FF")
sys <- c("NS","ES","CS","IS","RS",
         "DS","FRS","MFI","MRS","US",
		 "MSS","Others")
		 
names(panel_sys) <- sys
xTextCol <- panel_sys[sysOrg$System_abbr]


for(clt in c(22,26,48)){
    test2 <- test[!duplicated(test[,c(1:4,8)]),c(1:4,8)]
    test2$ratio <- 1e-9
    test2$counts <- 1
    test2$seurat_clusters <- clt
    test2 <- test2[,c(colnames(test2)[1:4],"seurat_clusters","counts","ratio","allcounts")]
    test <- rbind(test,test2)
}

# ----  define a function to draw the barplot ----------------

barplotprop <- function(celltype='LEC',cl=48){
#    celltype="LEC"
#    cl = 48
    d <- subset(test,celltype==celltype & seurat_clusters==cl & allcounts > 50)
    nm <- paste0("C",cl)
    col = c("Adult"='#DD3497','Fetus'='#4EB3D3')
    maxV <- max(abs(d$ratio))
    d[d$Individual=='Fetus','ratio'] <- -d[d$Individual=='Fetus','ratio']
#    d$Organ <- as.factor(d$Organ)
#    levels(d$Organ) <- as.character(sysOrg$Organ)
    d$Organ <- factor(d$Organ,levels = sysOrg$Organ)
    p <- ggplot(d,aes(x=Organ,y=ratio,fill=Individual)) + geom_bar(stat='identity',colour='black') +
            scale_fill_manual(values=col) + geom_hline(yintercept=0) +
            labs(x = '',y = 'Proportion', title=nm) + theme_bw() + 
            theme(panel.grid=element_blank(),
            text = element_text(size=16),
            axis.text.x = element_text(angle = 45,hjust = 1,vjust=1,size = 16,face='bold',colour = xTextCol),
            axis.text.y = element_text(size = 16,colour = "black"),
            legend.text = element_text(size=16, colour = "black"),
            legend.title = element_text(size=16,colour = "black"), 
            plot.margin = margin(10,10,10,100)) + ylim(-max(maxV),max(maxV))
    ggsave(paste0(nm,"_prop.pdf"),p,width=20,height=5)
    ggsave(paste0(nm,"_prop.png"),p,width=20,height=5)
}
barplotprop(celltype='VEC',cl=22)
barplotprop(celltype='VEC',cl=26)
barplotprop(celltype='LEC',cl=48)

