library(ggdendro)
library(Seurat)
library(tidyverse)
library(ggplot2)
library(cowplot)
library(RColorBrewer)

wd <- "~/endothelial/mnn_02"
outdir <- "~/endothelial/mnn_02/cluster_organ_barplot/" 
setwd(wd); options(stringsAsFactors=FALSE)

avgs <- read.csv("mnn_res1_datExpr.csv",row.names=1,header=T)
meta_data <- read.csv("endothelial1_mnn_cellInfo.csv",row.names=1,header=T)
avgs <- t(avgs)
avgs <- avgs[rowSums(avgs)>0,]
meanRank <- rank(rowSums(avgs))
avgs <- avgs[meanRank > 0.5*nrow(avgs),]
sdRank <- rank(apply(avgs,1,sd))
avgs <- avgs[sdRank>max(sdRank-2000),]

organ <- meta_data[,"Organ"]
unique(organ)
organ[organ=='Corpus callosum'] <- 'Brain'
organ[organ=='Hippocampus'] <- 'Brain'
organ[organ=='Cerebral cortex'] <- 'Brain'
organ[organ=='Pituitary'] <- 'Brain'
organ[organ=='Olfactory bulb'] <- 'Brain'
organ[organ=='Thalamus'] <- 'Brain'
organ[organ=='Cingulate gyrus'] <- 'Brain'
organ[organ=='Striatum'] <- 'Brain'
organ[organ=='Cerebellum'] <- 'Brain'
organ[organ=='Pineal body'] <- 'Brain'
organ[organ=='Hypothalamus'] <- 'Brain'
organ[organ=='Brain stem'] <- 'Brain'
organ[organ=='Retina'] <- 'Brain'
organ[organ=='Spinal cord'] <- 'Brain'
meta_data[,"Organ"] <- organ

meta_data[meta_data[,"Individual_2"] %in% c("Sow1_E1_1","Sow1_E1_2"),"Individual_2"] <- "Sow1_E1"
meta_data[meta_data[,"Individual_2"] %in% c("Sow1_1","Sow1_2"),"Individual_2"] <- "Sow1"
inds <- c("Sow1","Sow2","Boar","Sow1_E1","Sow1_E2","Sow1_E9","Sow2_E1")
ind_cols <- c(brewer.pal(9,"YlOrRd")[7:9],brewer.pal(9,'Blues')[6:9])
names(ind_cols) <- inds
meta_data[,"Individual_2"] <- factor(meta_data[,"Individual_2"],levels=inds)

Ind_organ <- paste0(meta_data$Individual, "_", meta_data$Organ)
meta_data <- data.frame(meta_data,Ind_organ)

hclt <- hclust(dist(MinMax(scale(t(avgs)), min = -3, max = 3)),method="ward.D2")
orders <- as.character(hclt$order-1)
meta_data[,"seurat_clusters"] <- factor(meta_data[,"seurat_clusters"],levels=orders)
meta_data[,"integrated_snn_res.1"] <- factor(meta_data[,"integrated_snn_res.1"],levels=orders)

groupby = "seurat_clusters"

# build the dendrogram
                        
theme_empty <- function(x) {
    theme_classic() + theme(axis.line.y=element_blank(),axis.ticks.y=element_blank(),axis.text.y=element_blank(),axis.title=element_blank(),legend.position="none")
}                                                                                                                                                          

# Rectangular lines
dendro <- dist(MinMax(scale(t(avgs)), min = -3, max = 3)) %>% hclust(., method = "ward.D") %>% as.dendrogram()
dendro <- reorder(dendro,ncol(avgs):1, agglo.FUN=mean)
pdf(paste0(outdir,"cluster_dendrogram.pdf"), width = 7, height = 22)
    ggdendrogram(dendro,rotate = TRUE)
dev.off()
           
ddata <- dendro_data(dendro,type='rectangle')
den_p <- ggplot(segment(ddata)) + 
             geom_segment(aes(x=x,y=y,xend=xend,yend=yend)) +
             coord_flip() + 
             scale_y_reverse(expand = c(0,0)) +  theme_empty()

# plot the region contribution
clssize_data  <- meta_data %>% 
                 group_by(seurat_clusters) %>% 
                 summarize(ncells=n()) %>% 
                 mutate(ncells2 = sqrt(ncells)) %>% 
                 ungroup()

indOrgan_data <-  meta_data %>%
                  group_by(seurat_clusters,Ind_organ) %>%
                  summarize(ncells=n()) %>% ungroup() %>% group_by(seurat_clusters) %>%
                  mutate(ratio = ncells/sum(ncells))
write.table(indOrgan_data,file="indOrganContri.txt",sep="\t",col.names=T,row.names=F,quote=F)

regcontri_data <- meta_data %>% 
                  group_by(seurat_clusters,Organ) %>% 
                  summarize(ncells=n()) %>% ungroup() %>% group_by(Organ) %>%
                  mutate(ratio = ncells/sum(ncells)) %>% ungroup() %>% 
                  group_by(seurat_clusters) %>% mutate(rratio = ncells/sum(ncells))

indcontri_data <- meta_data %>% group_by(seurat_clusters,Individual_2) %>% 
                  summarize(ncells=n()) %>% ungroup() %>% group_by(Individual_2) %>%
                  mutate(ratio = ncells/sum(ncells)) %>% ungroup() %>% 
                  group_by(seurat_clusters) %>% mutate(rratio = ratio * 100/sum(ratio))

size_bks <- seq(0,200,50)
xsize = 18
ysize = 18
cls_order <- orders
size_p <- ggplot(clssize_data,aes_string(x = "ncells2",y = "seurat_clusters")) +
          geom_bar(position = 'stack',color=NA, stat = 'identity') + 
          scale_y_discrete(limits=rev(cls_order)) + 
          scale_x_continuous(position='top',breaks=size_bks,labels=as.character(size_bks^2))+
          theme_empty() +
          theme(axis.text.y=element_text(size=ysize,face="bold"),
                axis.text.x=element_text(size=xsize,face="bold",angle=45,hjust=0),
                axis.ticks.x = element_line(size = 0.2), 
                axis.line.x = element_line(size = 0.2))

reg_p <- ggplot(regcontri_data,aes_string(x = "rratio",y="seurat_clusters",fill="Organ")) +
		geom_bar(position = "stack", color = NA, stat = "identity") +
		scale_fill_manual(values = colorRampPalette(brewer.pal(8,'Accent'))(length(unique(organ))) %>% setNames(.,unique(meta_data$Organ))) +
		scale_y_discrete(limits = rev(cls_order))+
		scale_x_continuous(position = "top") +
		theme_empty() +
		theme(axis.text.x=element_text(size = xsize,angle = 45,hjust = 0,face="bold"), 
                      axis.ticks.x=element_line(size = 0.2),axis.line.x=element_line(size=0.2))

ind_p <- ggplot(indcontri_data,aes_string(x = "rratio", y = "seurat_clusters", fill = "Individual_2")) +
		geom_bar(positio = "stack", color = NA, stat = "identity") +
		scale_fill_manual(values=ind_cols) +
		scale_y_discrete(limits = rev(cls_order))+
		scale_x_continuous(position="top") +
		theme_empty() +
		theme(axis.text.x=element_text(size = xsize,angle = 45,hjust = 0,face="bold"),
                      axis.ticks.x=element_line(size = 0.2),axis.line.x=element_line(size=0.2))

pdf(paste0(outdir, "MF1.Cluster_rratio.pdf"), width = 15, height = 20)
plot_grid(size_p,reg_p,ind_p,size_p, nrow = 1, ncol = 3, rel_widths = c(0.2, 0.3, 0.2), align = "h") %>% print()
dev.off()

reg_p <- ggplot(regcontri_data,aes_string(x = "rratio",y="seurat_clusters",fill="Organ")) +
                geom_bar(positio = "stack", color = NA, stat = "identity") +
                scale_fill_manual(values = colorRampPalette(brewer.pal(8,'Accent'))(length(unique(organ))) %>% setNames(.,unique(meta_data$Organ))) +
                scale_y_discrete(limits = rev(cls_order))+
                scale_x_continuous(position = "top") +
                theme_classic() +
                theme(axis.text.x=element_text(size = xsize,angle = 45,hjust = 0,face="bold"),
                      axis.ticks.x=element_line(size = 0.2),axis.line.x=element_line(size=0.2)) +
		guides(fill = guide_legend(ncol=2, byrow = TRUE,reverse = T)) +
		theme(legend.key.size = unit(18, "pt"))
pdf(paste0(outdir,"regcontri.pdf"), width = 12, height = 20)
reg_p
dev.off()

ind_p <- ggplot(indcontri_data,aes_string(x = "rratio", y = "seurat_clusters", fill = "Individual_2")) +
                geom_bar(positio = "stack", color = NA, stat = "identity") +
                scale_fill_manual(values=ind_cols) +
                scale_y_discrete(limits = rev(cls_order))+
                scale_x_continuous(position="top") +
                theme_classic() +
                theme(axis.text.x=element_text(size = xsize,angle = 45,hjust = 0,face="bold"),
                      axis.ticks.x=element_line(size = 0.2),axis.line.x=element_line(size=0.2))
                guides(fill = guide_legend(ncol=1, byrow = TRUE,reverse = T)) +
                theme(legend.key.size = unit(18, "pt"))
pdf(paste0(outdir,"indcontri.pdf"), width = 12, height = 20)
ind_p
dev.off()




