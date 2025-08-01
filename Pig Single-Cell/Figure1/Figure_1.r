library(reshape2)
library(ggplot2) 
library(ggsci)
library(gridExtra)
library(stringr)
library(tidyverse)
library(aplot)

#------------------ Set parameters ------------------#
setwd('Pathway')
Plotdata = read.table(file = 'All_sample_info.csv', sep = ',',header = T)
#------------------ Preprocessing ------------------#
for (i in 1:14) {
  tmp = aggregate(Plotdata[,i], by=list(type=Plotdata[,1]),sum)
  ylim=max(tmp[,2])+100
  p=ggplot(Plotdata)+
  geom_bar(aes(y=TissueName,x=Plotdata[,i],colour=Platform,fill=Platform),stat = 'identity',width=0.5,size = 0.05)+
  scale_colour_manual(values=c("#f85f73","#283c63")) +
  coord_fixed(ratio=2000) +
  scale_fill_manual(values=c("#f85f73","#283c63"))+
  theme_bw()+
  theme(panel.border = element_blank(), 
        axis.ticks.y = element_blank(),
        panel.grid.major =element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text.y =  element_blank(),
        axis.title =  element_blank(),
        axis.text.x = element_text(size=10,angle =45,hjust =0,vjust = 0),  
        plot.margin=unit(rep(3,4),'lines'),
        legend.position = 'none')+      
  scale_x_continuous(breaks=seq(0, 33000, 5000), labels = seq(0, 33000, 5000),expand = c(0,0),limits = c(0,ylim), position="top")+
  xlab('Cell Number before QC')+
  ylab('Anatomic sites')
#p
ggsave(paste0(colnames(Plotdata)[i],"_all_tissue_cell_number_after_QC_220530.pdf"),p, width = 20, height = 20)
}

#Samples are classified according to the system
Plotdata$TissueName <- gsub("_", " ", Plotdata$TissueName)
sys <- read.table(file = 'Tissue_system_class_info.csv', sep = ',',header = T)
sys <- sys[order(sys$order,decreasing = TRUE),]
Plotdata$TissueName <- factor(Plotdata$TissueName, levels = sys$TissueName)

sys <- sys %>%
  mutate(p="System")

sys$TissueName <- factor(sys$TissueName, levels = sys$TissueName)
sys <- sys[order(sys$order,decreasing = F),]
sys$System_abbr <- factor(sys$System_abbr, levels = unique(sys$System_abbr))

type <- ggplot(sys,aes(p,TissueName,fill=System_abbr))+
  geom_tile() +
  #geom_tile(color = "black",size=0.1) +
  coord_fixed(ratio=1)+
  scale_fill_d3(palette = "category20")+
  scale_y_discrete(position="left") +
  xlab(NULL) + ylab(NULL) +
  theme(#axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid = element_blank())+
  labs(fill = "System")
type
ggsave("tissue_type.pdf",type, width = 20, height = 20)

p2=p %>% insert_left(type, width=.25)
p2
ggsave(paste0(colnames(Plotdata)[i],"_all_tissue_cell_number_after_QC_220530_bar2.pdf"),p2, width = 20, height = 20)
