rm(list=ls());options(stringsAsFactors=FALSE)
library(tidyverse)
library(dplyr)
library(ggplot2)
library(ggsci)
library(RColorBrewer)
wd = '/inputfile/'
setwd(wd)
dat <- read.csv("amino_acid_composition.csv")

long_data <- dat %>% 
  pivot_longer(cols = c(ACTN3, MYO18B, MYL1,ACTC1,ACTB,ACTG1), 
               names_to = "Genes", 
               values_to = "Value")

long_data$Amino_acids=factor(long_data$Amino_acids,levels=c('Leu','Ile','Lys','Met','Phe','Thr','Tyr','Val','Ala','Arg','Asn','Asp','Cys','Gln','Glu','Gly','His','Pro','Ser','Trp','Pyl','Sec'))

#datM$Age=factor(datM$Age,levels=c('FP75','FP90','NP1','GP33','GP75','BW','AW','PP','OA'))
long_data$Genes=factor(long_data$Genes,levels=c('ACTN3','MYL1','MYO18B','ACTC1','ACTB','ACTG1'))

#colorPallete <- c("#FFFFCC","#FFEDA0","#FED976","#FEB24C","#9E9AC8","#807DBA","#C6DBEF","#9ECAE1","#6BAED6","#4292C6","#2171B5","#08519C","#08306B",'#FFE1FF','#EED2EE','#FFBBFF','#CD96CD','#FF83FA','#D15FEE','#CD69C9','#7A378B','#8B008B')#SampleID_new
colorPallete <- c('#FFCECE','#F8A8A8','#EC6D6D','#E25858','#D93D3D','#B22222','#7A0F0F','#5E0A0A','#F0F8FF','#CCE4FF','#B3D9FF','#99CCFF','#80BFFF','#66AFFF','#4D94FF','#0066FF','#0052CC','#0047AB','#003D99','#003366','#001933','#000C1A')


pdf("/outputfile/AA_composition.pdf",width=6.5,height=4.5)

ggplot(long_data,aes(x=factor(Genes),y=Value)) +
   geom_bar(aes(fill=Amino_acids,order=desc(Amino_acids)),stat = 'identity',position="fill",width = 0.6)+
   scale_fill_manual(values = colorPallete)+
   labs(y="The amino acids composition \n of actin and myosin",x="Genes")+
   theme(text=element_text(size=16)) +
   theme(panel.background = element_rect(fill='white', colour='black'), 
                     panel.grid=element_blank(), 
                     axis.title = element_text(color='black',size=16),axis.ticks.length = unit(0.4,"lines"),   #0.4表示x轴的胡须线
                     axis.line = element_line(colour = "black"), 
                     axis.title.x=element_text(colour='black', size=18),  
                     axis.title.y=element_text(colour='black', size=18),   
                     axis.text=element_text(colour='black',size=14),      
                     legend.title=element_blank(),    
                     legend.text=element_text(size=13),  
                     axis.text.x = element_text(angle = 50,hjust = 1,vjust=1))+ 
    guides(colour = guide_legend(override.aes = list(size=5)))
dev.off()



