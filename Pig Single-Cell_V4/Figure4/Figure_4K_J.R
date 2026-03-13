## in python
import scanpy as sc
import pandas as pd
import os

workdir = "./"
os.chdir(workdir)
dat = sc.read('bbknn_ratio_1_resolution_mirco_addcelltype1.h5ad')
gene_list = ['AIF1', 'ALOX5', 'ARRB1', 'ATF4', 'B2M', 'C3', 'CD33', 'CD81', 'CSF1R', 'CSK', 'CX3CR1', 'CYBA', 
'FGR', 'HFE', 'HSP90AA1', 'HSP90AB1', 'HSPD1', 'IL6ST', 'IRF8', 'KLF2', 'LAPTM5', 'MYD88', 'NDRG2', 'PLD4', 'POLR3B', 'RNF128', 
'SASH3', 'SRGN', 'TGFB1', 'TMSB4X', 'TNFRSF1B', 'TREM2']

sc.tl.score_genes(dat, gene_list=gene_list)
dat.obs["celltype1-State"] = dat.obs["celltype1"].str.cat(dat.obs["State"], sep="-")
result = pd.DataFrame({"celltype1-State": dat.obs["celltype1-State"], "score": dat.obs["score"]})
pd.DataFrame(result).to_csv("regulation_of_cytokine_production_score_mat.csv")



## in R
rm(list = ls())
library(Seurat)
library(dplyr)
library(ggplot2)
library(tidyverse)

workdir = "./"
setwd(workdir)
mat = read.csv("regulation_of_cytokine_production_score_mat.csv")
mat2 = separate(data = mat, col = celltype1.State, into = c("celltype1","State"), sep = "-")

mat3 = mat2[mat2$celltype1 == "Microglia1",]
mat3$State = factor(mat3$State, levels = c("WT","KO1","KO2"))
e <- ggplot(mat3, aes_string(x = "State", y = "score")) + 
  #geom_violin(width = 0.5,aes(fill = disease_state), trim = FALSE) + 
  scale_fill_manual(values=c("#ff6b31","#135f9b","#0088c0"))+
  labs(x="",y="")+
  geom_boxplot(width = 0.6,aes(fill = State))+
     theme(panel.background = element_rect(fill='white', colour='black'), 
                     panel.grid=element_blank(), 
                     axis.title = element_text(color='black',
                     size=16),axis.ticks.length = unit(0.4,"lines"),   #0.4表示x轴的胡须线
                     axis.line = element_line(colour = "black"), 
                     axis.title.x=element_text(colour='black', size=16),  
                     axis.title.y=element_text(colour='black', size=16),   
                     axis.text=element_text(colour='black',size=16),      
                     legend.title=element_blank(), 
                     legend.text=element_text(size=16),
                     axis.text.x = element_text(angle = 45,hjust = 1,vjust=1),
                     plot.margin = margin(10, 10, 10, 100))+ #plot.margin
     stat_compare_means(comparisons = list(c("WT", "KO1"),c("WT","KO2")))
ggsave('regulation_of_cytokine_production_genes_state_M1.pdf',e,height=4,width=6)

