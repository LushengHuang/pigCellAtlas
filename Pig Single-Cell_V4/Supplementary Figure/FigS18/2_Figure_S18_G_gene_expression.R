# Clear the current R environment
rm(list=ls())

# --------------------- Load required libraries ---------------------#
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggpubr)
library(rstatix)
library(stringr)
library(ggbreak)

# --------------------- Load and filter dataaset ---------------------#
subdata<- read.csv('/file/TypeIImyofiber_MYH_protein_synthesis_1_Group_4v4.csv', head=T, sep=",")
setwd("/inputfile/")

#------- Convert the dataaset from wide to long format for ggplot -------#
long_data <- subdata %>% 
  pivot_longer(cols = c(ACTN3,MYO18B,MYL1,ACTC1,MYO9A,ACTG1,ACTB,MYH8,MYH3,MYL2,MYH14,MYO1H,MYO1E,MYO6,MYO5A,MYO7A,MYO19,MYO1B,MYL10,MYL3,MYH7,MYO9B),names_to = "Genes", 
               values_to = "Value")

long_data$Genes=factor(long_data$Genes,levels=c('ACTN3','MYO18B','MYL1','ACTC1','MYO9A','ACTG1','ACTB','MYH8','MYH3','MYL2','MYH14','MYO1H','MYO1E','MYO6','MYO5A','MYO7A','MYO19','MYO1B','MYL10','MYL3','MYH7','MYO9B'))
               

long_data <- subdata %>% 
  pivot_longer(cols = c(MYO1F,ACTA2,ACTN1,MYO1D,MYO5C,MYO10,MYO7B,MYH15,MYL5,MYO1A,ACTN2,ACTG2,MYH13,MYO1C,MYL9,MYL4,MYL12A,ACTN4,MYO1G,MYH7B,MYO3B,MYH11),names_to = "Genes", 
               values_to = "Value")

long_data$Genes=factor(long_data$Genes,levels=c('MYO1F','ACTA2','ACTN1','MYO1D','MYO5C','MYO10','MYO7B','MYH15','MYL5','MYO1A','ACTN2','ACTG2','MYH13','MYO1C','MYL9','MYL4','MYL12A','ACTN4','MYO1G','MYH7B','MYO3B','MYH11'))
     


#--------------------- Calculate statistical significance ---------------------#
compare_means(Value~GROUP3,data=long_data,group.by='Genes',method="t.test")

# --------------------- Create ggplot: visualize gene expression across groups ---------------------#
plt <- ggplot(long_data, aes(x = GROUP3, y = Value, colour = GROUP3)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.8, dodge.width = 0.2),size=0.5) + 
  stat_summary(aes(group = GROUP3), fun = mean, geom = "point", 
               shape = 18, size = 3, color = "black", position = position_dodge(width = 0.75)) + 
  stat_compare_means(aes(group = GROUP3), label = "p.signif", method = "t.test")+
  #scale_y_break(breaks =c(2,15),scales = "free")+
  scale_colour_manual(values =  c("#939393","#8B0000"))+
  labs(x = "", y = "The expression of \n actin and myosin genes") +
  facet_grid(. ~ Genes)+
  theme_classic()
ggsave(paste0("22-44_gene_plot.pdf"),plt,width=13,height=3.5,onefile=FALSE)

