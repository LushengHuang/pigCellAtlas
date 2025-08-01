##Compare the proportion of type II myofiber C1 cluster across different groups

rm(list=ls())
library(dplyr)
library(ggpubr)
library(ggplot2)

setwd('res0.3')

#Calculate the proportion of type II myofiber C1 subtype in each sample
dat_all = read.csv('harmony_sampleID_HumanMuscleCell_res0.3_celltype_cellInfo.csv')
df1 <- dat_all %>%
       group_by(SampleID) %>% summarise(Nb = n())
colnames(df1)[2] <- 'nCells'

data=read.csv('Ingest_Human-pig_myofiberIIsubtypes_cellInfo.csv')
df=data[data$Project %in% 'Human',]
df=df[df$celltype %in% '1',]

df$Group <- df$SampleID
df$Group[df$Group == "SMCSN_1mg_1"] <- "1mg"
df$Group[df$Group == "SMCSN_1mg_2"] <- "1mg"
df$Group[df$Group == "SMCSN_1mg_3"] <- "1mg"
df$Group[df$Group == "SMCSN_10mg_1"] <- "10mg"
df$Group[df$Group == "SMCSN_10mg_2"] <- "10mg"
df$Group[df$Group == "SMCSN_10mg_3"] <- "10mg"
df$Group[df$Group == "SMCSN_105mg_1"] <- "105mg"
df$Group[df$Group == "SMCSN_105mg_2"] <- "105mg"
df$Group[df$Group == "SMCSN_105mg_3"] <- "105mg"

df2 <- df %>%
       group_by(Group,SampleID,celltype) %>% summarise(Nb = n())
dat_merge = merge(df1,df2,by='SampleID')
dat_merge <- dat_merge %>% mutate(Percent = Nb/nCells*100)

group_colors <- c("1mg" = "#4DBBD5FF", "10mg" = "#00A087FF", "105mg" = "#F39B7FFF")
comparisons <- list(c("1mg", "10mg"), c("10mg", "105mg"), c("1mg", "105mg"))
dat_merge$Group <- factor(dat_merge$Group, levels = c("1mg", "10mg", "105mg"))
#Visualization
p <- ggplot(dat_merge, aes(x = Group, y = Percent, fill = Group, color = Group)) +
  geom_boxplot(width = 0.5, outlier.shape = NA, alpha = 0.4) +  
  geom_jitter(width = 0.15, size = 2, shape = 21, stroke = 1) + 
  scale_fill_manual(values = group_colors) +
  scale_color_manual(values = group_colors) +
  stat_compare_means(comparisons = comparisons, method = "t.test", label = "p.format",bracket.size = 0.6,size=4) +
  labs(title = '',x = "Group", y = "Proportion of type II C1 (%)") +
  #ylim(-2.5, 0) +
  theme_bw(base_size = 14) +
  theme(panel.grid = element_blank(),      
        axis.text = element_text(color = "black"))+
  theme(axis.line = element_line(linewidth = 0.7), 
        axis.ticks = element_line(linewidth = 0.7),legend.position="none")+ 
  theme(axis.text.x= element_text(size=14,colour="black"),axis.text.y = element_text(size=14,colour="black"),axis.title.y = element_text(size=16))
ggsave('Boxplot_Human_MyofiberII_C1_Compare_P1-P10-P105.pdf',p,width =3,height = 4)