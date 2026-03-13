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
subdata <- read.csv('input_file.csv', sep=",")

#------- Convert the dataaset from wide to long format for ggplot -------#
subdata_long <- subdata %>%
  pivot_longer(cols = -c(GROUP3, broadcelltype1), names_to = "Gene", values_to = "expression")

#--------------------- Calculate statistical significance ---------------------#
compare_means(exression~GROUPs,data=subdat_long,group.by='Gene',method="t.test")

# --------------------- Create ggplot: visualize gene expression across groups ---------------------#
plt <- ggplot(subdata_long, aes(x = GROUP3, y = expression, colour = GROUP3)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.8, dodge.width = 0.2),size=0.5) + 
  stat_summary(aes(group = GROUP3), fun = mean, geom = "point", 
               shape = 18, size = 3, color = "black", position = position_dodge(width = 0.75)) + 
  stat_compare_means(aes(group = GROUP3), label = "p.signif", method = "t.test")+
  scale_y_break(breaks =c(2,15),scales = "free")+
  scale_colour_manual(values =  c("#939393","#C82D35"))+
  labs(x = "", y = "expression") +
  facet_grid(. ~ Gene) +
  theme_classic()
ggsave(paste0("Trophoblasts_SLC_gene_plot.pdf"),plt,width=4,height=3,onefile=FALSE)

