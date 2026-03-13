# Clear current R environment
rm(list = ls())

# --------------------- Load required libraries ---------------------#
library(tidyverse)
library(dplyr)
library(ggplot2)
library(ggsci)
library(rstatix)
library(ggpubr)
library(ggtext)

# --------------------- Load and filter dataset ---------------------#
subdata <- read.csv('input_files.csv')
subdata$GROUP3 <- factor(subdata$GROUP3, levels=c('Control','GR'))

# ----------- Create ggplot: visualize the percentage of type II myofiber C1 subset across groups ------------#
p1 <-ggplot(subdata, aes(x = GROUP3,y = Percent))+
  geom_bar(position = 'dodge', stat = 'summary', fun = 'mean',fill = "#CDCCCD") +
  geom_point(aes(x = GROUP3,fill = Family,group=Family), shape = 21,size=3)+
  labs(x = "",y=paste0("Percentage of Type II (C1)"))+
  geom_line(aes(color=Family, group=Family))+
  scale_fill_manual(values =  c('#ea1c1c','#356db1','#4bb049','#4bb049'))+
  scale_y_continuous(expand=c(0,0.11),limits=c(0,26))+
  scale_color_manual(values = colorPallete)+
  theme(panel.background = element_rect(fill='white', colour='black'),
        panel.grid=element_blank(),
        axis.text = element_text(color="black",size = 15),
        axis.title.y = element_text(color="black",size = 15),
        legend.title=element_text(color='black',size = 15),
        legend.text=element_text(color='black', size = 15))
ggsave(paste0("TypeII c1_pair.pdf"),plot = p1,width = 4,height = 3)
