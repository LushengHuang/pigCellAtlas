# Clear the current R environment
rm(list=ls());options(stringsAsFactors=FALSE)

# --------------------- Load required libraries ---------------------#
library(ggpubr)
library(patchwork)
library(ggsci)
library(ggbreak)
library(tidyr)

# --------------------- Load dataset ---------------------#
dfb <- read.csv("input_file.csv",header=T,sep=",")

#------- Convert the dataaset from wide to long format for ggplot -------#
long_data <- dfb %>% 
  pivot_longer(cols = c(ACTN3, MYO18B, MYL1), 
               names_to = "Genes", 
               values_to = "Value")
dfb$Group = factor(dfb$Group, levels=c('1mg','10mg','105mg'))

# --------------------- Create ggplot: visualize gene expression across groups ---------------------#
for(i in 2:3){
  tmpdata=long_data[which(long_data$Genes==paste0(colnames(dfb[,3:5]))[i]),]
  tmpdata$Group = factor(tmpdata$Group, levels=c('1mg','10mg','105mg'))
  my_comparisons <- list( c("1mg","10mg"),c("10mg","105mg"),c("1mg","105mg"))

plt <- ggplot(tmpdata, aes(x = Group, y = Value, colour = Group)) +
     geom_point(position = position_jitterdodge(jitter.width = 0.8, dodge.width = 0.2),size=0.5) + 
    #facet_grid(. ~ Genes) +
      stat_summary(aes(group = Group), fun = mean, geom = "point", 
               shape = 18, size = 5, color = "black", position = position_dodge(width = 0.75)) + 
   #stat_compare_means(aes(group = Group), label = "p.signif", method = "t.test")+  #p.format
  #scale_y_break(breaks =c(10,20),scales = "free")+
  scale_colour_manual(values =  c("#4DBBD5FF", "#00A087FF", "#F39B7FFF")) +
  labs(x = "", y = "expression") +
  theme_classic()+
  theme(axis.line = element_line(linewidth = 0.7), 
        axis.ticks = element_line(linewidth = 0.7),legend.position="none")+ 
  theme(axis.text.x= element_text(size=14,colour="black"),axis.text.y = element_text(size=14,colour="black"),axis.title.y = element_text(size=16))+
  stat_compare_means(comparisons=my_comparisons,method = "t.test",label ="p.format",bracket.size = 0.6,size=4)
  
ggsave(file=paste0(unique(tmpdata$Genes),"_plot.pdf"),width =3,height = 4)

}
