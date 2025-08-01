library(readxl)
library(reshape2)
library(ggpubr)

setwd('Analysis_results')
data=read_excel("XC02329B2PRM_Results.xlsx", sheet = "Protein_absolute_quantification")
df=data[,colnames(data) %in% c("Protein Gene","P1_1\n (nmol/μL)","P1_2\n (nmol/μL)","P1_3\n (nmol/μL)","P10_1\n (nmol/μL)",         
                               "P10_2\n (nmol/μL)","P10_3\n (nmol/μL)","P105_1\n (nmol/μL)",        
                               "P105_2\n (nmol/μL)","P105_3\n (nmol/μL)")]
colnames(df)=c("Protein Gene","P1_1","P1_2","P1_3","P10_1", 
               "P10_2","P10_3","P105_1","P105_2","P105_3")
df=as.data.frame(df)
df=melt(df)
df$Group <- sub("_.*", "", df$variable)

#Visualization
mydata = df
for (i in c('ACTN3','MYL1','MYO18B')){
df = mydata[mydata$`Protein Gene` %in% i,]
group_colors <- c("P1" = '#4DBBD5FF', "P10" = '#00A087FF', "P105" = '#F39B7FFF')
comparisons <- list(c("P1", "P10"), c("P10", "P105"), c("P1", "P105"))

p <- ggplot(df, aes(x = Group, y = value, fill = Group, color = Group)) +
  geom_boxplot(width = 0.5, outlier.shape = NA, alpha = 0.4) +  
  geom_point(size = 2, shape = 21, stroke = 1) + 
  scale_fill_manual(values = group_colors) +
  scale_color_manual(values = group_colors) +
  stat_compare_means(comparisons = comparisons, method = "t.test", label = "p.format",bracket.size = 0.6,size=4) +
  labs(title = i,x = "Group", y = paste("Absolute concentration of",i,"(nmol/μL)")) +
  theme_bw(base_size = 14) +
  theme(panel.grid = element_blank(),
        axis.text = element_text(color = "black"))+
		theme(axis.line = element_line(linewidth = 0.7), 
        axis.ticks = element_line(linewidth = 0.7),legend.position="none")+ 
  theme(axis.text.x= element_text(size=14,colour="black"),axis.text.y = element_text(size=14,colour="black"),axis.title.y = element_text(size=16))  
ggsave(paste0(i,'_Human_Compare_P1-P10-P105_PRM_ProteinAbsoluteQua.pdf'),p,width=3,height=4)
}