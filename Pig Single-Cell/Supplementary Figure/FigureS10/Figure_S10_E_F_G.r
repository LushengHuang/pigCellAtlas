
# Clear the current R environment
rm(list=ls());options(stringsAsFactors=FALSE)
# --------------------- Load required libraries ---------------------#
library(ggpubr)
library(patchwork)
library(ggsci)
#------------------ Load data ------------------#
dfb <- read.csv("input_file.csv",header=T,sep=",")
dfb$Group = factor(dfb$Group, levels=c('0mg/L','1mg/L','10mg/L','50mg/L','105mg/L'))

# --------- Visualization of RT-PCR results for genes across different groups in muscle myogenic progenitors--------------#
my_comparisons <- list( c("0mg/L","1mg/L"),c("0mg/L","10mg/L"),c("0mg/L","50mg/L"),c("0mg/L","105mg/L"))
for(i in 2:ncol(dfb)){
  p <- ggplot(dfb,aes(x=Group,y=dfb[,i],fill=Group))+
              geom_bar(stat="summary",fun=mean,width = 0.5)+
              scale_color_manual(values = pal_npg('nrc')(9))+scale_fill_manual(values = pal_npg('nrc')(9))+
              stat_summary(fun.data ='mean_sd', geom ="errorbar",colour ="black",width = 0.15,size=1,position = position_dodge(0.3))+
              #geom_jitter( size =2,alpha = 0.3,shape =20)+
          theme_bw()+  
          theme_classic()+
          theme(axis.line = element_line(linewidth = 0.7), 
                axis.ticks = element_line(linewidth = 0.7)) +
          ylab(paste0("Relative mRNA levels \n of ",colnames(dfb)[i]))+
          xlab(NULL)+
          theme(axis.text.x= element_text(size=14,colour="black",angle=45,vjust=1,hjust=1),axis.text.y = element_text(size=14,colour="black"),axis.title.y = element_text(size=16),legend.position = "none")+
          stat_compare_means(comparisons=my_comparisons,method = "t.test",label ="p.format",bracket.size = 0.6,size=4)
  p
ggsave(file=paste0(colnames(dfb)[i],"_barplot_20250531_t_test.pdf"),width =3,height = 4)
}


