# Clear the current R environment
rm(list=ls());options(stringsAsFactors=FALSE)

# --------------------- Load required libraries ---------------------#
library(ggpubr)# 继承ggplot语法
library(patchwork)# 拼图包
library(ggsci)#配色包

# --------------------- Load dataset ---------------------#
dfb <- read.csv("CK_final_value.csv",header=T,sep=",")
dfb$Group = factor(dfb$Group, levels=c('0mg/L','1mg/L','10mg/L','50mg/L','105mg/L'))

# ---------- Create ggplot: visualize creatine kinase activity of human myogenic progenitors across groups ----------#
my_comparisons <- list(c("0mg/L","1mg/L"),c("0mg/L","10mg/L"),c("0mg/L","50mg/L"),c("0mg/L","105mg/L"))
p <- ggplot(dfb,aes(x=Group,y=Value,color=Group))+ # 绘制箱线图
            geom_boxplot(aes(fill=Group),alpha=0.1)+ # 设置透明度# 绘制散点
            geom_jitter()+ #设置颜色
            scale_color_manual(values = pal_npg('nrc')(9))+scale_fill_manual(values = pal_npg('nrc')(9))+# 设置主题
            theme_bw()+  #去除网格线
            theme_classic()+
            theme(axis.line = element_line(linewidth = 0.7), # 设置坐标轴线条的粗细
                  axis.ticks = element_line(linewidth = 0.7)) +# 设置坐标轴刻度线条的粗细
            ylab("Creatine Kinase activity \n OD (340nm)")+
            xlab("Leucine concentration")+
#            ggtitle("Creatine Kinase activity")+
            theme(plot.title = element_text(size=15,colour="black",hjust=0.5))+
            theme(axis.text.x= element_text(size=12,colour="black",angle=45,vjust=1,hjust=1),axis.text.y = element_text(size=12,colour="black"),axis.title.y = element_text(size=16),axis.title.x = element_text(size=16),legend.position = "none")+ stat_compare_means(comparisons = my_comparisons,
                       #label.y = c(0.292,0.335,0.39,0.40,0.415),
                       method="t.test",label ="p.format",bracket.size = 0.6,size=4)                      
p
ggsave(file=paste0('CK_plot_final_t_test.pdf'),width =3,height = 4)

