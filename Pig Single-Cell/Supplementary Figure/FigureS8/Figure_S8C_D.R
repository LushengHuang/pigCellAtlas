# Fig S8C-D
rm(list=ls())
tab <- read.csv("cm_module_score.csv")
tab$stage <- gsub("Pregnancy","Late_pregnancy",tab$stage)
tab$stage <- gsub("After_pregnancy","Post_pregnancy",tab$stage)
data <- tab

library(dplyr)
library(ggplot2)
library(ggsignif)

data$stage<-factor(data$stage,levels = c("Non_pregnancy","Late_pregnancy","Post_pregnancy"))#factor
data$ShortName<-factor(data$ShortName,levels=c("RV","LV","AP","IVS"))#factor
data$stage_abbr <- data$stage
data$stage_abbr <- gsub("Non_pregnancy","NP",data$stage_abbr)
data$stage_abbr <- gsub("Late_pregnancy","LP",data$stage_abbr)
data$stage_abbr <- gsub("Post_pregnancy","PP",data$stage_abbr)
data$stage_abbr <- factor(data$stage_abbr,levels=c("NP","LP","PP"))
head(data)

my.palette<-c("#CC6600","#99CC33","#336699")
my.comparisons<-list(c("NP","LP"),c("LP","PP"),c("NP","PP"))


p1 <- ggplot(data,aes(x=stage_abbr,y=CM.M2))+#1200*800
  geom_violin(aes(fill=stage_abbr,color=stage_abbr),notch=TRUE,outlier.alpha=0.2,alpha=1)+
  ggsignif::geom_signif(comparisons=my.comparisons,test="t.test",map_signif_level=FALSE,
                        step_increase=0.2,size=0.2,textsize=4,margin_top=-0.1,
                        tip_length=0.03,vjust=0)+
  facet_wrap(~ShortName,scales="free")+
  theme_bw()+
  scale_fill_manual(values=my.palette)+
  scale_color_manual(values=my.palette)+
  theme(panel.grid=element_blank(),
        axis.text=element_text(size=10,color="black"),
        axis.title=element_text(size=15),
        axis.line=element_line(color="black",size=1),
        strip.text=element_text(size=12,color="black"),
        legend.position="none")+
  labs(x="Stage",y="Module M3 Eigengene Score")

ggsave("CM_M3.t.test.pdf",p1,width=6,height=6)



