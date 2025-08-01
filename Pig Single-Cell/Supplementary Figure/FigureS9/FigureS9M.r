rm(list=ls())
library(ggpubr)

#Phenotypic data
mydata = read.csv('phe_2024_5_16.csv')
rownames(mydata) <- mydata$ID
res1 <- residuals(lm(log2(weight)~as.factor(gender)+as.numeric(age)+as.factor(motherID),data=mydata))
outphe1 <- as.data.frame(res1)
colnames(outphe1) <- 'Weight'
outphe1$ID <- rownames(outphe1)

#Cell type proportion
prop = read.csv("LD_sn_deconvolution-snRNAreal_142_Ctproportion_celltype11_20240717.csv")
colnames(prop)[1]='ID'
prop$ID = gsub('E......','E',prop$ID)
prop$ID = gsub('E.....','E',prop$ID)

df3 = merge(mydata,prop,by='ID',all=FALSE)
rownames(df3)=df3$ID

outphe3 = df3[,c(1,2)]
i=29
res2 <- residuals(lm(log2(df3[,i]+1)~as.factor(gender)+as.numeric(age)+as.factor(motherID),data=df3))
outphe2 <- as.data.frame(res2)
colnames(outphe2) <- colnames(df3)[i]
outphe2$ID <- rownames(outphe2)

outphe3=merge(outphe3,outphe2,by='ID',all=TRUE)
outphe3=outphe3[,-2]

dff = merge(outphe3,outphe1,by='ID',all=FALSE)

cor.test <- cor.test(dff$X1,dff$Weight)
p=as.data.frame(cor.test$p.value)  
r=as.data.frame(cor.test$estimate)
		 		 
p1=ggplot(data=dff,aes(x=X1 ,y=Weight,size=2.5))+
   geom_point(color = 'blue',size=2.5)+
   labs(title='',x= paste('Proportion of type II C1'),y='Body weight')+
   ggtitle(paste0(paste0('P=',p)," ",paste0('Peason r=',round(r,3))))+
   theme_bw()+     
   theme(panel.grid =element_blank()) +    
   theme(plot.title = element_text(vjust=0.5,hjust = 0.5,size=15,colour="black"),           
         title = element_text(size=12),                                                       
		 axis.text=element_text(size=12),                   
		 legend.background=element_rect(size=0.2,color=1))+
		 geom_smooth(method = "lm", formula = y ~ x, color = "black",se=FALSE, size = .8)		 
ggsave("C1_proportion-Weight.pdf",width=4,height=4,p1)