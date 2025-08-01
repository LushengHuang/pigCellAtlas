rm(list=ls())
library(ggplot2)

#Metabolomics data
df11=read.csv('Merge_Leucine_117samples_20250522.csv',check.names = FALSE)

#Phenotypic data
mydata = read.csv('phe_2024_5_16.csv')

df2 = merge(mydata,df11,by='ID',all=FALSE)
rownames(df2) = df2$ID
colnames(df2)[29]
#[1] "L-Leucine(Leu)"
colnames(df2)[29] = 'Leucine'

res <- residuals(lm(log2(Leucine)~as.factor(gender)+as.numeric(age)+as.factor(motherID),data=df2))
outphe1 <- as.data.frame(res)
colnames(outphe1) <- 'Leucine'
outphe1$ID <- rownames(outphe1)

#Outlier Removal Based on Mean ± 3 SD(D48E and D45E)
mean <- mean(outphe1$Leucine)
sd <- sd(outphe1$Leucine)
min <- mean-3*sd
max <- mean+3*sd
outphe1 <- outphe1[outphe1$Leucine > min & outphe1$Leucine < max,] 

#Cell type proportion
prop = read.csv("LD_sn_deconvolution-snRNAreal_142_Ctproportion_celltype11_20240717.csv")
prop[1:3,1:5]
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

cor.test <- cor.test(dff$X1,dff$Leucine)
p=as.data.frame(cor.test$p.value)  
r=as.data.frame(cor.test$estimate)		 
		 
p1=ggplot(data=dff,aes(x=X1 ,y=Leucine,size=2.5))+
   geom_point(color = '#221714',size=2.5)+
   labs(title='',x='Proportion of type II C1',y='L-leucine')+
   ggtitle(paste0(paste0('P=',p)," ",paste0('Peason r=',round(r,3))))+
   theme_bw()+     
   theme(panel.grid =element_blank()) +     
   theme(plot.title = element_text(vjust=0.5,hjust = 0.5,size=15,colour="black"),           
         title = element_text(size=12),                                                          
		 axis.text=element_text(size=12),                        
		 legend.background=element_rect(size=0.2,color=1))+
		 geom_smooth(method = "lm", formula = y ~ x, color = "#2850a2",se=FALSE, size = .8)		 
ggsave("C1_proportion-Leucine.pdf",width=4,height=4,p1)