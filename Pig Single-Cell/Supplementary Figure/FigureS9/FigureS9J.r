library(ggplot2)

#Leucine
df11=read.csv('Merge_Leucine_117samples_20250522.csv',check.names = FALSE)

#Phenotypic data
mydata = read.csv('phe_2024_5_16.csv')

dff=merge(df11,mydata,by='ID',all=FALSE)
colnames(dff)[3]='Leu'
merge_dat = dff 
rownames(merge_dat)=merge_dat$ID
res <- residuals(lm(log2(Leu)~as.factor(gender)+as.numeric(age)+as.factor(motherID),data=merge_dat))
outphe <- as.data.frame(res)
colnames(outphe) <- 'Leu'
outphe$ID <- rownames(outphe)

rownames(mydata) <- mydata$ID
res1 <- residuals(lm(log2(weight)~as.factor(gender)+as.numeric(age)+as.factor(motherID),data=mydata))
outphe1 <- as.data.frame(res1)
colnames(outphe1) <- 'Weight'
outphe1$ID <- rownames(outphe1)

dat=merge(outphe,outphe1,all=FALSE)
##Outlier Removal Based on Mean ± 3 SD
dat <- dat[!dat$ID %in% c("D48E","D45E"),]

dat = dat[!dat$ID %in% c('D109E','D70E','D73E','D98E','D67E','D69E','D43E','D46E','D49E','D60E'),]
cor.test <- cor.test(dat$Weight,dat$Leu)
p=as.data.frame(cor.test$p.value)
r=as.data.frame(cor.test$estimate)
p1=ggplot(data=dat,aes(x=Weight ,y=Leu,size=2.5))+
   geom_point(color = '#221714',size=2.5)+
   labs(title='',x= 'Body weight',y='L-Leucine')+
   ggtitle(paste0(paste0('P=',round(p,3))," ",paste0('Peason r=',round(r,3))))+
   theme_bw()+     
   theme(panel.grid =element_blank()) +
   theme(plot.title = element_text(vjust=0.5,hjust = 0.5,size=15,colour="black"),           
         title = element_text(size=12),                                                        
		 axis.text=element_text(size=12),                     
		 legend.background=element_rect(size=0.2,color=1))+
		 geom_smooth(method = "lm", formula = y ~ x, color = "#2850a2",se=FALSE, size = .8)		 
ggsave("Cor_Weight_leucine_20250528_n=105.pdf",width=4,height=4,p1)
