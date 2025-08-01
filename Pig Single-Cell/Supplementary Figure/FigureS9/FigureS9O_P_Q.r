##FigureS9O(ACTN3)
library(reshape2)
library(stats)

#Bulk RNA-seq data
Bulk_dat <- read.table('TPM_normalization_142LD.txt',header=TRUE)

Bulk_dat <- Bulk_dat[rowSums(Bulk_dat) > 0, ]
data <- as.data.frame(t(Bulk_dat))
data$ID <- rownames(data)
data$ID <- gsub('E......','E',data$ID)
data$ID <- gsub('E.....','E',data$ID)
data <- data[,colnames(data) %in% c('ID','ACTN3')]

#Phenotypic data
phe <- read.csv('phe_2024_5_16.csv')

dat_merge <- merge(data,phe,by='ID',all=FALSE)
df1=subset(dat_merge,select=c(ID,ACTN3,age,gender,motherID))
rownames(df1) <- df1$ID
fit1 = as.data.frame(residuals(lm(log2(ACTN3+1)~as.factor(gender)+as.numeric(age)+as.factor(motherID),data=df1)))
fit1 <- as.data.frame(fit1)
colnames(fit1)='ACTN3'
fit1$ID <- rownames(fit1)

prop=read.csv('LD_sn_deconvolution-snRNAreal_142_Ctproportion_celltype11_20240717.csv')
prop=subset(prop,select=c(X,X1))
prop$X<- gsub('E......','E',prop$X)
prop$X <- gsub('E.....','E',prop$X)
colnames(prop)[1] <- 'ID'
dat_merge1 <- merge(prop,phe,by='ID',all=FALSE)
df2=subset(dat_merge1,select=c(ID,X1,age,gender,motherID))
rownames(df2) <- df2$ID
fit2 = as.data.frame(residuals(lm(log2(X1+1)~as.factor(gender)+as.numeric(age)+as.factor(motherID),data=df2)))
fit2 <- as.data.frame(fit2)
colnames(fit2)='X1'
fit2$ID <- rownames(fit2)

df=merge(fit1,fit2,by='ID',all=FALSE)
cor.test <- cor.test(df$ACTN3,df$X1)
p=as.data.frame(cor.test$p.value)  #6.294329e-33
r=as.data.frame(cor.test$estimate) #0.8021323

#Visualization
dat=df
p=ggplot(data=dat,aes(x= ACTN3,y= X1),size=2.5)+geom_point(size=2.5,colour='blue')+  
   labs(title='',x= paste('The expression of ACTN3'),y='Proportion of type II C1')+theme_bw()+      
   theme(panel.grid =element_blank()) +     
   theme(plot.title = element_text(vjust=0.5,hjust = 0.5,size=15,family="bold",colour="black"),            
         title = element_text(size=12,face='bold'),                                                        
		 axis.text=element_text(size=12,face='bold'),                           
		 legend.background=element_rect(size=0.2,color=1))+
		 geom_smooth(method = "lm", formula = y ~ x, color = "black",se=FALSE, size = .8)
ggsave("Scatterplot_LD-ACTN3_TypeIImyofibers(C1).pdf",width=4,height=4,p)



##FigureS9P(MYL1)
library(reshape2)
library(stats)

Bulk_dat <- read.table('TPM_normalization_142LD.txt',header=TRUE)

Bulk_dat <- Bulk_dat[rowSums(Bulk_dat) > 0, ]
data <- as.data.frame(t(Bulk_dat))
data$ID <- rownames(data)
data$ID <- gsub('E......','E',data$ID)
data$ID <- gsub('E.....','E',data$ID)
data <- data[,colnames(data) %in% c('ID','MYL1')]

#Phenotypic data
phe <- read.csv('phe_2024_5_16.csv')

dat_merge <- merge(data,phe,by='ID',all=FALSE)
df1=subset(dat_merge,select=c(ID,MYL1,age,gender,motherID))
rownames(df1) <- df1$ID
fit1 = as.data.frame(residuals(lm(log2(MYL1+1)~as.factor(gender)+as.numeric(age)+as.factor(motherID),data=df1)))
fit1 <- as.data.frame(fit1)
colnames(fit1)='MYL1'
fit1$ID <- rownames(fit1)

prop=read.csv('LD_sn_deconvolution-snRNAreal_142_Ctproportion_celltype11_20240717.csv')
prop=subset(prop,select=c(X,X1))
prop$X<- gsub('E......','E',prop$X)
prop$X <- gsub('E.....','E',prop$X)
colnames(prop)[1] <- 'ID'
dat_merge1 <- merge(prop,phe,by='ID',all=FALSE)
df2=subset(dat_merge1,select=c(ID,X1,age,gender,motherID))
rownames(df2) <- df2$ID
fit2 = as.data.frame(residuals(lm(log2(X1+1)~as.factor(gender)+as.numeric(age)+as.factor(motherID),data=df2)))
fit2 <- as.data.frame(fit2)
colnames(fit2)='X1'
fit2$ID <- rownames(fit2)

df=merge(fit1,fit2,by='ID',all=FALSE)  
cor.test <- cor.test(df$MYL1,df$X1)
p=as.data.frame(cor.test$p.value)  #4.16517e-16
r=as.data.frame(cor.test$estimate) #0.6162182

#Visualization 
dat=df
p=ggplot(data=dat,aes(x= MYL1,y= X1),size=2.5)+geom_point(size=2.5,colour='blue')+  
   labs(title='',x= paste('The expression of MYL1'),y='Proportion of type II C1')+theme_bw()+      
   theme(panel.grid =element_blank()) +     
   theme(plot.title = element_text(vjust=0.5,hjust = 0.5,size=15,family="bold",colour="black"),            
         title = element_text(size=12,face='bold'),                                                        
		 axis.text=element_text(size=12,face='bold'),                            
		 legend.background=element_rect(size=0.2,color=1))+
		 geom_smooth(method = "lm", formula = y ~ x, color = "black",se=FALSE, size = .8)
ggsave("Scatterplot_LD-MYL1_TypeIImyofibers(C1).pdf",width=4,height=4,p)




##FigureS9Q(MYO18B)
library(reshape2)
library(stats)

Bulk_dat <- read.table('TPM_normalization_142LD.txt',header=TRUE)

Bulk_dat <- Bulk_dat[rowSums(Bulk_dat) > 0, ]
data <- as.data.frame(t(Bulk_dat))
data$ID <- rownames(data)
data$ID <- gsub('E......','E',data$ID)
data$ID <- gsub('E.....','E',data$ID)
data <- data[,colnames(data) %in% c('ID','MYO18B')]

#Phenotypic data
phe <- read.csv('phe_2024_5_16.csv')

dat_merge <- merge(data,phe,by='ID',all=FALSE)
df1=subset(dat_merge,select=c(ID,MYO18B,age,gender,motherID))
rownames(df1) <- df1$ID
fit1 = as.data.frame(residuals(lm(log2(MYO18B+1)~as.factor(gender)+as.numeric(age)+as.factor(motherID),data=df1)))
fit1 <- as.data.frame(fit1)
colnames(fit1)='MYO18B'
fit1$ID <- rownames(fit1)

prop=read.csv('LD_sn_deconvolution-snRNAreal_142_Ctproportion_celltype11_20240717.csv')
prop=subset(prop,select=c(X,X1))
prop$X<- gsub('E......','E',prop$X)
prop$X <- gsub('E.....','E',prop$X)
colnames(prop)[1] <- 'ID'
dat_merge1 <- merge(prop,phe,by='ID',all=FALSE)
df2=subset(dat_merge1,select=c(ID,X1,age,gender,motherID))
rownames(df2) <- df2$ID
fit2 = as.data.frame(residuals(lm(log2(X1+1)~as.factor(gender)+as.numeric(age)+as.factor(motherID),data=df2)))
fit2 <- as.data.frame(fit2)
colnames(fit2)='X1'
fit2$ID <- rownames(fit2)

df=merge(fit1,fit2,by='ID',all=FALSE)
cor.test <- cor.test(df$MYO18B,df$X1)
p=as.data.frame(cor.test$p.value)  #5.982345e-15
r=as.data.frame(cor.test$estimate) #0.5964669

#Visualization 
dat=df
p=ggplot(data=dat,aes(x= MYO18B,y= X1),size=2.5)+geom_point(size=2.5,colour='blue')+  
   labs(title='',x= paste('The expression of MYO18B'),y='Proportion of type II C1')+theme_bw()+      
   theme(panel.grid =element_blank()) +  
   theme(plot.title = element_text(vjust=0.5,hjust = 0.5,size=15,family="bold",colour="black"),            
         title = element_text(size=12,face='bold'),                                                        
		 axis.text=element_text(size=12,face='bold'),                          
		 legend.background=element_rect(size=0.2,color=1))+
		 geom_smooth(method = "lm", formula = y ~ x, color = "black",se=FALSE, size = .8)
ggsave("Scatterplot_LD-MYO18B_TypeIImyofibers(C1).pdf",width=4,height=4,p)