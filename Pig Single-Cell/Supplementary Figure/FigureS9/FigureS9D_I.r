rm(list=ls())
library(ggpubr)

#Metabolomics data
df11=read.csv('Merge_Leucine_117samples_20250522.csv',check.names = FALSE)

#Phenotypic data
mydata = read.csv('phe_2024_5_16.csv')

dff=merge(df11,mydata,by='ID',all=FALSE) #117
colnames(dff)[3]='Leu'

merge_dat = dff 
rownames(merge_dat)=merge_dat$ID
res <- residuals(lm(log2(Leu)~as.factor(gender)+as.numeric(age)+as.factor(motherID),data=merge_dat))
outphe <- as.data.frame(res)
colnames(outphe) <- 'Leu'
outphe$ID <- rownames(outphe)

dat=merge(outphe,mydata,all=FALSE)

#Outlier Removal Based on Mean ± 3 SD
dat <- dat[!dat$ID %in% c("D48E","D45E"),] #115

##FigureS8B
dat1 = dat[dat$ID %in% c('D109E','D70E','D73E','D98E','D67E','D69E','D43E','D46E','D49E','D60E'),]
p <- ggboxplot(dat1,x = "GROUP3",y = "L-leucine",color = "GROUP3",palette = c("Control" = "#919191", "GR" = "#ab2d32"),add = "jitter")+ stat_compare_means(method = "t.test")
ggsave("Leu-t.test_C5-GR5_GROUP3_discovery.pdf",width=4,height=4,p)

##FigureS8E
dat = dat[!dat$ID %in% c('D109E','D70E','D73E','D98E','D67E','D69E','D43E','D46E','D49E','D60E'),]
p <- ggboxplot(dat,x = "GROUP3",y = "L-leucine",color = "GROUP3",palette = c("Control" = "#919191", "GR" = "#ab2d32"),add = "jitter")+ stat_compare_means(method = "t.test")
table(dat$GROUP3)
ggsave("Leu-t.test_C99-GR6_GROUP3.pdf",width=4,height=4,p)