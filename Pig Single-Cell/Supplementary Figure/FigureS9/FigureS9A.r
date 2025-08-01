rm(list=ls())
library(ggpubr)

mydata = read.csv('phe_2024_5_16.csv')
rownames(mydata) <- mydata$ID
res <- residuals(lm(log2(weight)~as.factor(gender)+as.numeric(age)+as.factor(motherID),data=mydata))
outphe <- as.data.frame(res)
colnames(outphe) <- 'Corr_weight'
outphe$ID <- rownames(outphe)

dat=merge(outphe,mydata,by='ID',all=FALSE)

p <- ggboxplot(dat,x = "GROUP3",y = "Body weight",color = "GROUP3",add = "jitter")+ stat_compare_means(method = "t.test") 
ggsave("weight-t.test_C132-GR19_GROUP3.pdf",width=4,height=4,p)