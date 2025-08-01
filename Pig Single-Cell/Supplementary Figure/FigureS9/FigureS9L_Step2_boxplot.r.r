rm(list=ls())
library(ggpubr)

mydata = read.csv('phe_2024_5_16.csv')

prop = read.csv("LD_sn_deconvolution-snRNAreal_142_Ctproportion_celltype11_20240717.csv")
prop[1:3,1:5]
colnames(prop)[1]='ID'
prop$ID = gsub('E......','E',prop$ID)
prop$ID = gsub('E.....','E',prop$ID)

df = merge(mydata,prop,by='ID')

merge_dat=df
rownames(merge_dat)=merge_dat$ID

#Type II myofiber C1 cluster
i=33
res <- residuals(lm(log2(merge_dat[,i]+1)~as.factor(gender)+as.numeric(age)+as.factor(motherID),data=merge_dat))
outphe <- as.data.frame(res)
colnames(outphe) <- colnames(merge_dat)[i]
outphe$ID <- rownames(outphe)

dat=merge(outphe,mydata,by='ID',all=FALSE)

p <- ggboxplot(dat,x = "GROUP3",y = colnames(merge_dat)[i],color = "GROUP3",palette=c("#C82D35","#939393"),add = "jitter")+ stat_compare_means(method = "t.test") 
ggsave(paste0(colnames(merge_dat)[i],"_proportion-t.test_C124-GR18_allsample_GROUP3.pdf"),width=4,height=4,p)