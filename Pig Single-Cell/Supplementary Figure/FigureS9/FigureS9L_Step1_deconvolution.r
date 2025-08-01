rm(list=ls())
library(Biobase)
library(xbioc)
library(BisqueRNA)
library(reshape2)
library(gdata)
library(Seurat)

#single-cell data
mydata = readRDS('merge_sn_8LD_rawcounts_Ct11.rds')
dat <- ExpressionSet(assayData = data.matrix(mydata@assays$RNA@counts),
                      phenoData =  new("AnnotatedDataFrame", data = mydata@meta.data))

#Top100 signature genes of cell types calculated by cellid
markers <- read.table('8snLD_normalized_celltype11_mk_use_0717.csv',header=TRUE,sep=',')

##Bulk RNA-seq data 
Bulk_dat <- read.table('143LD_Count_table_merge.txt',header=TRUE)
rownames(Bulk_dat) <- Bulk_dat$geneID
Bulk_dat <- Bulk_dat[,-1]
Bulk_dat <- Bulk_dat[,-7]

##Deconvolution--Bisque
Bulk_dat2 <- Bulk_dat[,!colnames(Bulk_dat) %in% c('D120E63a1R','D69E63a1RF','D107E63a1RF','D66E63a1RF','D117E63a1R','D102E63a1RF','D114E63a1R','D125E63a1R')]
Bulk <- ExpressionSet(assayData = data.matrix(Bulk_dat2))
res <- BisqueRNA::ReferenceBasedDecomposition(Bulk, dat, markers = markers, use.overlap=FALSE, subject.names="SampleID", cell.types="celltype11")
Bulk_prop <- res$bulk.props
dat_prop <- as.data.frame(t(Bulk_prop))
dat_prop <- round(dat_prop,4)

#real proportion of each sample celltype
table <- table(mydata@meta.data$celltype11,mydata@meta.data$SampleID)
table.prop <- apply(table,2,function(x) {x/sum(x)})
table.prop=as.data.frame(table.prop)
table.prop <- round(table.prop,4)
table.prop$celltype=rownames(table.prop)

dat_prop2=as.data.frame(t(dat_prop))
dat_prop2$celltype=rownames(dat_prop2)

df=merge(dat_prop2,table.prop,by='celltype')
df2=df
rownames(df2)=df2$celltype
df2=df2[,!colnames(df2) %in% 'celltype']
df2=as.data.frame(t(df2))
write.csv(df2,'LD_sn_deconvolution-snRNAreal_142_Ctproportion_celltype11_20240717.csv')
