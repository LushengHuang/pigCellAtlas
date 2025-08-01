library(Biobase)
library(xbioc)
library(BisqueRNA)
library(reshape2)
library(gdata)
library(Seurat)

#single-cell data
mydata=readRDS('merge_sn_8LD_rawcounts_Ct11.rds')

dat <- ExpressionSet(assayData = data.matrix(mydata@assays$RNA@counts),
                      phenoData =  new("AnnotatedDataFrame", data = mydata@meta.data))

#Top100 signature genes of cell types calculated by cellid
markers <- read.table('8snLD_normalized_celltype11_mk_use_0717.csv',header=TRUE,sep=',')

##Bulk RNA-seq data
Bulk_dat <- read.table('143LD_Count_table_merge.txt',header=TRUE)
rownames(Bulk_dat) <- Bulk_dat$geneID
Bulk_dat <- Bulk_dat[,-1]

#Extract the corresponding samples in the snRNA-seq data from the bulk RNA-seq data
Bulk_dat1 <- Bulk_dat[,colnames(Bulk_dat) %in% c('D120E63a1R','D69E63a1RF','D107E63a1RF','D66E63a1RF','D117E63a1R','D102E63a1RF','D114E63a1R','D125E63a1R')]
Bulk1 <- ExpressionSet(assayData = data.matrix(Bulk_dat1))

##Deconvolution--Bisque
res <- BisqueRNA::ReferenceBasedDecomposition(Bulk1, dat, markers = markers, use.overlap=FALSE, subject.names="SampleID", cell.types="celltype11")
est_prop <- res$bulk.props
dat_est <- melt(est_prop)
colnames(dat_est) <- c('celltype','sampleID','estation')
dat_est$sampleID = gsub('E......','E',dat_est$sampleID)
dat_est$sampleID = gsub('E.....','E',dat_est$sampleID)

#real proportion of each sample celltype in snRNA-seq data
table <- table(mydata@meta.data$celltype11,mydata@meta.data$SampleID)
table.prop <- apply(table,2,function(x) {x/sum(x)})
dat_real <- melt(table.prop)
colnames(dat_real) <- c('celltype','sampleID','real')
dat_real$sampleID = gsub('E......','E',dat_real$sampleID)
dat_real$sampleID = gsub('E.....','E',dat_real$sampleID)

dat_prop <- merge(dat_est,dat_real,by=c('celltype','sampleID'),all=FALSE)
cor(dat_prop$estation,dat_prop$real,method = "pearson")
#[1] 0.9533706   

p1=ggplot(data=dat_prop,aes(x=real,y=estation),size=2.5)+
   geom_point(size=2.5,color='blue')+ 
   labs(title='',x='The real proportion in snRNA-seq',y='The estimated proportion in bulk RNA-seq')+
   theme_bw()+     
   theme(panel.grid =element_blank()) +     
   theme(plot.title = element_text(vjust=0.5,hjust = 0.5,size=15,family="bold",colour="black"),            
         title = element_text(size=12,face='bold'),                                                        
		 axis.text=element_text(size=12,face='bold'),                          
		 legend.background=element_rect(size=0.2,color=1))+
		 geom_smooth(method = "lm", formula = y ~ x, color = "black",se=FALSE, size = .8)
ggsave("Scatter_8LDbulk_real-estimated_cell-proportion_celltype11_0717.pdf",width=4.5,height=4.5,p1)