rm(list=ls())
#------------------ Load data ------------------#
#dat <- read.csv(file="input_file.csv",header=T,sep=",")
setwd("/Users/zhangqing/Downloads/01_work/03_singlecell/01_cell_atlas/18_胎儿生长受限/GR/生长受限个体筛选/")
dat <- read.table(file="/Users/zhangqing/Downloads/01_work/03_singlecell/01_cell_atlas/18_胎儿生长受限/QTT_tau/GHR_gene_GWAS/phe2.txt",header=TRUE)

colnames(dat)
unique(dat$motherID)
#[1] "306130"   "311875"   "335129"   "333588"   "331083"   "331941"   "YL341621"
#[8] "YL340479" "YF340868" "YF340858"
#idx <- dat$motherID %in% c('YF340868','331083','333588')
#dat=dat[!idx,]
dim(dat)

dat <- dat[!is.na(dat[,"weight"]),]
rownames(dat) <- as.character(1:nrow(dat))
dat[,"motherID"] <- as.factor(dat[,"motherID"])
#------------------ Adjudst body weight for gender, age and mother‘s identify ------------------#
fit <- lm(weight~1+gender+age+motherID,data=dat)
res <- residuals(fit)
dat[names(res),"Weight_adjSow"] <- res
SD <- sd(dat[,"Weight_adjSow"],na.rm=T)
MEAN <- mean(dat[,"Weight_adjSow"],na.rm=T)
#------------------Define normal and FGR individuals among the 151 fetal pigs  ------------------#
idx <- dat[,"Weight_adjSow"] < MEAN - SD
dat[idx,"GROUP"] <- "GR"
dat[!idx,"GROUP"] <- "Control"
    
#------------------Define normal and FGR individuals among fetal pigs within the same litter  ------------------#
dat[,"GROUP2"] <- "Control"
dat2 <- NULL
for(i in unique(dat$motherID)){
    tmp <- dat[dat$motherID==i,]
    SD <- sd(tmp[,"weight"],na.rm=T)
    MEAN <- mean(tmp[,"weight"],na.rm=T)
    idx <- tmp[,"weight"] < MEAN - SD
    tmp[idx,"GROUP2"] <- "GR"
    dat2 <- rbind(dat2,tmp)
}

#------Define normal and FGR individuals among 151 fetal pigs and across individuals within the same litter  --------#
dat2[,"GROUP3"] <- "Control"
dat2[dat2$GROUP=="GR" & dat2$GROUP2=="GR","GROUP3"] <- "GR"

#------------------ Plotting  ------------------#
ggplot(data=dat2,aes(x=motherID,y=Weight_adjSow)) + 
geom_jitter(aes(color=GROUP3), position=position_jitter(0.05)) + theme_bw() +
scale_color_manual(values=c("grey","darkred")) + 
ylab("Corrected Body weight") +theme(legend.position="top")+
theme(axis.text.x=element_text(angle=90,vjust=0.5,hjust=1,color="black"))+
theme(plot.title=element_text(size=14),
axis.title.x=element_text(size=14),
axis.title.y=element_text(size=14))
ggsave(paste0("20250722weight_3",".png"),width=4,height=4.5)

#save to files
write.csv(dat2,file="FGR_outputfile.csv")