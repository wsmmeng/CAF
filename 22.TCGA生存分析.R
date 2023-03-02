rm(list=ls())
load('./数据集/OVfpkm.Rdata')
rownames(OVfinalfpkm) <- OVfinalfpkm[,1]
OVfinalfpkm <- OVfinalfpkm[,-1]
OVclinical <- OVfinalfpkm[,c(1,2)]
OVexpr <- OVfinalfpkm[,-c(1:5)]
OVexpr[1:4,1:4]
OVexpr <- as.data.frame(t(OVexpr))


marker <- gene
ingene <- intersect(marker,colnames(OVfinalfpkm))
nogene <- marker[-which(marker %in% ingene)]
nogene
marker <- marker[-which(nogene %in% marker)]
#marker[which(marker=='GSDME')]='DFNA5'
#marker[which(marker=='POLR1H')]='HTEX6'
#ingene <- intersect(marker,colnames(OVfinalfpkm))
#nogene <- marker[-which(marker %in% ingene)]
#nogene
#marker <- marker[-which(marker %in% nogene)]
setgenexpr <- OVexpr[marker,]
setgenexpr <- na.omit(setgenexpr)
set <- c(rep('DC',689))
##ssGSEA分数计算
genemarker <- data.frame(marker,set)
list<- split(as.matrix(genemarker)[,1], genemarker[,2])
library(GSVA)
ssgsea.res<-gsva(as.matrix(OVexpr),
                 list,
                 method = "ssgsea",
                 kcdf = "Gaussian",
                 abs.ranking = T)
ssgsea <- data.frame(ssgsea.res)
ssgsea <- as.data.frame(t(ssgsea ))
identical(rownames(ssgsea),rownames(OVclinical))
OS <- cbind(OVclinical,ssgsea)

##KM生存曲线
library(survminer) # 加载包
library(survival) # 加载包

##对score画生存曲线
Mscore <- median(OS$DC,na.rm = FALSE)
group <- as.vector(ifelse((OS$DC)>Mscore,"high","low"))
OSall <- cbind(OS,group)
fit <- survfit(Surv(OS.time,OS) ~ group,  # 创建生存对象 
               data = OSall) # 数据集来源
ggsurvplot(fit,pval = TRUE,title='                  Desmoplastic fibroblast',
           conf.int = TRUE,
           xlab = "Follow up time(day)",
           ylab = "survival probability",
           surv.median.line = "hv",
           palette = "aaas",
           risk.table = T)
##确实是有差异的



marker <- gene
ingene <- intersect(marker,colnames(OVfinalfpkm))
nogene <- marker[-which(marker %in% ingene)]
nogene
marker <- marker[-which(nogene %in% marker)]
#marker[which(marker=='GSDME')]='DFNA5'
#marker[which(marker=='POLR1H')]='HTEX6'
#ingene <- intersect(marker,colnames(OVfinalfpkm))
#nogene <- marker[-which(marker %in% ingene)]
#nogene
#marker <- marker[-which(marker %in% nogene)]
setgenexpr <- OVexpr[marker,]
setgenexpr <- na.omit(setgenexpr)
set <- c(rep('IC',256))
##ssGSEA分数计算
genemarker <- data.frame(marker,set)
list<- split(as.matrix(genemarker)[,1], genemarker[,2])
library(GSVA)
ssgsea.res<-gsva(as.matrix(OVexpr),
                 list,
                 method = "ssgsea",
                 kcdf = "Gaussian",
                 abs.ranking = T)
ssgsea <- data.frame(ssgsea.res)
ssgsea <- as.data.frame(t(ssgsea ))
identical(rownames(ssgsea),rownames(OVclinical))
OS <- cbind(OVclinical,ssgsea)

##KM生存曲线
library(survminer) # 加载包
library(survival) # 加载包

##对score画生存曲线
Mscore <- median(OS$IC,na.rm = FALSE)
group <- as.vector(ifelse((OS$IC)>Mscore,"high","low"))
OSall <- cbind(OS,group)
fit <- survfit(Surv(OS.time,OS) ~ group,  # 创建生存对象 
               data = OSall) # 数据集来源
ggsurvplot(fit,pval = TRUE,title='                  Desmoplastic fibroblast',
           conf.int = TRUE,
           xlab = "Follow up time(day)",
           ylab = "survival probability",
           surv.median.line = "hv",
           palette = "aaas",
           risk.table = T)
##确实是有差异的

marker<- read.csv('MCmarker.csv',header=T)
marker <- marker[,1]
ingene <- intersect(marker,colnames(OVfinalfpkm))
nogene <- marker[-which(marker %in% ingene)]
nogene
marker <- marker[-which(nogene %in% marker)]
#marker[which(marker=='GSDME')]='DFNA5'
#marker[which(marker=='POLR1H')]='HTEX6'
#ingene <- intersect(marker,colnames(OVfinalfpkm))
#nogene <- marker[-which(marker %in% ingene)]
#nogene
#marker <- marker[-which(marker %in% nogene)]
setgenexpr <- OVexpr[marker,]
setgenexpr <- na.omit(setgenexpr)
set <- c(rep('MC',47))
##ssGSEA分数计算
genemarker <- data.frame(marker,set)
list<- split(as.matrix(genemarker)[,1], genemarker[,2])
library(GSVA)
ssgsea.res<-gsva(as.matrix(OVexpr),
                 list,
                 method = "ssgsea",
                 kcdf = "Gaussian",
                 abs.ranking = T)
ssgsea <- data.frame(ssgsea.res)
ssgsea <- as.data.frame(t(ssgsea ))
identical(rownames(ssgsea),rownames(OVclinical))
OS <- cbind(OVclinical,ssgsea)

##KM生存曲线
library(survminer) # 加载包
library(survival) # 加载包

##对score画生存曲线
Mscore <- median(OS$MC,na.rm = FALSE)
group <- as.vector(ifelse((OS$MC)>Mscore,"high","low"))
OSall <- cbind(OS,group)
fit <- survfit(Surv(OS.time,OS) ~ group,  # 创建生存对象 
               data = OSall) # 数据集来源
ggsurvplot(fit,pval = TRUE,title='                  Desmoplastic fibroblast',
           conf.int = TRUE,
           xlab = "Follow up time(day)",
           ylab = "survival probability",
           surv.median.line = "hv",
           palette = "aaas",
           risk.table = T)
##确实是有差异的


marker<- read.csv('PCmarker.csv',header=T)
marker <- marker[,1]
ingene <- intersect(marker,colnames(OVfinalfpkm))
nogene <- marker[-which(marker %in% ingene)]
nogene
marker <- marker[-which(nogene %in% marker)]
#marker[which(marker=='GSDME')]='DFNA5'
#marker[which(marker=='POLR1H')]='HTEX6'
#ingene <- intersect(marker,colnames(OVfinalfpkm))
#nogene <- marker[-which(marker %in% ingene)]
#nogene
#marker <- marker[-which(marker %in% nogene)]
setgenexpr <- OVexpr[marker,]
setgenexpr <- na.omit(setgenexpr)
set <- c(rep('PC',46))
##ssGSEA分数计算
genemarker <- data.frame(marker,set)
list<- split(as.matrix(genemarker)[,1], genemarker[,2])
library(GSVA)
ssgsea.res<-gsva(as.matrix(OVexpr),
                 list,
                 method = "ssgsea",
                 kcdf = "Gaussian",
                 abs.ranking = T)
ssgsea <- data.frame(ssgsea.res)
ssgsea <- as.data.frame(t(ssgsea ))
identical(rownames(ssgsea),rownames(OVclinical))
OS <- cbind(OVclinical,ssgsea)

##KM生存曲线
library(survminer) # 加载包
library(survival) # 加载包

##对score画生存曲线
Mscore <- median(OS$PC,na.rm = FALSE)
group <- as.vector(ifelse((OS$PC)>Mscore,"high","low"))
OSall <- cbind(OS,group)
fit <- survfit(Surv(OS.time,OS) ~ group,  # 创建生存对象 
               data = OSall) # 数据集来源
ggsurvplot(fit,pval = TRUE,title='                  Desmoplastic fibroblast',
           conf.int = TRUE,
           xlab = "Follow up time(day)",
           ylab = "survival probability",
           surv.median.line = "hv",
           palette = "aaas",
           risk.table = T)
##确实是有差异的