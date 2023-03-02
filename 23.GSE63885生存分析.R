rm(list=ls())
load('./数据集/data63885.Rdata')

###表达数据
##如果用这个芯片做训练集的话就执行下面的代码
anno <- anno570
##直接写一个函数
#将过滤得到的探针名字与表达矩阵对应
matchGEO <- function(exprSet){
  tmp = rownames(exprSet)%in% anno[,1] 
  exprSet = exprSet[tmp,] 
  dim(exprSet) 
  match(rownames(exprSet),anno$ID) 
  anno = anno[match(rownames(exprSet),anno$ID),] 
  match(rownames(exprSet),anno$ID) 
  dim(exprSet) 
  dim(anno) 
  tail(sort(table(anno[,2])), n = 12L) 
  #注释到相同基因的探针，保留表达量最大的那个探针 
  { 
    MAX = by(exprSet, anno[,2],  
             function(x) rownames(x)[ which.max(rowMeans(x))]) 
    MAX = as.character(MAX) 
    exprSet = exprSet[rownames(exprSet) %in% MAX,] 
    rownames( exprSet ) = anno[ match( rownames( exprSet ), anno[ , 1 ] ), 2 ] 
  } 
  return(exprSet)
}
clinial <- data63885[,c(1:15)]
dat63885 <- data63885[,-c(1:15)]
dat63885 <- as.data.frame(t(dat63885))##去除临床数据
dat63885deal <- matchGEO(dat63885)
dat63885deal$genename <- rownames(dat63885deal)
dat63885deal <- dat63885deal[,-102]
clinical <- clinial[-which(clinial$os.time=='notknown'),]
sample <- rownames(clinical)
dat63885deal <- dat63885deal[,sample]
expr <- as.data.frame(t(dat63885deal))


sensi <- c(rep('resistant',34),rep('sensitive',41))
clinical$sensi <- sensi
all <- merge(expr,clinical,by='row.names',all = F)
rownames(all) <- all[,1]
all <- all[,-1]
OVexpr <- as.data.frame(t(all))
OVexpr2 <- apply(OVexpr,2,as.numeric)
rownames(OVexpr2) <- rownames(OVexpr)


marker <- gene
ingene <- intersect(marker,colnames(all))
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
set <- c(rep('DC',661))
##ssGSEA分数计算
genemarker <- data.frame(marker,set)
list<- split(as.matrix(genemarker)[,1], genemarker[,2])
library(GSVA)
ssgsea.res<-gsva(as.matrix(OVexpr2),
                 list,
                 method = "ssgsea",
                 kcdf = "Gaussian",
                 abs.ranking = T)
ssgsea <- data.frame(ssgsea.res)
ssgsea <- as.data.frame(t(ssgsea ))

clinicalselect <- clinical[,c(11,12,14,15,16)]
allexpr <- cbind(clinicalselect,ssgsea)


##把耐药的患者提取出来看看
#allexpr <- subset(allexpr,sensi=='resistant')
M1 <- median(allexpr[,6],na.rm = FALSE)
expr <- as.vector(ifelse((allexpr[,6])>M1,"high","low"))
M1all <- cbind(allexpr,expr)
str(M1all)
M1all$os.time <- as.numeric(M1all$os.time)
M1all$os.status.1.sur.0.dead. <- as.integer(M1all$os.status.1.sur.0.dead.)
fit <- survfit(Surv(os.time,os.status.1.sur.0.dead.) ~ expr,  # 创建生存对象 
               data = M1all) # 数据集来源

ggsurvplot(fit, # 创建的拟合对象
           data =M1all,  # 指定变量数据来源
           conf.int = F, # 显示置信区间
           pval = TRUE, # 添加P值
           title ='Desmoplastic CAFs',
           risk.table = TRUE,xlab='days') 


marker<- read.csv('ICmarker.csv',header=T)
marker <- marker[,1]
ingene <- intersect(marker,colnames(all))
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
set <- c(rep('IC',239))
##ssGSEA分数计算
genemarker <- data.frame(marker,set)
list<- split(as.matrix(genemarker)[,1], genemarker[,2])
library(GSVA)
ssgsea.res<-gsva(as.matrix(OVexpr2),
                 list,
                 method = "ssgsea",
                 kcdf = "Gaussian",
                 abs.ranking = T)
ssgsea <- data.frame(ssgsea.res)
ssgsea <- as.data.frame(t(ssgsea ))

clinicalselect <- clinical[,c(11,12,14,15,16)]
allexpr <- cbind(clinicalselect,ssgsea)


##把耐药的患者提取出来看看
#allexpr <- subset(allexpr,sensi=='resistant')
M1 <- median(allexpr[,6],na.rm = FALSE)
expr <- as.vector(ifelse((allexpr[,6])>M1,"high","low"))
M1all <- cbind(allexpr,expr)
str(M1all)
M1all$os.time <- as.numeric(M1all$os.time)
M1all$os.status.1.sur.0.dead. <- as.integer(M1all$os.status.1.sur.0.dead.)
fit <- survfit(Surv(os.time,os.status.1.sur.0.dead.) ~ expr,  # 创建生存对象 
               data = M1all) # 数据集来源

ggsurvplot(fit, # 创建的拟合对象
           data =M1all,  # 指定变量数据来源
           conf.int = F, # 显示置信区间
           pval = TRUE, # 添加P值
           title ='Inflammatory CAFs',
           risk.table = TRUE,xlab='days') 


marker<- read.csv('MCmarker.csv',header=T)
marker <- marker[,1]
ingene <- intersect(marker,colnames(all))
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
set <- c(rep('IC',42))
##ssGSEA分数计算
genemarker <- data.frame(marker,set)
list<- split(as.matrix(genemarker)[,1], genemarker[,2])
library(GSVA)
ssgsea.res<-gsva(as.matrix(OVexpr2),
                 list,
                 method = "ssgsea",
                 kcdf = "Gaussian",
                 abs.ranking = T)
ssgsea <- data.frame(ssgsea.res)
ssgsea <- as.data.frame(t(ssgsea ))

clinicalselect <- clinical[,c(11,12,14,15,16)]
allexpr <- cbind(clinicalselect,ssgsea)


##把耐药的患者提取出来看看
#allexpr <- subset(allexpr,sensi=='resistant')
M1 <- median(allexpr[,6],na.rm = FALSE)
expr <- as.vector(ifelse((allexpr[,6])>M1,"high","low"))
M1all <- cbind(allexpr,expr)
str(M1all)
M1all$os.time <- as.numeric(M1all$os.time)
M1all$os.status.1.sur.0.dead. <- as.integer(M1all$os.status.1.sur.0.dead.)
fit <- survfit(Surv(os.time,os.status.1.sur.0.dead.) ~ expr,  # 创建生存对象 
               data = M1all) # 数据集来源

ggsurvplot(fit, # 创建的拟合对象
           data =M1all,  # 指定变量数据来源
           conf.int = F, # 显示置信区间
           pval = TRUE, # 添加P值
           title ='Myofibroblast-like CAFs',
           risk.table = TRUE,xlab='days')

marker<- read.csv('NCmarker.csv',header=T)
marker <- marker[,1]
ingene <- intersect(marker,colnames(all))
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
set <- c(rep('NC',104))
##ssGSEA分数计算
genemarker <- data.frame(marker,set)
list<- split(as.matrix(genemarker)[,1], genemarker[,2])
library(GSVA)
ssgsea.res<-gsva(as.matrix(OVexpr2),
                 list,
                 method = "ssgsea",
                 kcdf = "Gaussian",
                 abs.ranking = T)
ssgsea <- data.frame(ssgsea.res)
ssgsea <- as.data.frame(t(ssgsea ))

clinicalselect <- clinical[,c(11,12,14,15,16)]
allexpr <- cbind(clinicalselect,ssgsea)


##把耐药的患者提取出来看看
#allexpr <- subset(allexpr,sensi=='resistant')
M1 <- median(allexpr[,6],na.rm = FALSE)
expr <- as.vector(ifelse((allexpr[,6])>M1,"high","low"))
M1all <- cbind(allexpr,expr)
str(M1all)
M1all$os.time <- as.numeric(M1all$os.time)
M1all$os.status.1.sur.0.dead. <- as.integer(M1all$os.status.1.sur.0.dead.)
fit <- survfit(Surv(os.time,os.status.1.sur.0.dead.) ~ expr,  # 创建生存对象 
               data = M1all) # 数据集来源

ggsurvplot(fit, # 创建的拟合对象
           data =M1all,  # 指定变量数据来源
           conf.int = F, # 显示置信区间
           pval = TRUE, # 添加P值
           title ='Normal fibroblast',
           risk.table = TRUE,xlab='days')



marker<- read.csv('PCmarker.csv',header=T)
marker <- marker[,1]
ingene <- intersect(marker,colnames(all))
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
ssgsea.res<-gsva(as.matrix(OVexpr2),
                 list,
                 method = "ssgsea",
                 kcdf = "Gaussian",
                 abs.ranking = T)
ssgsea <- data.frame(ssgsea.res)
ssgsea <- as.data.frame(t(ssgsea ))

clinicalselect <- clinical[,c(11,12,14,15,16)]
allexpr <- cbind(clinicalselect,ssgsea)


##把耐药的患者提取出来看看
#allexpr <- subset(allexpr,sensi=='resistant')
M1 <- median(allexpr[,6],na.rm = FALSE)
expr <- as.vector(ifelse((allexpr[,6])>M1,"high","low"))
M1all <- cbind(allexpr,expr)
str(M1all)
M1all$os.time <- as.numeric(M1all$os.time)
M1all$os.status.1.sur.0.dead. <- as.integer(M1all$os.status.1.sur.0.dead.)
fit <- survfit(Surv(os.time,os.status.1.sur.0.dead.) ~ expr,  # 创建生存对象 
               data = M1all) # 数据集来源

ggsurvplot(fit, # 创建的拟合对象
           data =M1all,  # 指定变量数据来源
           conf.int = F, # 显示置信区间
           pval = TRUE, # 添加P值
           title ='Proliferating CAF',
           risk.table = TRUE,xlab='days')