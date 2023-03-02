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
clinicalselect <- clinical[,c(11,12,14,15,16)]
allexpr <- cbind(clinicalselect,ssgsea)
allexpr <- subset(allexpr,sensi=='resistant')