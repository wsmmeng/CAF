###提取对象中的count数并合并
OVPPP <- readRDS('OVPPP.rds')
OVI <- readRDS('OVI.rds')
OVP <- readRDS('OVP.rds')
OVH <- readRDS('OVH.rds')
OV <- readRDS('OV.rds')


###OVPPP
a <- rownames(OVPPP@meta.data)
length(a)##21761
m <- grep('EOC1005',a)#3071
length(m)
n <- grep('EOC136',a)#2501
length(n)
k <- grep('EOC153',a)#1380
length(k)
z <- grep('EOC227',a)#2178
length(z)
j <- grep('EOC349',a)#1483
length(j)
o <- grep('EOC372',a)#711
length(o)
l <- grep('EOC3_pPer1',a)#3142
length(l)
t <- grep('EOC443',a)#2122
length(t)
d <- grep('EOC540',a)#1921
length(d)
e <- grep('EOC733',a)#1742
length(e)
g <- grep('EOC87',a)#1510
length(g)


###OVI
a <- rownames(OVI@meta.data)
length(a)#30025
m <- grep('EOC1005',a)#1288
length(m)
n <- grep('EOC136',a)#3602
length(n)
k <- grep('EOC153',a)#2401
length(k)
z <- grep('EOC227',a)#2015
length(z)
j <- grep('EOC349',a)#1471
length(j)
o <- grep('EOC372',a)#4671
length(o)
l <- grep('EOC3_iOme2',a)#2503
length(l)
t <- grep('EOC443',a)#4463
length(t)
d <- grep('EOC540',a)#2043
length(d)
e <- grep('EOC733',a)
length(e)#4077
g <- grep('EOC87',a)
length(g)#1491

##GSE147082
##在meta.data加上一列
a <- rownames(OVP@meta.data)
length(a)#9885
m <- grep('p1',a)#1909
length(m)
n <- grep('p2',a)#1102
length(n)
k <- grep('p3',a)#3451
length(k)
z <- grep('p4',a)#1108
length(z)
x <- grep('p5',a)#1244
length(x)
l <- grep('p6',a)#1071
length(l)

##GSE158937
##在meta.data加上一列
a <- rownames(OVH@meta.data)
length(a)#15202
m <- grep('H1',a)#7123
length(m)
n <- grep('H2',a)#1533
length(n)
h <- grep('H3',a)#6546
length(h)




OV@meta.data$sample <- c(rep("EOC1005_primary_Peritoneum", 3071), 
                         rep("EOC136_primary_Mesentery",2501 ),
                         rep("EOC153_primary_Omentum",1380 ),
                         rep("EOC227_primary_Omentum", 2178),
                         rep("EOC349_primary_Peritoneum", 1483),
                         rep("EOC372_primary_Peritoneum", 711),
                         rep("EOC3_primary_Peritoneum", 3142),
                         rep("EOC443_primary_Omentum", 2122),
                         rep("EOC540_primary_Omentum", 1921),
                         rep("EOC733_primary_Peritoneum",1742 ),
                         rep("EOC87_primary_Peritoneum", 1510),
                         rep("EOC1005_interval_Tumor",1288 ),
                         rep("EOC136_interval_Omentum",3602 ),
                         rep("EOC153_interval_Omentum", 2401),
                         rep("EOC227_interval_Omentum", 2015),
                         rep("EOC349_interval_Omentum",1471 ),
                         rep("EOC372_interval_Peritoneum", 4671),
                         rep("EOC3_interval_Omentum", 2503),
                         rep("EOC443_interval_Omentum",4463 ),
                         rep("EOC540_interval_Omentum", 2043),
                         rep("EOC733_interval_Omentum", 4077),
                         rep("EOC87_interval_Omentum", 1491),
                         rep("GSE4416534", 1909),
                         rep("GSE4416535", 1102),
                         rep("GSE4416536", 3451),
                         rep("GSE4416537", 1108),
                         rep("GSE4416538", 1244),
                         rep("GSE4416539", 1071),
                         rep("GSM4816045", 7123),
                         rep("GSM4816046", 1533),
                         rep("GSM4816047", 6546)
                         
                        )#赋值条件变量
saveRDS(OV,'OV.rds')
