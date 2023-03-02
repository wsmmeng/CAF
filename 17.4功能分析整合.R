rm(list=ls())
load('ego_BHresultIC.Rdata')
load('ego_BHresultMC.Rdata')
load('ego_BHresultDC.Rdata')


###用p.adjust画功能分析的热图,这个需要整合,先DC，再IC,最后MC
pathway <- c('GO:0030198','GO:0043062','GO:0045229',"GO:0030199",
             'GO:0031960','GO:0006979',"GO:0009636",'GO:0062197',
             "GO:0007015",'GO:0010810',"GO:0032956",'GO:0006936')
DC <- ego_BPresultDC[pathway,]
p.adjustDC <- DC[,c(2,6)]
rownames(p.adjustDC) <- p.adjustDC[,1]
colnames(p.adjustDC)[2] <- 'dCAF'
p.adjustDC$dCAF <- -log10(p.adjustDC$dCAF)

IC <- ego_BPresultIC[pathway,]
p.adjustIC <- IC[,c(2,6)]
rownames(p.adjustIC) <- p.adjustIC[,1]
colnames(p.adjustIC)[2] <- 'iCAF'
p.adjustIC$iCAF <- -log10(p.adjustIC$iCAF)


MC <- ego_BPresultMC[pathway,]
p.adjustMC <- MC[,c(2,6)]
rownames(p.adjustMC) <- p.adjustMC[,1]
colnames(p.adjustMC)[2] <- 'mCAF'
p.adjustMC$mCAF <- -log10(p.adjustMC$mCAF)

p.adjust <- cbind(p.adjustDC,p.adjustIC,p.adjustMC)
p.adjust <- p.adjust[,c(2,4,6)]

library(ComplexHeatmap)
##先进行标准化
data <- p.adjust
exp <- apply(data, 2, scale)##1行2列
rownames(exp) <- rownames(data)
library(circlize)

col_fun <- colorRamp2(
  c(-2, 0, 2), 
  c("white", "white", "#CC0033")
)

Heatmap(exp,
        name = 'P-value',
        col = col_fun,
        cluster_columns = F,
        cluster_rows = F,
        border = 'black',
        row_names_gp = gpar(
          col = "black",
          fontsize = 10),
        column_names_gp = gpar(
          col = "black",
          fontsize = 15),
        column_names_rot = 45)

