rm(list=ls())
{library(multtest)
  if(!require(multtest))install.packages("multtest")
  if(!require(Seurat))install.packages("Seurat")
  if(!require(dplyr))install.packages("dplyr")
  if(!require(mindr))install.packages("mindr")
  if(!require(mindr))install.packages("tidyverse")}

OVfall <- readRDS('OVfall.rds')
unique(OVfall$Sample)
table(OVfall$celltype)

##饼图
library(plotrix)
library(dplyr)
library(ggsci)
mynames <-   table(OVfall$celltype) %>% names()
myratio <-  table(OVfall$celltype) %>% as.numeric()
pielabel <- paste0(mynames," (", round(myratio/sum(myratio)*100,2), "%)")

cols <-pal_npg("nrc")(10)#
col

pie(myratio, labels=pielabel,
    radius = 1.0,clockwise=T,
    main = "Celltype",col = cols)

# 绘制 3D 图，family 要设置你系统支持的中文字体库
pie3D(myratio,labels = pielabel,explode = 0.1, main = "Cell Proption",
      height = 0.3)

##甜甜圈图
# 并没有直接画甜甜圈图的R包，所以在饼图源代码的基础上改改
doughnut <- function (x, labels = names(x), edges = 200, outer.radius = 0.8,
                      inner.radius=0.6, clockwise = FALSE,
                      init.angle = if (clockwise) 90 else 0, density = NULL,
                      angle = 45, col = NULL, border = FALSE, lty = NULL,
                      main = NULL, ...)
{
  if (!is.numeric(x) || any(is.na(x) | x < 0))
    stop("'x' values must be positive.")
  if (is.null(labels))
    labels <- as.character(seq_along(x))
  else labels <- as.graphicsAnnot(labels)
  x <- c(0, cumsum(x)/sum(x))
  dx <- diff(x)
  nx <- length(dx)
  plot.new()
  pin <- par("pin")
  xlim <- ylim <- c(-1, 1)
  if (pin[1L] > pin[2L])
    xlim <- (pin[1L]/pin[2L]) * xlim
  else ylim <- (pin[2L]/pin[1L]) * ylim
  plot.window(xlim, ylim, "", asp = 1)
  if (is.null(col))
    col <- if (is.null(density))
      palette()
  else par("fg")
  col <- rep(col, length.out = nx)
  border <- rep(border, length.out = nx)
  lty <- rep(lty, length.out = nx)
  angle <- rep(angle, length.out = nx)
  density <- rep(density, length.out = nx)
  twopi <- if (clockwise)
    -2 * pi
  else 2 * pi
  t2xy <- function(t, radius) {
    t2p <- twopi * t + init.angle * pi/180
    list(x = radius * cos(t2p),
         y = radius * sin(t2p))
  }
  for (i in 1L:nx) {
    n <- max(2, floor(edges * dx[i]))
    P <- t2xy(seq.int(x[i], x[i + 1], length.out = n),
              outer.radius)
    polygon(c(P$x, 0), c(P$y, 0), density = density[i],
            angle = angle[i], border = border[i],
            col = col[i], lty = lty[i])
    Pout <- t2xy(mean(x[i + 0:1]), outer.radius)
    lab <- as.character(labels[i])
    if (!is.na(lab) && nzchar(lab)) {
      lines(c(1, 1.05) * Pout$x, c(1, 1.05) * Pout$y)
      text(1.1 * Pout$x, 1.1 * Pout$y, labels[i],
           xpd = TRUE, adj = ifelse(Pout$x < 0, 1, 0),
           ...)
    }      
    Pin <- t2xy(seq.int(0, 1, length.out = n*nx),
                inner.radius)
    polygon(Pin$x, Pin$y, density = density[i],
            angle = angle[i], border = border[i],
            col = "white", lty = lty[i])
  }
  
  title(main = main, ...)
  invisible(NULL)
}

# 绘图
df <- table(OVfall$celltype) %>% as.data.frame()
labs <- paste0(df$Var1," (", round(df$Freq/sum(df$Freq)*100,2), "%)")

p.circle <- doughnut(
  df$Freq,
  labels=labs, 
  init.angle=90,     # 设置初始角度
  col = cols , # 设置颜色 
  border="white",    # 边框颜色 
  inner.radius= 0.4, # 内环大小
  cex = 0.5)           # 字体大小
## Warning in rep(lty, length.out = nx): 'x' is NULL so the result will be NULL
## Warning in rep(density, length.out = nx): 'x' is NULL so the



##堆积柱状图，注意group这个参数
library(ggplot2)
## Warning: package 'ggplot2' was built under R version 4.1.3
cellnum <- table(OVfall$celltype,OVfall$Group)
cell.prop<-as.data.frame(prop.table(cellnum))
colnames(cell.prop)<-c("Celltype","Group","Proportion")


p.bar <- ggplot(cell.prop,aes(Group,Proportion,fill=Celltype))+
  geom_bar(stat="identity",position="fill")+
  scale_fill_manual(values=cols)+#自定义fill的颜色
  ggtitle("cell proportion")+
  theme_bw()+
  theme(axis.ticks.length=unit(0.5,'cm'))+
  guides(fill=guide_legend(title=NULL))
p.bar


##箱线图--主要用于可视化单个细胞，注意sample这个参数
cellnum <- table(OVfall$Sample,OVfall$celltype)
cellnum[1:4,1:4]
##             
##              B lymphocytes Endothelial cells Epithelial cells Fibroblasts
##   BRONCHO_11           421                 0               67           2
##   BRONCHO_58            51                 0              190          37
##   EBUS_06               52                 0              520           2
##   EBUS_10              435                 0               70           1
for (i in 1:nrow(cellnum)) {
  cellnum[i,] <- cellnum[i,]/sum(cellnum[i,])  
}
cellnum <- as.data.frame(cellnum)

library(reshape2)
colnames(cellnum) <- c('Sample','Celltype','Freq')

for(i in 1:nrow(cellnum)){
  cellnum$Group[i] <- strsplit(as.character(cellnum$Sample[i]),'_')[[1]][1]
}
unique(cellnum$Group)
## [1] "BRONCHO"  "EBUS"     "EFFUSION" "LN"       "LUNG"     "NS"
library(tidyverse)
## ── Attaching packages ─────────────────────────────────────── tidyverse 1.3.1 ──
## ✔ tibble  3.1.7     ✔ purrr   0.3.4
## ✔ tidyr   1.2.0     ✔ stringr 1.4.0
## ✔ readr   2.1.2     ✔ forcats 0.5.1
## Warning: package 'tibble' was built under R version 4.1.3
## Warning: package 'tidyr' was built under R version 4.1.2
## Warning: package 'readr' was built under R version 4.1.2
## ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
## ✖ dplyr::filter() masks stats::filter()
## ✖ dplyr::lag()    masks stats::lag()

cellnum$Group <- factor(cellnum$Group,levels = c("LUNG" ,"BRONCHO" ,"EBUS",
                                                 "EFFUSION", "LN" , "NS" ))

library(ggplot2)
library(ggsci)
library(ggsignif)
## Warning: package 'ggsignif' was built under R version 4.1.1
unique(cellnum$Group)
## [1] BRONCHO  EBUS     EFFUSION LN       LUNG     NS      
## Levels: LUNG BRONCHO EBUS EFFUSION LN NS
compaired <- list(c('LUNG','BRONCHO'),
                  c('LUNG',"EBUS"),
                  c("LUNG","EFFUSION"))

myplot<- list()
for (i in 1:length(unique(cellnum$Celltype))) {
  myplot[[i]] <-ggplot(data=cellnum[cellnum$Celltype==unique(cellnum$Celltype)[i],],
                       aes(x=Group,y=Freq,fill=Group))+
    geom_violin()+
    geom_boxplot(width=0.2,
                 position = position_dodge(0.9))+
    theme_bw()+
    theme(panel.grid = element_blank(),
          axis.text.x.bottom = element_text(angle = 45,hjust = 1,vjust = 1))+
    scale_fill_manual(values = pal_npg("nrc")(10))+facet_wrap(~Celltype)+
    # stat_compare_means(aes(group=Group))+
    geom_signif(comparisons = compaired,step_increase = 0.1,map_signif_level = F,test = 't.test')#t.test, wilcox.test
}

p.box <- Seurat::CombinePlots(myplot,ncol = 3,legend='right')
## Warning: CombinePlots is being deprecated. Plots should now be combined using
## the patchwork system.
## Warning: Computation failed in `stat_signif()`:
## missing value where TRUE/FALSE needed
## Computation failed in `stat_signif()`:
## missing value where TRUE/FALSE needed
p.box