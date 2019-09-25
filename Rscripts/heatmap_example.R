rm(list=ls(all=TRUE))
setwd("~/Dropbox/R/LCMS")
library(heatmap3)
#example(heatmap3)
library(gplots)
library(RColorBrewer)
library(stats)
library(ggplot2)
library(reshape2)
library(plyr)
library(ggpubr)
library(ggthemes)

data<-read.csv("RTE-foldchg.csv", header=T)
data<-data[order(data$NN),]

row.names(data)<-data$GO
data1<-data[,-1]
data.m<-data.matrix(data1)
#hmcol<-brewer.pal(5,"RdBu")
rte_heatmap <- heatmap(data.m, Rowv=TRUE, Colv=TRUE, col=redgreen(75), scale="column",margins=c(5,10))
#heatmap(data.m,col=hmcol)
heatmap.2(data.m,col =colorRampPalette(c("blue","#6668F1","white","#EE8984","red"))(100))

heatmap.2(data.m,col =colorRampPalette(c("blue","white","red"))(50),margins=c(5,16), colsep=c(1:3),
          rowsep=c(1:22),sepcolor="gray",sepwidth=c(0.01,0.01),cexCol=1, trace="none",density.info="none", keysize=1, key.xlab="Mean fold-change (log2)",
          key.title="")
heatmap3(data.m,col=colorRampPalette(c("blue","white","red"))(256),margins=c(5,18), cexRow=1,cexCol=1,scale=c("none"))

"#6668F1","white","#F2A098",
data.go<-melt(data)
ggplot(data.go, aes(variable, GO), scale="none") + geom_tile(aes(fill = value),colour = "white")+
        scale_fill_gradient(low = "blue", high = "red")

### use foldchg2.csv file
fc<-read.csv("RTE-foldchg2.csv", header=T)
fc<-fc[order(fc$NN),]

row.names(fc)<-fc$GO
fc1<-fc[,-1]
fc.m<-data.matrix(fc1)
heatmap.2(fc.m,col =colorRampPalette(c("blue","white","red"))(50),margins=c(5,16), colsep=c(1:3),
          rowsep=c(1:22),sepcolor="gray",sepwidth=c(0.01,0.01),cexCol=1, trace="none",density.info="none", keysize=1, key.xlab="Mean fold-change (log2)",
          key.title="")



##################################################
### fold changes of selected stress proteins -scatter plot
fc<-read.csv("foldchg_stressresponses.csv", header=T)
plot( FoldChg_O~FoldChg_N, data=fc,xlim=c(-1.7,1.7), ylim=c(-1.7,1.7), col=c('red','blue','green'))
abline(0,1, lty=3)

with(fc, plot(FoldChg_O ~ FoldChg_N, col=GO))
abline(0,1, lty=3)

n1<-fc$FoldChg_N    
o1<-fc$FoldChg_O
foldchg1<-data.frame(pop=c(rep("N",times=length(n1)),rep("O",times=length(o1))), log2=c(n1,o1))
boxplot(log2~pop, data=foldchg1)
mean(n1)  #-0.05238889
mean(o1)  #0.1058056

########### Fold change of all proteins (unedited)
fc2<-read.csv("foldchg_NvsO_all.csv", header=T)
plot(LogFoldCh.O ~ LogFoldCh.N, data=fc2, xlim=c(-4,4), ylim=c(-4,4), xlab="Nearshore coral protein fold change (log2)",ylab="Offshore coral protein fold change (log2)",
     pch=1,cex=0.5, col="gray30") 
abline(0,1, lty=3)

#boxplot
nf<-fc2$LogFoldCh.N    
of<-fc2$LogFoldCh.O
foldchg<-data.frame(pop=c(rep("N",times=4650),rep("O",times=4650)), log2=c(fc2$LogFoldCh.N,fc2$LogFoldCh.O))
boxplot(log2~pop, data=foldchg)

ggboxplot(foldchg,x="pop",y="log2", xlab="population", ylab="Log fold change in protein abundance",color="pop",palette=c("#00FF00","#0000FF"),show.legend=FALSE )

mean(of) #-0.04283548
mean(nf) #-0.09720538

wilcox.test(log2~pop, data=foldchg, alternative='two.sided', paired = FALSE)
#	Wilcoxon rank sum test with continuity correction

#data:  log2 by pop
#W = 16359000, p-value < 2.2e-16
#alternative hypothesis: true location shift is not equal to 0

# take the absolute change

nf2<-abs(nf)
of2<-abs(of)
foldchg2<-data.frame(pop=c(rep("N",times=4650),rep("O",times=4650)), log2=c(nf2,of2))
boxplot(log2~pop, data=foldchg2)

ggboxplot(foldchg2,x="pop",y="log2", xlab="Population", ylab="Log fold change",color="pop",palette=c("#00FF00","#0000FF"),show.legend=FALSE)

wilcox.test(of2,nf2, alternative='greater', paired = FALSE)
#data:  of2 and nf2
#W = 14912000, p-value < 2.2e-16
#alternative hypothesis: true location shift is greater than 0

mean(of2) #0.3417047
mean(nf2) #0.237866

## positive vs neagative mean
nfp<-subset(nf, nf>0)
nfn<-subset(nf, nf<0)
nfn2<-subset(nf, nf<=0)

ofp<-subset(of, of>0)
ofn<-subset(of, of<0)
ofn2<-subset(of, of<=0)

mean(nfn)   #-0.2201303
mean(nfn2)  #-0.2199438
mean(nfp)   #0.2951588
mean(ofn)   #-0.2727444
mean(ofp)   #0.5075756
mean(ofn2)   #-0.272495


#scatter plot
gs<-ggplot(foldchg,aes(x=pop,y=log2,color=pop)) + geom_point() +labs(x="Population",y="Log2 fold change in protein abundance")
gs+scale_color_manual(values=c("#00FF0080","#0000FF80"))+theme_bw()+geom_jitter(width = 0.1, height = 0.03)


#check assumption for parametric test -> don't meet. Don't use.

#Levene's Test 
library(Rcmdr)
leveneTest(log2~pop, data=foldchg2)
#Levene's Test for Homogeneity of Variance (center = median)
#Df F value    Pr(>F)    
#group    1  353.82 < 2.2e-16 ***

#"Shapiro-Wilk normality test".
shapiro.test(of2) 
#data:  nf2
#W = 0.75449, p-value < 2.2e-16

##can't use t-test


##########################################
#Volcano plot
#removed non present proteins
fc3<-read.csv("foldchg_NvsO2.csv", header=T)

with(fc3, plot(LogFoldCh.O, abs(Zscore.O), pch=4, cex=.4, col="#0000FF80", ylim=c(0,12),xlim=c(-4.2,4.2), xlab="Log2 fold change", ylab='|Zstat value|'))
with(fc3, points(LogFoldCh.N, abs(Zscore.N), pch=1, cex=.4, col="#00FF0080"))
abline(h=2, col='grey50', lty=3,lwd=1.3)

plot(LogFoldCh.O ~ LogFoldCh.N, data=fc3, xlim=c(-4,4), ylim=c(-4,4), xlab="Nearshore coral protein fold change (log2)",ylab="Offshore coral protein fold change (log2)",
     pch=1,cex=0.5, col="gray30") 
abline(0,1, lty=3)

nf3<-fc3$LogFoldCh.N
of3<-fc3$LogFoldCh.O
nf3<-abs(nf)
of3<-abs(of)
foldchg3<-data.frame(pop=c(rep("N",times=length(nf3)),rep("O",times=length(of3))), log2=c(nf3,of3))
boxplot(log2~pop, data=foldchg3)
ggboxplot(foldchg3,x="pop",y="log2", xlab="Population", ylab="Log fold change",color="pop",palette=c("#00FF00","#0000FF"),show.legend=FALSE)


mean(fc3$LogFoldCh.N, na.rm=TRUE)  ## -0.08767355
mean(fc3$LogFoldCh.O)   #-0.01240317


wilcox.test(of3,nf3, alternative='greater', paired = FALSE)
#	Wilcoxon rank sum test with continuity correction
#data:  of3 and nf3
#W = 7699600, p-value < 2.2e-16
#alternative hypothesis: true location shift is greater than 0

### Calculate the value for positive and negative values separately:
#1. Nearshore corals
nn<-fc3$LogFoldCh.N
nnp<-subset(nn,nn>0)   #1092 
mean(nnp)  #0.2992866

nnn<-subset(nn,nn<0)   #1092 
mean(nnn)  #-0.2938284
nnn2<-subset(nn,nn<=0)
mean(nnn2) #-0.2933992

#2. OFfshore corals
oo<-fc3$LogFoldCh.O
oop<-subset(oo,oo>0)   
mean(oop) #0.5145119

oon<-subset(oo,oo<0)   
mean(oon) #-0.3930198

oon2<-subset(oo,oo<=0)   
mean(oon2) #-0.39239


wilcox.test(oop,nnp, alternative='greater', paired = FALSE)
#W = 979880, p-value < 2.2e-16
wilcox.test(oon,nnn, alternative='less', paired = FALSE)
#W = 1747500, p-value = 8.8e-07
wilcox.test(oon2,nnn2, alternative='less', paired = FALSE)
#W = 1753700, p-value = 9.655e-07#


#col="#53BD3480"
#00FF0033"
#rgb(0,1,0,alpha=0.2)



#########################
### scatter plot / box plot
## eliminate the proteins that are all zero in 4 samples:
fc2<-read.csv("foldchg_NvsO_3635.csv", header=T)
plot(LogFoldCh.O ~ LogFoldCh.N, data=fc2, xlim=c(-4,4), ylim=c(-4,4), xlab="Nearshore coral protein fold change (log2)",ylab="Offshore coral protein fold change (log2)",
     pch=1,cex=0.5, col="gray30") 
abline(0,1, lty=3)

nf<-fc2$LogFoldCh.N    
of<-fc2$LogFoldCh.O
foldchg<-data.frame(pop=c(rep("N",times=3635),rep("O",times=3635)), log2=c(fc2$LogFoldCh.N,fc2$LogFoldCh.O))
boxplot(log2~pop, data=foldchg)
mean(of)# -0.02394938
mean(nf) #-0.09170124

wilcox.test(oop,nnp, alternative='greater', paired = FALSE)




#### replace the fold change to zero for 000 000 proteins
fc2<-read.csv("foldchg_NvsO_3635_replace).csv", header=T)
plot(LogFoldCh.O ~ LogFoldCh.N, data=fc2, xlim=c(-4,4), ylim=c(-4,4), xlab="Nearshore coral protein fold change (log2)",ylab="Offshore coral protein fold change (log2)",
     pch=1,cex=0.5, col="gray30") 
abline(0,1, lty=3)

#boxplot
nf<-fc2$LogFoldCh.N    
of<-fc2$LogFoldCh.O
foldchg<-data.frame(pop=c(rep("N",times=3635),rep("O",times=3635)),log2=c(fc2$LogFoldCh.N,fc2$LogFoldCh.O))
boxplot(log2~pop, data=foldchg)
mean(of)# -0.01099395
mean(nf) #-0.07587923

nf2<-abs(nf)
of2<-abs(of)
foldchg2<-data.frame(pop=c(rep("N",times=length(nf2)),rep("O",times=length(of2))), log2=c(nf2,of2))
boxplot(log2~pop, data=foldchg2)
ggboxplot(foldchg2,x="pop",y="log2", xlab="Population", ylab="Log fold change",color="pop",palette=c("#00FF00","#0000FF"),show.legend=FALSE)

wilcox.test(of2, nf2, alternative='greater', paired = FALSE)
#W = 7693900, p-value < 2.2e-16

########################


