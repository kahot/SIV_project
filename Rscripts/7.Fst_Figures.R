library(ggplot2)
library(ggpubr)
library(ggthemes)
library(plotrix)
library(grid)
library(purrr)
library(zoo)
library(tidyverse)
library(gplots)

source("Rscripts/baseRscript2.R")

# read the files saved in Overview_output:
FstFiles<-list.files("Output/Fst/",pattern=".csv")

FstList<-list()
for (i in 1:length(FstFiles)){ 
        File<-read.csv(paste0("Output/Fst/",FstFiles[i]),stringsAsFactors=FALSE)
        File<-File[,-1]
        FstList[[i]]<-File
        names(FstList)[i]<-substr(paste(FstFiles[i]),start=1,stop=11)
}

#####################
FstList2<-list()
for (j in 1:length(FstList)){
        data<-FstList[[j]]
        splitted<-strsplit(data$comb, split=" ")
        for (i in 1:nrow(data)){
                data$Sample1[i]<-splitted[[i]][1]
                data$Sample2[i]<-splitted[[i]][2]
        }
        FstList2[[j]]<-data
        names(FstList2)[j]<-names(FstList)[j]
}


samples<-read.csv('Data/SampleID_combined.csv',stringsAsFactors = F)
samples<-samples[order(samples$SampleID),]
monkeys<-unique(samples$Monkey)
for (i in 1:length(FstList2)){
        dat<-FstList2[[i]]
        dat<-dat[,4:6]
        for (j in 1:nrow(dat)){
                dat$Sample.1[j]<-paste0(substr(paste(dat$Sample1[j]), start=2,stop=2),".", samples$Sample2[samples$SampleID==dat$Sample1[j]])
                dat$Sample.2[j]<-paste0(substr(paste(dat$Sample2[j]), start=2,stop=2),".",samples$Sample2[samples$SampleID==dat$Sample2[j]])
                
        }
        s<-unique(c(dat$Sample.1,dat$Sample.2))
        for (j in 1:length(s)){
                newrow<-c(0,"","",s[j],s[j])
                dat<-rbind(dat,newrow)
        }
        dat$Fst<-as.numeric(dat$Fst)
        dat$Fst<-as.numeric(format(round(dat$Fst,4),  nsmall=4))
        
        title<-paste0("Macaque #",monkeys[i])
        p1<-ggplot(data=dat, aes(Sample.1,Sample.2, fill=Fst)) +geom_tile(color="gray40")+ scale_fill_gradient2(low="white", high="blue", midpoint=0,limit=c(-.2,.2))+
                theme_light()+ggtitle(title)+
                theme(axis.title.x = element_blank(), axis.title.y = element_blank(),panel.grid.major = element_blank(),
                      panel.border = element_blank(),
                      panel.background = element_blank(),
                      axis.text.x=element_text(size=9, angle=45,hjust=1))+
                theme(plot.title = element_text(hjust = 0.5))+
                geom_text(aes(Sample.1, Sample.2, label = Fst), color = "black", size = 3) 
        
        ggsave(file=paste0("Output/Fst/",names(FstList2)[i],".pdf"),plot=p1,width=7, height=7, units='in',device='pdf')
}






library(gplots)
library(heatmap3)
dat<-dat[,c(1,4,5)]

mat<-with(dat,tapply(Fst,list(Sample.2,Sample.1),"[[",1)) 
heatmap.2(mat, col=colorRampPalette(c("white","blue"))(50) )
heatmap3(mat, col=colorRampPalette(c("white","blue"))(50) )





mat2<-rbind(a1=c(0,NA,NA,NA),mat)
mat3<-cbind(mat2,a5=c(NA,NA,NA,NA,0))
diag(mat3)<-0
mat3m<-melt(mat3, na.rm=T)
mat3m<-mat3m[!is.na(mat3m$value),]
ggplot(data=mat3m, aes(X1,X2, fill=value)) +geom_tile(color="white")+ scale_fill_gradient2(low="blue", high="red",mid="white",midpoint=0,limit=c(-.1,.1))+
        theme_minimal()


ggplot(data = mlet_cor, aes(X2, X1, fill = value))+
        geom_tile(color = "white")+
        scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                             midpoint = 0, limit = c(-1,1), space = "Lab", 
                             name="Pearson\nCorrelation") +
        theme_minimal()+ 
        theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                         size = 12, hjust = 1))+
        coord_fixed()

library(heatmap3)
list(dat$Sample1,dat$Sample2)
M <- matrix(,5,5) 
M[cbind(dat$Sample1,dat$Sample2)]<-dat$Fst


        
        