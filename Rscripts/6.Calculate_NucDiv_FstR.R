library(ggplot2)
library(reshape)
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
SIVFiles_overview<-list.files("Output/Overview/",pattern="overview.csv")

Overview_summary<-list()
for (i in 1:length(SIVFiles_overview)){ 
        overviews<-read.csv(paste0("Output/Overview/",SIVFiles_overview[i]),stringsAsFactors=FALSE)
        overviews<-overviews[,-1]
        Overview_summary[[i]]<-overviews
        names(Overview_summary)[i]<-substr(paste(SIVFiles_overview[i]),start=1,stop=2)
}

########### No Filtering based on mut.freq -Use Minor variant freq (based on MajNT)
Overview.1<-list()
pi<-data.frame(SampleID=names(Overview_summary))
for (i in 1:length(Overview_summary)){
        dat<-Overview_summary[[i]]
        filename<-names(Overview_summary)[i]
        low_reads<-which(dat$TotalReads<1000) 
        dat[low_reads,c(10:ncol(dat))]<-NA
        
        #just transition info
        dat<-dat[,c(1:16,25,26,33,42,45)]
        
        for (k in 1:nrow(dat)){
                if (is.na(dat$MajNt[k]))  dat$D[k]<-NA
                else{
                        ac<-dat[k, "a"]*dat[k,"c"]
                        ag<-dat[k, "a"]*dat[k,"g"]
                        at<-dat[k, "a"]*dat[k,"t"]
                        cg<-dat[k, "c"]*dat[k,"g"]
                        ct<-dat[k, "c"]*dat[k,"t"]
                        gt<-dat[k, "g"]*dat[k,"t"]
                        m<-(dat[k,"TotalReads"]^2-dat[k,"TotalReads"])/2
                        dat$D[k]<-(ac+ag+at+cg+ct+gt)/m
                }
        }
        Overview.1[[i]]<-dat
        names(Overview.1)[i]<-filename
        
        n<-nrow(dat[!is.na(dat$D),])
        pi$SampleID[i]<-filename
        pi$Pi[i]<-(sum(dat$D, na.rm=T))/n
        write.csv(dat,paste0("Output/Overview_D/",filename,"_overviewD.csv"))
}


## Plot the results:   

## Sample ID information ###
#sp1<-read.csv("Data/SampleID_1.csv", stringsAsFactors = F)
#sp2<-read.csv("Data/SampleID_2.csv", stringsAsFactors = F)
#samples<-rbind(sp1,sp2)
#samples<-samples[order(samples$SampleID),]

samples<-read.csv('Data/SampleID_combined.csv',stringsAsFactors = F)
samples<-samples[order(samples$SampleID),]
s_counts<-count(samples, vars=samples$AnimalID)

s_counts$count<-as.integer("")
s_counts$count[1]<-s_counts$n[1]
for (i in 2:nrow(s_counts)) s_counts$count[i]<-s_counts$count[i-1]+s_counts$n[i]

linebreak<-s_counts$count[1:10]+0.5

s_counts$idplace<-as.integer("")
s_counts$idplace[1]<-1
for (i in 2:nrow(s_counts)) s_counts$idplace[i]<-s_counts$n[i-1]+s_counts$idplace[i-1]

## plot with different colors
pi$SIVonly<-"N"
n<-which(samples$copies.ul=="SIV only")
pi$SIVonly[n]<-"Y"
pi$SIVonly[48]<-"Y"

animalid<-unique(samples$Monkey)

xlabel<-as.character(samples$Sample2)
xlabel[48]<-"SIV stock"
xloc<-seq(1,48, by=1)

xlabel.n<-xlabel
xlabel.y<-xlabel
xlabel.y[which(pi$SIVonly=="N")]<-""
xlabel.n[which(pi$SIVonly=="Y")]<-""


pdf(file="Output/Nuc.Div_Allsamples2.pdf", width=13,height=6)
par(mar=c(9.1,4.1,4.1,1.1))
plot(pi$Pi,xaxt="n", xlab='', ylab='Nucleotide diversity',ylim=c(0.002,0.009), type="n" )
points(x=which(pi$SIVonly=="N"), pi$Pi[pi$SIVonly=="N"], pch=16,ylim=c(0.002,0.009), col="#EE6677")
points(x=which(pi$SIVonly=="Y"), pi$Pi[pi$SIVonly=="Y"], pch=16, col="#4477AA",ylim=c(0.002,0.009))
        
#mtext(xlabel, side=1, line=0.2, at=xloc, las=2, cex=0.8)
mtext(xlabel.n, col="#EE6677",side=1, line=0.2, at=xloc, las=2, cex=0.8)
mtext(xlabel.y, col="#4477AA",side=1, line=0.2, at=xloc, las=2, cex=0.8)


abline(v=linebreak, col="gray70")
par(xpd=T)
for (i in 1:length(linebreak)){
        x<-c(linebreak[i],linebreak[i])
        lines(x=x, y=c(0.002, -0.0017), col="gray70")
}
par(xpd=F)
mtext(animalid[c(1:5,7:10)],side=1, line= 6, at=s_counts$idplace[c(1:5,7:10)], cex=0.9,adj=0)
mtext(animalid[6],side=1, line= 6, at=s_counts$idplace[6]-0.4, cex=0.9,adj=0)
legend('topleft',c('Coinfected','SIV only'), bty='n', col=c("#EE6677","#4477AA"), pch=c(19,19), cex=0.8)
mtext("Macaque ID", side =1, line=6, at=-2, cex=0.9)


dev.off()

#wk<-c("3wk","3wk","3wk","8wk","5wk","3wk","8wk")
#mtext(wk, side=1, line=2.2, at=c(1,6,11,12,16,21,22), cex=0.5, col="gray30")




##############################################################
###### Calculate Fst between the samples

names(Overview.1)
#[1] "a1" "a2" "a3" "a4" "a5" "a6" "a7" "a8" "b1" "b2" "b3" "b4" "b5" "b6" "b7" "c1" "c2" "c3" "c4" "c5" "c6" "c7" "d1"
#[24] "d2" "d3" "d4" "d5" "e1" "e2" "e3" "e4" "f1" "f2" "f3" "f4" "g1" "g2" "h1" "h2" "h3" "h4" "j1" "j2" "j3" "k1" "k2"
#[47] "k3" "SI"
animals<-unique(samples$AnimalID)
Comb_list<-list()
for (i in 1:10){
        animal<-animals[i]
        A.names<-samples$SampleID[samples$AnimalID==animal]
        A.comb<-t(combn(A.names,2))
        name<-paste0("Acomb.",animal)
        Comb_list[[i]]<-A.comb
        names(Comb_list)[i]<-name
        assign(name,A.comb)
}
        


for (j in 1:10){
        Comb<-Comb_list[[j]]
        FstDF<-data.frame(matrix(ncol=0,nrow=nrow(Comb)))
        
        for (i in 1:nrow(Comb)){
                samp1<-Comb[i,1]
                samp2<-Comb[i,2]
                df1<-Overview.1[[samp1]]
                df2<-Overview.1[[samp2]]
                fname<-paste(samp1,samp2)
                FstDF$comb[i]<-fname
                
                for (k in 1:nrow(df1)){
                        if (is.na(df1$MajNt[k])|is.na(df2$MajNt[k])){
                                df1$DT[k]<-NA
                        }
                        else{
                                ac<-(df1[k, "a"]+df2[k, "a"])*(df1[k,"c"]+df2[k,"c"])
                                ag<-(df1[k, "a"]+df2[k, "a"])*(df1[k,"g"]+df2[k,"g"])
                                at<-(df1[k, "a"]+df2[k, "a"])*(df1[k,"t"]+df2[k,"t"])
                                cg<-(df1[k, "c"]+df2[k, "c"])*(df1[k,"g"]+df2[k,"g"])
                                ct<-(df1[k, "c"]+df2[k, "c"])*(df1[k,"t"]+df2[k,"t"])
                                gt<-(df1[k, "g"]+df2[k, "g"])*(df1[k,"t"]+df2[k,"t"])
                                m<-((df1[k,"TotalReads"]+df2[k,"TotalReads"])^2-df1[k,"TotalReads"]-df2[k,"TotalReads"])/2
                                df1$DT[k]<-(ac+ag+at+cg+ct+gt)/m
                        }
                }
                n<-nrow(df1[!is.na(df1$DT),])
                FstDF$piT[i]<-sum(df1$DT, na.rm=T)/n
                FstDF$piS.ave[i]<-mean(c(pi$Pi[pi$SampleID==samp1],pi$Pi[pi$SampleID==samp2]),na.rm=T )
                FstDF$Fst[i]<-(FstDF$piT[i]- FstDF$piS.ave[i])/ FstDF$piT[i]
        }
        tname<-paste0("Fst_",names(Comb_list)[j])
        assign(tname,FstDF)
        write.csv(FstDF,paste0("Output/Fst/Fst_",names(Comb_list)[j],".csv"))
}






##############
error.bar <- function(x, y, upper, lower=upper, length=0.1,...){
        if(length(x) != length(y) | length(y) !=length(lower) | length(lower) != length(upper))
                stop("vectors must be same length")
        arrows(x,y+upper, x, y-lower, angle=90, code=3, length=length, ...)
}
