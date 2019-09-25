library(ggplot2)
library(reshape)
library(ggpubr)
library(ggthemes)
library(plotrix)
library(grid)

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
FilteredOverview.1<-list()
for (i in 1:length(Overview_summary)){
        dat<-Overview_summary[[i]]
        filename<-names(Overview_summary)[i]
        low_reads<-which(dat$TotalReads<1000) 
        dat[low_reads,c(10:50)]<-NA
        
        #mf_high<-which(dat$freq.Ts.ref>=0.2|dat$freq.transv.ref>=0.2) # row numbers of high mut freq
        #dat[mf_high,c("freq.Ts.ref","freq.transv1.ref","freq.transv2.ref","freq.transv.ref","freq.mutations.ref")]<-NA
        #dat[mf_high,c("EstSelCoeff","EstSelCoeff_trans1","EstSelCoeff_trans2","EstSelCoeff_transv")]<-NA
        FilteredOverview.1[[i]]<-dat
        names(FilteredOverview.1)[i]<-filename
        #write.csv(dat,paste0("Output/Overview_filtered/",filename,"_overview3.csv"))
        
        
}
SIV<-FilteredOverview1[[48]]
mean(SIV$freq.Ts[SIV$Type.mj=="syn"],na.rm=T)
mean(SIV$freq.Ts[SIV$Type.mj=="nonsyn"],na.rm=T)


# combine all samples
MutFreq_all<-list()
for (i in 1:length(FilteredOverview.1)){
        dat<-FilteredOverview.1[[i]]
        filename<-names(FilteredOverview.1)[i]
        MutFreq_all[[i]]<-dat[,c("pos","freq.Ts")] 
        names(MutFreq_all)[i]<-filename
}
#assign column names for the list
for (i in 1:length(MutFreq_all)) {
        colnames(MutFreq_all[[i]])<-c("pos",paste0(names(MutFreq_all[i])))
}

AllMutFreq<-MutFreq_all %>% purrr::reduce(full_join, by='pos') #8537 sites
#AllMutFreq.maj<-MutFreq_all %>% purrr::reduce(full_join, by='pos') 
#write.csv(AllMutFreq, "Output/TsMutFreq_Summary_Ref.csv")
write.csv(AllMutFreq.maj, "Output/TsMutFreq_Summary_Maj.csv")

colMeans(AllMutFreq[,2:49],na.rm=T)
#        a1_         a2_         a3_         a4_         a5_         a6_         a7_         a8_         b1_         b2_ 
#0.002681703 0.003896694 0.001991666 0.002104200 0.002026578 0.001914841 0.001993632 0.002180669 0.002702152 0.002134920 
#b3_         b4_         b5_         b6_         b7_         c1_         c2_         c3_         c4_         c5_ 
#0.002684906 0.002150785 0.003040363 0.002254933 0.002099501 0.001868721 0.003543730 0.002031938 0.002238565 0.002157295 
#c6_         c7_         d1_         d2_         d3_         d4_         d5_         e1_         e2_         e3_ 
#0.001607114 0.001828359 0.003613817 0.003554644 0.002957782 0.003971563 0.002927706 0.002946612 0.002173432 0.002033086 
#e4_         f1_         f2_         f3_         f4_         g1_         g2_         h1_         h2_         h3_ 
#0.002345832 0.002662798 0.002226592 0.002212702 0.001820323 0.002194178 0.001699814 0.001868692 0.001884567 0.001823278 
#h4_         j1_         j2_         j3_         k1_         k2_         k3_         SI_ 
#0.002296713 0.001649991 0.002001614 0.001980859 0.001696106 0.002004724 0.001840418 0.002639172 


#look at Syn vs Nsyn
mut<-data.frame("ID"=names(FilteredOverview.1))
mut$syn<-""
mut$nonsyn<-""
for (i in 1:length(FilteredOverview.1)){
        dat<-FilteredOverview.1[[i]]
        #filename<-names(FilteredOverview.1)[i]
        mut$syn[i]<-mean(dat$freq.Ts[dat$Type.mj=="syn"],na.rm=T)
        mut$nonsyn[i]<-mean(dat$freq.Ts[dat$Type.mj=="nonsyn"],na.rm=T)
}

write.csv(mut,"Output/MutFreq/SynvsNonsyn_nofilter.csv")

y1<-dat$freq.Ts[dat$Type.mj=="syn"]
y1<-y1[!is.na(y1)]
y2<-dat$freq.Ts[dat$Type.mj=="nonsyn"]
y2<-y2[!is.na(y2)]
plot(dat$freq.Ts, type='n',ylim=c(0,0.5),ylab="Mutation frequency",xlab="")
points(x=which(dat$Type.mj=="syn"),y=y1,pch=16,col="pink",cex=0.7)
points(x=which(dat$Type.mj=="nonsyn"),y=y2,pch=16,col="lightblue",cex=0.7)
# Onw NS site has a very high mut freq (>0.25)

#plot all of them
for (i in 1:length(FilteredOverview.1)){
        dat<-FilteredOverview.1[[i]]
        dname<-names(FilteredOverview.1)[i]
        
        y1<-dat$freq.Ts[dat$Type.mj=="syn"]
        y1<-y1[!is.na(y1)]
        y2<-dat$freq.Ts[dat$Type.mj=="nonsyn"]
        y2<-y2[!is.na(y2)]
        
        pdf(file=paste0("Output/mut.freq.plots/",dname,"_plt.pdf"),height =5,width=11)
        par(mfrow=c(1,2))
        plot(dat$freq.Ts, type='n',ylim=c(0,0.5),ylab="Mutation frequency",xlab="",main=dname)
        points(x=which(dat$Type.mj=="syn"),y=y1,pch=16,col="pink",cex=0.7)
        points(x=which(dat$Type.mj=="nonsyn"),y=y2,pch=16,col="lightblue",cex=0.7)
        
        plot(dat$freq.Ts, type='n',ylim=c(0,0.1),ylab="Mutation frequency",xlab="",main=dname)
        points(x=which(dat$Type.mj=="syn"),y=y1,pch=16,col="pink",cex=0.7)
        points(x=which(dat$Type.mj=="nonsyn"),y=y2,pch=16,col="lightblue",cex=0.7)
        
        
        dev.off()
}
        
### Remove the extremes

mut.filter<-data.frame("ID"=names(FilteredOverview.1))
mut.filter$syn<-""
mut.filter$nonsyn<-""
for (i in 1:length(FilteredOverview.1)){
        dat<-FilteredOverview.1[[i]]
        #filename<-names(FilteredOverview.1)[i]
        mut.filter$syn[i]<-mean(dat$freq.Ts[dat$Type.mj=="syn" & dat$freq.Ts<0.2],na.rm=T)
        mut.filter$nonsyn[i]<-mean(dat$freq.Ts[dat$Type.mj=="nonsyn" & dat$freq.Ts<0.2],na.rm=T)
}

write.csv(mut.filter,"Output/MutFreq/SynvsNonsyn_filtered.csv")




AveMut<-data.frame(colMeans(AllMutFreq[,2:25],na.rm=T)) 




colnames(AveMut)[1]<-"Mean.Mut.Freq"
AveMut$SampleID<-rownames(AveMut)

for ( i in 2:25){
        x<-AllMutFreq[,i]
        k=i-1
        AveMut$SE[k]=std.error(x,na.rm=T)
}

id<-read.csv("Data/IDs.csv", stringsAsFactors = F)
id$no<-c(1:24)
id<-id[order(id$SampleID),c(1,2,3)]
id$SampleID <- AveMut$SampleID
AveMut_id<-merge(id, AveMut, all.x=T)
AveMut_id<-AveMut_id[order(AveMut_id$no),]


error.bar <- function(x, y, upper, lower=upper, length=0.1,...){
        if(length(x) != length(y) | length(y) !=length(lower) | length(lower) != length(upper))
                stop("vectors must be same length")
        arrows(x,y+upper, x, y-lower, angle=90, code=3, length=length, ...)
}



plot(AveMut_id$Mean.Mut.Freq, xaxt="n", xlab='', ylab='Average mutation freq.' , pch=16, ylim=c(0.0015,0.004))

xlabel<-AveMut_id$SampleID
xloc<-seq(1,24, by=1)
mtext(xlabel, side=1, line=0.8, at=xloc, las=2, cex=0.8)
abline(v=c(5.5,10.5,15.5,20.5), col="gray80")
error.bar(xloc, AveMut_id$Mean.Mut.Freq, AveMut_id$SE)

wk<-c("3wk","3wk","3wk","5wk","3wk")
mtext(wk, side=1, line=2.2, at=c(1,6,11,16,21), cex=0.7, col="gray50")

write.csv(AveMut, "Output/Ave.Ts.Mut.Freq.summary_nofilter.csv")


################  remove mut.freq >0.2)
FilteredOverview1<-list()
for (i in 1:length(Overview_summary)){
        dat<-Overview_summary[[i]]
        filename<-names(Overview_summary)[i]
        low_reads<-which(dat$TotalReads<1000) 
        dat[low_reads,c("freq.Ts.ref","freq.transv1.ref","freq.transv2.ref","freq.transv.ref","freq.mutations.ref","EstSelCoeff","EstSelCoeff_trans1","EstSelCoeff_trans2","EstSelCoeff_transv")]<-NA
        
        mf_high<-which(dat$freq.Ts.ref>=0.2|dat$freq.transv.ref>=0.2) # row numbers of high mut freq
        dat[mf_high,c("freq.Ts.ref","freq.transv1.ref","freq.transv2.ref","freq.transv.ref","freq.mutations.ref")]<-NA
        dat[mf_high,c("EstSelCoeff","EstSelCoeff_trans1","EstSelCoeff_trans2","EstSelCoeff_transv")]<-NA
        FilteredOverview1[[i]]<-dat
        names(FilteredOverview1)[i]<-filename
        #write.csv(dat,paste0("Output/Overview_filtered/",filename,"_overview3.csv"))
}

MutFreq_all<-list()
for (i in 1:length(FilteredOverview1)){
        dat<-FilteredOverview1[[i]]
        filename<-names(FilteredOverview1)[i]
        MutFreq_all[[i]]<-dat[,c("pos","freq.Ts.ref")] 
        names(MutFreq_all)[i]<-filename
}
#assign column names for the list
for (i in 1:length(MutFreq_all)) {
        colnames(MutFreq_all[[i]])<-c("pos",paste0(names(MutFreq_all[i])))
}

AllMutFreq2<-MutFreq_all %>% purrr::reduce(full_join, by='pos') #8537 sites

AveMut2<-data.frame(colMeans(AllMutFreq2[,2:25],na.rm=T))
#        R1_         R10         R11         R12         R13         R14         R15         R16         R17         R18 
#0.002838399 0.002075283 0.001922203 0.002995062 0.002032554 0.001915765 0.002236195 0.002145149 0.002055727 0.001831365 
#R19         R2_         R20         R21         R22         R23         R24         R3_         R4_         R5_ 
#0.002288820 0.002113989 0.001763534 0.002629916 0.002194542 0.001997149 0.002385449 0.002030102 0.002211760 0.002099233 
#R6_         R7_         R8_         R9_ 
#0.002490347 0.002164985 0.002329531 0.002131304 

colnames(AveMut2)[1]<-"Mean.Mut.Freq"
AveMut2$SampleID<-rownames(AveMut2)

for ( i in 2:25){
        x<-AllMutFreq2[,i]
        k=i-1
        AveMut2$SE[k]=std.error(x,na.rm=T)
}

error.bar <- function(x, y, upper, lower=upper, length=0.1,...){
        if(length(x) != length(y) | length(y) !=length(lower) | length(lower) != length(upper))
                stop("vectors must be same length")
        arrows(x,y+upper, x, y-lower, angle=90, code=3, length=length, ...)
}



plot(AveMut2$Mean.Mut.Freq, xaxt="n", xlab='', ylab='average mutation freq.' , pch=16, ylim=c(0.001,0.005))


xlabel<-AveMut2$SampleID
xloc<-seq(1,24, by=1)
mtext(xlabel, side=1, line=0.8, at=xloc, las=2, xe=0.8)

error.bar(xloc, AveMut2$Mean.Mut.Freq, AveMut2$SE)

write.csv(AveMut2, "Output/Ave.Ts.Mut.Freq.summary_filtered.csv")

id<-read.csv("Data/IDs.csv", stringsAsFactors = F)
id$no<-c(1:24)
id<-id[order(id$SampleID),c(1,2,3)]
id$SampleID <- AveMut2$SampleID
AveMut_id<-merge(id, AveMut2, all.x=T)
AveMut_id<-AveMut_id[order(AveMut_id$no),]
plot(AveMut_id$Mean.Mut.Freq, xaxt="n", xlab='', ylab='Average mutation freq.' , pch=16, ylim=c(0.001,0.005))

xlabel<-AveMut_id$SampleID
xloc<-seq(1,24, by=1)
mtext(xlabel, side=1, line=0.8, at=xloc, las=2, cex=0.8)
abline(v=c(5.5,10.5,15.5,20.5), col="gray80")
error.bar(xloc, AveMut_id$Mean.Mut.Freq, AveMut_id$SE)

wk<-c("3wk","3wk","3wk","5wk","3wk")
mtext(wk, side=1, line=2.2, at=c(1,6,11,16,21), cex=0.7, col="gray50")

reads<-read.csv("Output/ReadsSummary.csv",stringsAsFactors = F)
reads<-reads[,-1]
reads_id<-merge(id, reads, by="SampleID")
reads_id<-reads_id[order(reads_id$no),]
plot(reads_id$AveDepth, xaxt="n", xlab='', ylab='Average depth' , pch=16)
mtext(xlabel, side=1, line=0.8, at=xloc, las=2, cex=0.8)
abline(v=c(5.5,10.5,15.5,20.5), col="gray80")
error.bar(xloc, reads_id$AveDepth, reads_id$SE)



