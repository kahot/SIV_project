#Script to analyse the frequency data and associate with features
library(dplyr)
source("Rscripts/baseRscript2.R")

#get the file name
SIVFiles_SeqData<-list.files("Output/SeqData/",pattern="SeqData")

Overview<-list()
for (i in 1:length(SIVFiles_SeqData)){   
        id<-substr(paste(SIVFiles_SeqData[i]),start=9,stop=10)
        print(id)
        OverviewDF<-read.csv(paste0("Output/SeqData/",SIVFiles_SeqData[i]),stringsAsFactors=FALSE)
        OverviewDF<-OverviewDF[,-1]
                
        TypeOfSite<-c() 
        TypeOfSite.mj<-c() 
        TypeOfSite.tv1<-c()
        TypeOfSite.tv2<-c()

        for (codon in 1:(nrow(OverviewDF)/3)) { #modify based on reading frame
                positions <- c(codon*3-2,codon*3-1, codon*3)  #modify based on reading frame
                WTcodon <- OverviewDF$ref[positions]  
                Majcodon<- OverviewDF$MajNt[positions]
                if (is.na(Majcodon[1])|is.na(Majcodon[2])|is.na(Majcodon[3])){ 
                        WTcodon<-c('n','n','n')
                        Majcodon<-c('n','n','n')
                        mutant1codon<-c('n','n','n')
                        mutant2codon<-c('n','n','n')
                        mutant3codon<-c('n','n','n')
                        mutant4codon<-c('n','n','n')
                        mutant5codon<-c('n','n','n')
                        mutant6codon<-c('n','n','n')
                        mutant1codon.tv1<-c('n','n','n')
                        mutant2codon.tv1<-c('n','n','n')
                        mutant3codon.tv1<-c('n','n','n')
                        mutant1codon.tv2<-c('n','n','n')
                        mutant2codon.tv2<-c('n','n','n')
                        mutant3codon.tv2<-c('n','n','n')
                        }
                
                else{                        
                        mutant1codon <- c(transition(WTcodon[1]), WTcodon[2:3])  
                        mutant2codon <- c(WTcodon[1],transition(WTcodon[2]), WTcodon[3])
                        mutant3codon <- c(WTcodon[1:2], transition(WTcodon[3]))
                        
                        #transversion mutation to 'a' or 'c'
                        mutant1codon.tv1 <- c(transv1(WTcodon[1]), WTcodon[2:3]) 
                        mutant2codon.tv1 <- c(WTcodon[1],transv1(WTcodon[2]), WTcodon[3])
                        mutant3codon.tv1 <- c(WTcodon[1:2], transv1(WTcodon[3]))
                        #transversion mutation to 'g' or 't'
                        mutant1codon.tv2 <- c(transv2(WTcodon[1]), WTcodon[2:3])  
                        mutant2codon.tv2 <- c(WTcodon[1],transv2(WTcodon[2]), WTcodon[3])
                        mutant3codon.tv2 <- c(WTcodon[1:2], transv2(WTcodon[3]))

                        mutant4codon <- c(transition(Majcodon[1]), Majcodon[2:3])  
                        mutant5codon <- c(Majcodon[1],transition(Majcodon[2]), Majcodon[3])
                        mutant6codon <- c(Majcodon[1:2], transition(Majcodon[3]))
                         }

                
                TypeOfSite<-c(TypeOfSite,typeofsitefunction(WTcodon,mutant1codon))
                TypeOfSite<-c(TypeOfSite,typeofsitefunction(WTcodon,mutant2codon))
                TypeOfSite<-c(TypeOfSite,typeofsitefunction(WTcodon,mutant3codon))
                
                TypeOfSite.mj<-c(TypeOfSite.mj,typeofsitefunction(Majcodon,mutant4codon))
                TypeOfSite.mj<-c(TypeOfSite.mj,typeofsitefunction(Majcodon,mutant5codon))
                TypeOfSite.mj<-c(TypeOfSite.mj,typeofsitefunction(Majcodon,mutant6codon))
                
                TypeOfSite.tv1<-c(TypeOfSite.tv1,typeofsitefunction(WTcodon,mutant1codon.tv1))
                TypeOfSite.tv1<-c(TypeOfSite.tv1,typeofsitefunction(WTcodon,mutant2codon.tv1))
                TypeOfSite.tv1<-c(TypeOfSite.tv1,typeofsitefunction(WTcodon,mutant3codon.tv1))

                TypeOfSite.tv2<-c(TypeOfSite.tv2,typeofsitefunction(WTcodon,mutant1codon.tv2))
                TypeOfSite.tv2<-c(TypeOfSite.tv2,typeofsitefunction(WTcodon,mutant2codon.tv2))
                TypeOfSite.tv2<-c(TypeOfSite.tv2,typeofsitefunction(WTcodon,mutant3codon.tv2))
                
                
                } 

        OverviewDF$Type<-TypeOfSite[1:length(OverviewDF$pos)]
        OverviewDF$Type.mj<-TypeOfSite.mj[1:length(OverviewDF$pos)]
        OverviewDF$Type.tv1<-TypeOfSite.tv1[1:length(OverviewDF$pos)]
        OverviewDF$Type.tv2<-TypeOfSite.tv2[1:length(OverviewDF$pos)]
        
        Overview[[i]]<-OverviewDF
        
        names(Overview)[i]<-id   
}



###############################
#Mut rates and sel coefficients from Abrams paper (for HIV)
mutrates1<-read.csv("~/programs/BachelerProject_New/Data/HIVMutRates/HIVMutRates.csv")

Overview_sum<-list()

for (i in 1:length(Overview)){

        DF<-Overview[[i]]
        id<-substr(paste(SIVFiles_SeqData[i]),start=9,stop=10)

        DF$TSmutrate[DF$ref=="a"]<-mutrates1$Probability[mutrates1$Nucleotide.substitution=="AG"]
        DF$TSmutrate[DF$ref=="c"]<-mutrates1$Probability[mutrates1$Nucleotide.substitution=="CU"]
        DF$TSmutrate[DF$ref=="g"]<-mutrates1$Probability[mutrates1$Nucleotide.substitution=="GA"]
        DF$TSmutrate[DF$ref=="t"]<-mutrates1$Probability[mutrates1$Nucleotide.substitution=="UC"]

        DF$TVSmutrate1[DF$ref=="a"]<-mutrates1$Probability[mutrates1$Nucleotide.substitution=="AC"]
        DF$TVSmutrate1[DF$ref=="c"]<-mutrates1$Probability[mutrates1$Nucleotide.substitution=="CA"]
        DF$TVSmutrate1[DF$ref=="g"]<-mutrates1$Probability[mutrates1$Nucleotide.substitution=="GC"]
        DF$TVSmutrate1[DF$ref=="t"]<-mutrates1$Probability[mutrates1$Nucleotide.substitution=="UA"]
        
        DF$TVSmutrate2[DF$ref=="a"]<-mutrates1$Probability[mutrates1$Nucleotide.substitution=="AU"]
        DF$TVSmutrate2[DF$ref=="c"]<-mutrates1$Probability[mutrates1$Nucleotide.substitution=="CG"]
        DF$TVSmutrate2[DF$ref=="g"]<-mutrates1$Probability[mutrates1$Nucleotide.substitution=="GU"]
        DF$TVSmutrate2[DF$ref=="t"]<-mutrates1$Probability[mutrates1$Nucleotide.substitution=="UG"]
        
        DF$TVSmutrate.tvs[DF$ref=="a"]<-mean(mutrates1$Probability[mutrates1$Nucleotide.substitution=="AU"],mutrates1$Probability[mutrates1$Nucleotide.substitution=="AC"])
        DF$TVSmutrate.tvs[DF$ref=="c"]<-mean(mutrates1$Probability[mutrates1$Nucleotide.substitution=="CG"],mutrates1$Probability[mutrates1$Nucleotide.substitution=="CA"])
        DF$TVSmutrate.tvs[DF$ref=="g"]<-mean(mutrates1$Probability[mutrates1$Nucleotide.substitution=="GU"],mutrates1$Probability[mutrates1$Nucleotide.substitution=="GC"])
        DF$TVSmutrate.tvs[DF$ref=="t"]<-mean(mutrates1$Probability[mutrates1$Nucleotide.substitution=="UG"],mutrates1$Probability[mutrates1$Nucleotide.substitution=="UA"])

        for (k in 1:length(DF$pos)){
                DF$EstSelCoeff[k] <- EstimatedS(DF$TSmutrate[k],DF[k,colnames(DF)=='freq.Ts.ref'])
                #DF$EstSelCoeff_hiv[k] <- EstimatedS(DF$TSmutrate.hiv[k],DF[k,colnames(DF)=='freq.Ts.ref'])
                DF$EstSelCoeff_transv[k] <- EstimatedS(DF$TVSmutrate.tvs[k],DF[k,colnames(DF)=='freq.transv.ref'])
                DF$EstSelCoeff_trans1[k] <- EstimatedS(DF$TVSmutrate1[k],DF[k,colnames(DF)=='freq.transv1.ref'])
                DF$EstSelCoeff_trans2[k] <- EstimatedS(DF$TVSmutrate2[k],DF[k,colnames(DF)=='freq.transv2.ref'])
                
                if (k%%3==1){
                        if (is.na(DF$MajNt[k])|is.na(DF$MajNt[k+1])|is.na(DF$MajNt[k+2])) { DF$MajAA[k]<-"NA"
                        DF$WTAA[k]<-"NA"
                        DF$MUTAA[k]<-"NA"
                        DF$TVS1_AA[k]<-"NA"
                        DF$TVS2_AA[k]<-"NA"}
                        else {  DF$MajAA[k] = seqinr::translate(DF$MajNt[c(k,k+1,k+2)])
                        DF$WTAA[k] = seqinr::translate(DF$ref[c(k,k+1,k+2)])
                        DF$MUTAA[k] = seqinr::translate(c(transition(DF$ref[k]),DF$ref[c(k+1,k+2)]))
                        DF$TVS1_AA[k] = seqinr::translate(c(transv1(DF$ref[k]),DF$ref[c(k+1,k+2)]))
                        DF$TVS2_AA[k] = seqinr::translate(c(transv2(DF$ref[k]),DF$ref[c(k+1,k+2)]))}
                } 
                if (k%%3==2){
                        if (is.na(DF$MajNt[k-1])|is.na(DF$MajNt[k])|is.na(DF$MajNt[k+1]))  {DF$MajAA[k]<-"NA"
                        DF$WTAA[k]<-"NA"
                        DF$MUTAA[k]<-"NA"
                        DF$TVS1_AA[k]<-"NA"
                        DF$TVS2_AA[k]<-"NA"}
                        else {  DF$MajAA[k] = seqinr::translate(DF$MajNt[c(k-1,k,k+1)])
                        DF$WTAA[k] = seqinr::translate(DF$ref[c(k-1,k,k+1)])
                        DF$MUTAA[k] = seqinr::translate(c(DF$ref[c(k-1)],transition(DF$ref[k]),DF$ref[c(k+1)]))
                        DF$TVS1_AA[k] = seqinr::translate(c(DF$ref[c(k-1)],transv1(DF$ref[k]),DF$ref[c(k+1)]))
                        DF$TVS2_AA[k] = seqinr::translate(c(DF$ref[c(k-1)],transv2(DF$ref[k]),DF$ref[c(k+1)]))}
                }
                if (k%%3==0){
                        if (is.na(DF$MajNt[k-2])|is.na(DF$MajNt[k-1])|is.na(DF$MajNt[k]))  {  DF$MajAA[k]<-"NA"
                        DF$WTAA[k]<-"NA"
                        DF$MUTAA[k]<-"NA"
                        DF$TVS1_AA[k]<-"NA"
                        DF$TVS2_AA[k]<-"NA"}
                        else {  DF$MajAA[k] = seqinr::translate(DF$MajNt[c(k-2,k-1,k)])
                        DF$WTAA[k] = seqinr::translate(DF$ref[c(k-2,k-1,k)])
                        DF$MUTAA[k] = seqinr::translate(c(DF$ref[c(k-2,k-1)],transition(DF$ref[k])))
                        DF$TVS1_AA[k] = seqinr::translate(c(DF$ref[c(k-2,k-1)],transv1(DF$ref[k])))
                        DF$TVS2_AA[k] = seqinr::translate(c(DF$ref[c(k-2,k-1)],transv2(DF$ref[k])))}
                }
                

        }
       #Add whether AA change is drastic & makes CpG
        DF$bigAAChange<-0
        DF$bigAAChange.tv1<-0
        DF$bigAAChange.tv2<-0
        DF$makesCpG <- 0
        DF$makesCpG.tvs <- 0
        DF$makesCpG.tv1 <- 0
        DF$makesCpG.tv2 <- 0
        
        for(j in 1:nrow(DF)){
                WT <- amCat(DF[j,'WTAA'])
                MUT <- amCat(DF[j,'MUTAA'])
                MUT1<-amCat(DF[j,'TVS1_AA'])
                MUT2<-amCat(DF[j,'TVS2_AA'])
                
                if (WT != MUT) DF$bigAAChange[j] <- 1
                if (WT != MUT1) DF$bigAAChange.tv1[j] <- 1
                if (WT != MUT2) DF$bigAAChange.tv2[j] <- 1
        
                trip <- DF$MajNt[c(j, j+1,j+2)]
                if (is.na(trip[1])|is.na(trip[2])|is.na(trip[3])) 
                        next
                        else{
                                if (trip[1] == "c" & trip[2] == "a" ) DF$makesCpG[j] <- 1 
                                if (trip[2] == "t" & trip[3] == "g")  DF$makesCpG[j] <- 1
                                if (trip[1] == "c" & (trip[2]=="c"|trip[2]=='t')) DF$makesCpG.tvs[j] <- 1
                                if (trip[3] == "g" & (trip[2]=="a"|trip[2]=="g")) DF$makesCpG.tvs[j] <- 1

                                if (trip[1] == "c" & (trip[2]=="c"|trip[2]=='t')) DF$makesCpG.tv2[j] <- 1                                
                                if (trip[3] == "g" & (trip[2]=="a"|trip[2]=="g")) DF$makesCpG.tv1[j] <- 1
                                
                        }
        } 
        
      
        write.csv(DF,paste0("Output/Overview/",id,"_overview.csv"))
        Overview_sum[[i]]<-DF
        names(Overview_sum)[i]<-id
        print(id)
}        

#####################################################

### Read depths for all files ###
SIVFiles_SeqData<-list.files("Output/SeqData/",pattern="SeqData")

ReadsSummary<-data.frame(SampleID=matrix(nrow=length(SIVFiles_SeqData)))
ReadsSummary$MaxDepth<-""
ReadsSummary$AveDepth<-""

for (i in 1:length(SIVFiles_SeqData)){
        id<-substr(paste(SIVFiles_SeqData[i]),start=9,stop=10)
        ReadsSummary$SampleID[i]<-id
        print(id)
        SeqData<-read.csv(paste("Output/SeqData/",SIVFiles_SeqData[i],sep=""))
        ReadsSummary$MaxDepth[i]<-max(SeqData$TotalReads,na.rm=T)
        ReadsSummary$AveDepth[i]<-mean(SeqData$TotalReads,na.rm=T)
        ReadsSummary$No.ofSites[i]<-nrow(SeqData[!is.na(SeqData$MajNt),])
        ReadsSummary$SE[i]<-std.error(SeqData$TotalReads, na.rm=T)
}

write.csv(ReadsSummary,"Output/ReadsSummary2.csv")      




