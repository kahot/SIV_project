library(tidyverse)


source("Rscripts/baseRscript2.R")

#start of the coding region in HCV is #342 and the end is #9378 for H77. THe data don't have all the way to the end.


###########################
#
SIVFiles<-list.files("Output/CSV/",pattern="csv")

freqPatTs<-data.frame(row.names=SIVFiles)
freqPatTsRef<-data.frame(row.names=SIVFiles)
freqPatTvs<-data.frame(row.names=SIVFiles)
freqPatTvsRef<-data.frame(row.names=SIVFiles)


#coding.start<-238
#coding.end<-238+572
start<-91
end<-810
no<-data.frame("pos"=c(start:end))

reference<-read.dna("Data/SIVmac251_ENV.fasta", format = "fasta",as.character=TRUE)
ref.code<-reference[start:end]
#replace y with t  at position 726
#ref.code[636] <-"t"

#write.dna(ref.code, "Data/reference.fasta")

for (i in 1:length(SIVFiles)){
        print(i)
        id<-substr(paste(SIVFiles[i]),start=1,stop=2)
        print(id)
        SeqData<-read.csv(paste("Output/CSV/",SIVFiles[i],sep=""))
        SeqData<-SeqData[,-c(1,2,9)]
        colnames(SeqData)[1]<-"pos"
        colnames(SeqData)[7:8]<-c("deletion","insertion")
        #remove the reads<100
        reads<-SeqData$TotalReads
        
        #SeqData<-SeqData[which(SeqData$pos==coding.start):which(SeqData$pos==coding.end),]
        #count the total reads with indels
        SeqData$TotalReads_indels<-rowSums(SeqData[c("A","C","G","T","deletion","insertion")])
        colnames(SeqData)[2:5]<-c("a","c","g","t")
        
        #determine the majority nucleotide base at each site
        SeqData$MajNt<-apply(SeqData[,2:5],1,function(x) c("a","c","g","t")[which.max(x)])
        
        #Add the refrence sequence:
        SeqData<-merge(no,SeqData,by="pos",all.x=T)
        SeqData$ref<-ref.code[1:length(SeqData[,1])]
        
        #check that the right position is read  in the right reading frame
        print(seqinr::translate(SeqData$MajNt[148:178]))
        
        SeqData$transition.maj<-NA
        SeqData$transition.ref<-NA
        for (j in 1:nrow(SeqData)) SeqData$transition.maj[j]<-transition(SeqData$MajNt[j])        
        for (j in 1:nrow(SeqData)) SeqData$transition.ref[j]<-transition(SeqData$ref[j])
        
        #rearrange the columns
        #SeqData<-SeqData[,c("a","c","g","t","deletion","insertion","N","pos","TotalReads","MajNt","ref","transition.maj","transition.ref")]
     
        #determine Transition mutation freq of every site.
        for (k in 1:nrow(SeqData)){
                if (is.na(SeqData$MajNt[k])) {
                        SeqData$freq.Ts[k]<-NA #transition mutations
                        SeqData$freq.Ts.ref[k]<-NA

                        SeqData$freq.transv[k]<-NA #transversion mutations
                        SeqData$freq.transv.ref[k]<-NA
                        SeqData$freq.transv1[k]<-NA
                        SeqData$freq.transv2[k]<-NA
                        SeqData$freq.transv1.ref[k]<-NA
                        SeqData$freq.transv2.ref[k]<-NA

                        SeqData$freq.mutations.ref[k]<-NA #all mutations
                        SeqData$freq.mutations[k]<-NA
                        
                        }
                else {MajNum <- SeqData [k,paste0(SeqData$MajNt[k])]
                      MutNum1<- SeqData [k,paste0(SeqData$transition.maj[k])]
                      WTNum <- SeqData [k,paste0(SeqData$ref[k])]
                      MutNum2<- SeqData [k,paste0(SeqData$transition.ref[k])]
                      
                      SeqData$freq.Ts[k]<-MutNum1/SeqData$TotalReads[k]
                      SeqData$freq.Ts.ref[k]<-MutNum2/SeqData$TotalReads[k]
                      
                      
                      #mutation frequencies of all transversion mutataions
                      if (SeqData$MajNt[k]=="a"|SeqData$MajNt[k]=='g'){
                                TrvMutNum<-SeqData[k,"c"]+SeqData[k,"t"]}
                      if (SeqData$MajNt[k]=="c"|SeqData$MajNt[k]=="t"){
                                TrvMutNum<-SeqData[k,"a"]+SeqData[k,"g"]}
                      SeqData$freq.transv[k]<-TrvMutNum/SeqData$TotalReads[k]
                      if (SeqData$ref[k]=="a"|SeqData$ref[k]=='g'){
                              TrvMutNum2<-SeqData[k,"c"]+SeqData[k,"t"]}
                      if (SeqData$ref[k]=="c"|SeqData$ref[k]=="t"){
                              TrvMutNum2<-SeqData[k,"a"]+SeqData[k,"g"]}
                      SeqData$freq.transv.ref[k]<-TrvMutNum2/SeqData$TotalReads[k]
                      
                      #Frequenceis for specific transversion mutations (1 & 2)
                      Tvs1Num<-SeqData[k,paste0(transv1(SeqData$MajNt[k]))]
                      Tvs2Num<-SeqData[k,paste0(transv2(SeqData$MajNt[k]))]
                      SeqData$freq.transv1[k]<-Tvs1Num/SeqData$TotalReads[k]
                      SeqData$freq.transv2[k]<-Tvs2Num/SeqData$TotalReads[k]
                      Tvs1rNum<-SeqData[k,paste0(transv1(SeqData$ref[k]))]
                      Tvs2rNum<-SeqData[k,paste0(transv2(SeqData$ref[k]))]
                      SeqData$freq.transv1.ref[k]<-Tvs1Num/SeqData$TotalReads[k]
                      SeqData$freq.transv2.ref[k]<-Tvs2Num/SeqData$TotalReads[k]

                      
                      #Frequencies of all SNPs (no indels)
                      AllMutNum<-SeqData$TotalReads[k]-MajNum
                      AllMutNum2<-SeqData$TotalReads[k]-WTNum
                      
                      SeqData$freq.mutations[k]<-AllMutNum/SeqData$TotalReads[k]
                      SeqData$freq.mutations.ref[k]<-AllMutNum2/SeqData$TotalReads[k]
                      
                      }
                
                freqPatTs[i,k]<-SeqData$freq.Ts[k]
                freqPatTsRef[i,k]<-SeqData$freq.Ts.ref[k]
                freqPatTvs[i,k]<-SeqData$freq.transv[k]
                freqPatTvsRef[i,k]<-SeqData$freq.transv.ref[k]
        }
                
        write.csv(SeqData,paste0("Output/SeqData/SeqData_",id,".csv"))
}

#export the frequency summary (non-filtered) to a file as .csv
#write.csv(freqPatTs,file="Output/Ts_NonFiltered_2018-10-3.csv")
#write.csv(freqPatTsRef,file="Output/freqPatTsRef_NonFiltered_2018-10-3.csv")
#write.csv(freqPatTvs,file="Output/freqPatTvs_NonFiltered_2018-10-3.csv")


