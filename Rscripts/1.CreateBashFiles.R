#create the bash files to run bbmap and bwa
# read the template command text file:
cmmd<-readLines("Data/template/Bashxx.sh")
cmmd<-readLines("Data/template/Bashxx2.sh")

#choose the fastq files to be prrocessed
fq<-list.files("Data/SIV_R21/",pattern="fastq.gz$") 

#create vector of odd numbers:
n<-seq(1, by = 2, len = (length(fq)/2))
fq2<-fq[n]


samples<-read.csv("Data/SampleID_1.csv", stringsAsFactors = F)
samples<-read.csv("Data/SampleID_2.csv", stringsAsFactors = F)

for (i in 1:length(fq2)){
  #choose the paired reads fastq files
  fa1<-fq2[i]
  fa2<-gsub(pattern="R1",replace="R2",x=fa1)
  MiseqID<-paste0("R",substr(fa1,start=1,stop=2))
  sampleid<-samples$SampleID[samples$MiseqID==MiseqID]
  new<-gsub(pattern="1_S1_L001_R1_001.fastq.gz", replace=paste0(fa1),x=cmmd)
  new<-gsub(pattern="1_S1_L001_R2_001.fastq.gz", replace=paste0(fa2),x=new)
  new<-gsub(pattern="R00",replace=paste0(sampleid),x=new)
  writeLines(new, con=paste0("Bashscripts/",sampleid,".sh"))

}

