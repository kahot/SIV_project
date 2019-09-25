#!/bin/bash


mkdir Output/R00

#0 adapter trimming
bbduk.sh in1=Data/SIV_R21/1_S1_L001_R1_001.fastq.gz in2=Data/SIV_R21/1_S1_L001_R2_001.fastq.gz  out=Output/R00/R00_adp.trimmed.fastq ref=/Users/kahotisthammer/programs/bbmap/resources/adapters.fa ktrim=r k=23 mink=11 hdist=1 stats=Output/R00/stats30_0.txt

#1 Trim reads at both ends at Score<30
bbduk.sh in=Output/R00/R00_adp.trimmed.fastq out=Output/R00/R00_trimmed.q30.fastq qtrim=rl trimq=30 stats=Output/R00/stats30_1.txt

#2. Kmer filtering
bbduk.sh in=Output/R00/R00_trimmed.q30.fastq out=Output/R00/R00_unmatched.q30.fq ref=/Users/kahotisthammer/programs/bbmap/resources/phix174_ill.ref.fa k=31 hdist=1 stats=Output/R00/stats30_2.txt

#3.Remove reads with ave score <30 
bbduk.sh in=Output/R00/R00_unmatched.q30.fq out=Output/R00/R00_clean.q30.fq maq=30 stats=Output/R00/stats30_3.txt


#4.deduplication
clumpify.sh in=Output/R00/R00_clean.q30.fq out=Output/R00/R00_clumped.fq dedupe subs=0

#5.merge#$bbmerge.sh in=Output/R00/R00_clumped.fq out=Output/merged/R00_merged.fq outu=Output/unmerged/R00_unmerged.fq ihist=Output/R00/R00_ihist.txt


#6. Align the file using bwa to the reference 

#bwa mem -t 4 -k 15 -a SIV Output/merged/R00_merged.fq  > Output/sam/R00_BWAmapped.me.sam
#bwa mem -t 4 -k 15 -a SIV Output/unmerged/R00_unmerged.fq  > Output/sam/R00_BWAmapped.un.sam
bwa mem -t 4 -k 15 -a SIV Output/R00/R00_clumped.fq  > Output/sam/R00_BWAmapped.sam


#5. convert sam to bam
samtools view -S -b Output/sam/R00_BWAmapped.sam > Output/bam/R00_BWAmapped.bam


#6. sort the bam file
samtools sort Output/bam/R00_BWAmapped.bam -o  Output/bam/R00_BWA.sort.bam

#7. index the bam file
samtools index  Output/bam/R00_BWA.sort.bam  Output/bam/R00_BWA.sort.bam.bai

rm Output/R00/R00_trimmed.q30.fastq Output/R00/R00_unmatched.q30.fq Output/bam/R00_BWAmapped.bam

