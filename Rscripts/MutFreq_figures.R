
####### check the third base transition mutation rates

#read the overview output file.
env<-read.csv("Output/Overview/5._overview.csv")
#E2freq<-env$freq.Ts

first<-mean(env$freq.Ts[seq(1, length(env$freq.Ts), 3)],na.rm=T)
second<-mean(env$freq.Ts[seq(2, length(env$freq.Ts), 3)],na.rm=T)
third<-mean(env$freq.Ts[seq(3, length(env$freq.Ts), 3)],na.rm=T)
nonsyn<-mean(env$freq.Ts[env$Type=="nonsyn"],na.rm=T)
syn<-mean(env$freq.Ts[env$Type=="syn"],na.rm=T)

first.sd<-std.error(env$freq.Ts[seq(1, length(env$freq.Ts), 3)],na.rm=T)
second.sd<-std.error(env$freq.Ts[seq(2, length(env$freq.Ts), 3)],na.rm=T)
third.sd<-std.error(env$freq.Ts[seq(3, length(env$freq.Ts), 3)],na.rm=T)
nonsyn.sd<-std.error(env$freq.Ts[env$Type=="nonsyn"],na.rm=T)
syn.sd<-std.error(env$freq.Ts[env$Type=="syn"],na.rm=T)

#create a mutation freq mean and SE data frame:
codon<-c('first','second','third','non-syn','syn')
freq<-data.frame(CodonPosition=c('1st','2nd','3rd','syn','non-syn'), MutFreq=c(first, second,third,nonsyn,syn),
                 se=c(first.sd,second.sd,third.sd,nonsyn.sd,syn.sd))

#boxplot
barx<-barplot(freq$MutFreq, names.arg = c("first", "second", "third","non-syn","syn"),ylim=c(0,0.025))

# function for error bars
error.bar <- function(x, y, upper, lower=upper, length=0.1,...){
        if(length(x) != length(y) | length(y) !=length(lower) | length(lower) != length(upper))
                stop("vectors must be same length")
        arrows(x,y+upper, x, y-lower, angle=90, code=3, length=length, ...)
}

error.bar(barx, freq$MutFreq, freq$se)

