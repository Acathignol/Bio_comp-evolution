rm(list=ls())
graphics.off()

setwd(dir="/home/aweber/TP/Computational/New/Bio_comp-evolution-master/BioComp/")

rawdata = read.table(file="SimulationResults.txt",header=F,sep="|")

# rawdata = read.table(file="/home/aweber/simulation11_01_17.txt",header=F,sep="|")


plot(rawdata$V2,main="fitness",xlab="fitness",type='l',lwd=2,col="darkorange1")
# points(rawdata$V3,main="fitness",xlab="fitness",col="royalblue2",add=T,pch=15)

origin = c(250,250,250,250,250,250,250,250,250,250)
target=c(25,250,400,22,300,44,500,230,145,957)
which(rawdata$V2==max(rawdata$V2))[1]

counts <- data.frame("Origin"=origin,"Actual"=as.numeric(rawdata[which(rawdata$V2==max(rawdata$V2))[1],5:14]),"Target"= target)
counts 

barplot(as.matrix(counts))


# r = rawdata$V3/rawdata$V2
# 
# proba = (exp(exp(exp(r)))-1)/(exp(exp(exp(1)))-1)
# 
# plot(proba,ylim=c(0,1),pch=15)
# 
# plot(proba,r)

