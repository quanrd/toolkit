library(beeswarm)

args <- commandArgs(TRUE)
input <- args[1]
outname <- args[2]

data <- read.delim(input)
hapnum <- length(levels(as.factor(data[,2])))

for(i in 3:length(colnames(data))){
  trait <- colnames(data)[i]
  plotname <- paste(outname,trait,"png",sep=".")
  png(plotname,width = 50 + hapnum * 150,height = 500)
  boxplot(data[,i] ~ as.factor(data[,2]),main = trait,xlab="haplotype",ylab=trait)
  beeswarm(data[,i] ~ as.factor(data[,2]),col=1:7,pch=16,add=TRUE)
  dev.off()
}