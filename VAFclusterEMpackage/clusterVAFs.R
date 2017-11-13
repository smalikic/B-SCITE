#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

datas = read.table(args[1], sep="\t", header=TRUE)/2
header <- colnames(datas)

library(VAFclusterEM)

coverage <- 1e4 # coverage of data

minK <- 1
if(length(args)>2){
  maxK <- min(as.numeric(args[3]),length(datas)) # go up to number of clusters specified as third argument
} else {
  maxK <- length(datas) # otherwise try all possibilities!
}

AICsearch<-bestAICsearch(dataVec = datas, minK = minK, maxK = maxK, coverage = coverage, startseed = 100, nIterations = 40, breakOnIncrease=TRUE, verbose=TRUE) # but exit the loop when the AIC increases

AICs<-rep(NA,length(AICsearch))

for(ii in 1:length(AICsearch)){
    AICs[ii]<-AICsearch[[ii]]$AIC
}

#print(AICs)

# in particular the final assignment of mutations to clusters is

clusterAssignment <- AICsearch[[which.min(AICs)]]$newclustermembership
output <- as.data.frame(t(as.matrix(clusterAssignment)))
colnames(output) <- header

write.table(output, args[2], sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)

#VAFclusterEM(datas,coverage,kclust=3,startseed = 100, nIterations = 100, verbose=TRUE)

