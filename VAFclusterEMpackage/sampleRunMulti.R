library(VAFclusterEM)

coverage <- 1e4 # coverage of data
startseed<-100
nIterations<-100

# create some data

set.seed(startseed)

freqs<-c(0.25,0.32,0.36)
nos <- c(5,3,1)

datas <- c()

for(ii in 1:length(nos)){
  datas<-c(datas,rnorm(nos[ii],mean=freqs[ii],sd=VAFclusterEM:::sdfromcoverage(freqs[ii],coverage)))
}

freqs<-c(0.12,0.28,0.34)

datas2 <- c()

for(ii in 1:length(nos)){
    datas2<-c(datas2,rnorm(nos[ii],mean=freqs[ii],sd=VAFclusterEM:::sdfromcoverage(freqs[ii],coverage)))
}

freqs<-c(0.30,0.32,0.33)

datas3 <- c()

for(ii in 1:length(nos)){
    datas3<-c(datas3,rnorm(nos[ii],mean=freqs[ii],sd=VAFclusterEM:::sdfromcoverage(freqs[ii],coverage)))
}

datas <- rbind(datas, datas2, datas3)

# run a single example of the clustering, with fixed number of clusters:

VAFclusterEM(datas,coverage,kclust=3,startseed = startseed, nIterations = nIterations, verbose=FALSE)

# if you start with far to many clusters, then some here end up identical or empty

VAFclusterEM(datas,coverage,kclust=8,startseed = startseed, nIterations = nIterations, verbose=FALSE)

# now try lots of different numbers of clusters

minK <- 1
maxK <- min(8,ncol(datas))

AICsearch<-bestAICsearch(dataVec = datas, minK = minK, maxK = maxK, coverage = coverage, startseed = 100, nIterations = 100, breakOnIncrease=TRUE, verbose=TRUE)

AICs<-rep(NA,length(AICsearch))

for(ii in 1:length(AICsearch)){
    AICs[ii]<-AICsearch[[ii]]$AIC
}

print(AICs)

# select the final result with the best AIC

AICsearch[[which.min(AICs)]]

# in particular the final assignment of mutations to clusters is

AICsearch[[which.min(AICs)]]$newclustermembership


