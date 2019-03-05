VAFclusterEM <- function(dataVec, coverage, kclust, startseed=100, nIterations=40, verbose=FALSE) {

    # Check for input arguments
    if (missing(dataVec)) stop("Need a full data vector as input to cluster.")
    if(missing(coverage)) stop('Need to provide a value for the coverage.')
    if (missing(kclust)) stop('Need to provide a value for kclust.')
    if (as.integer(nIterations) < 1) {
        stop('Need to specify a positive integer.')
    } else {
        nIterations <- as.integer(nIterations)
    }
    if (startseed < 0) stop("Need to specify a positive integer as startseed.")

    tmp <- lapply(seq(nIterations), doIterate, startseed, coverage, kclust, dataVec,verbose)

    idx <- which.min(sapply(tmp, '[', 'AIC'))

    output <- tmp[[idx]]

    return(output)
}

doIterate <- function(idx, startseed, coverage, kclust, datatocluster,verbose) {

    seednumber <- startseed + idx
    set.seed(seednumber)

    if(verbose==TRUE){
      print(paste("Seed", seednumber, "with", kclust, "clusters"))
    }

    clusterresults <- VAFclusterEMcore(kclust,coverage,datatocluster)

    # finally assign to the maximal cluster
    newclustermembership <- reassignsamples(clusterresults$relativeweights)
    # could also output the relative probability of being in that cluster?

    clustersizes <- table(newclustermembership)
    kfound <- length(which(clustersizes>0))

   # these values include the pseudocounts

    totalloglike <- calcloglike(clusterresults$scoresagainstclusters,clusterresults$tauvec)

    totalAIC <- (kclust-1)-2*totalloglike
    totalBIC <- log(length(datatocluster))*(kclust-1)-2*totalloglike

    #print(totalloglike)
    #print(totalAIC)
    #print(totalBIC)

   if(verbose==TRUE){
      print(paste("Log likelihood is",totalloglike))
      print(paste("AIC is",totalAIC))
      print(paste("BIC is",totalBIC))
    }

    return(list(AIC = totalAIC, seed = seednumber, kclust = kclust, newclustermembership = newclustermembership, relativeweights=clusterresults$relativeweights))
}

# This function samples an element from a vector properly

propersample <- function(x){if(length(x)==1) x else sample(x,1)}

# this function assigns sample to the cluster with the highest weight

reassignsamples <- function(sampleprobs){
    newclustermembership <-apply(sampleprobs,1,propersample(which.max)) # find the max per row
    return(newclustermembership)
}

# this function takes in log scores and returns normalised probabilities

allrelativeprobs <- function(samplescores){
    maxscorey<-apply(samplescores,1,max) # find the max of each row
    relativeprobs<-exp(samplescores-maxscorey) # remove max for numerical stability and exponentiate
    relativeprobs<-relativeprobs/rowSums(relativeprobs) # normalise
    return(relativeprobs)
}

# this function takes in probabilities, weights them by the vector tau
# and returns normalised probabilities

allrelativeweights <- function(sampleprobs,tau){
    relativeprobs<-t(t(sampleprobs)*tau)
    relativeprobs<-relativeprobs/rowSums(relativeprobs) # normalise
    return(relativeprobs)
}

# this function takes in log scores with the weight vector and returns the log likelihood

calcloglike <- function(samplescores,tau){
    maxscorey<-apply(samplescores,1,max) # find the max of each row
    loglike<-sum(log(colSums(t(exp(samplescores-maxscorey))*tau))+maxscorey) # remove max for numerical stability and exponentiate
    return(loglike)
}

sdfromcoverage <- function(freq,coverage){
    sqrt((freq)*(1-freq)/coverage)
}

VAFclusterEMcore <- function(kclust, coverage, datatocluster){

    multiFlag <- FALSE # if we have a single sample or not

    if(is.matrix(datatocluster)||is.data.frame(datatocluster)){ # if multiple, should be a matrix
        mbig<-ncol(datatocluster)
        if(nrow(datatocluster)>1){ # with more than one row
            multiFlag <- TRUE
        }
    } else if(is.vector(datatocluster)){
        mbig<-length(datatocluster)
    }

    scoresagainstclusters<-matrix(0,mbig,kclust)

    diffystart<-10
    diffy<-diffystart # to check for convergence
    county<-0 # to check how many loops we run
    countlimit<-1e4 # hard limit on the number of loops
    errortol<-1e-10# when to stop the assignment

    clustermembership<-sample.int(kclust,mbig,replace=TRUE) # start with random groupings

    for(k in 1:kclust){
        clustersamps<-which(clustermembership==k) # find the members of the cluster
        scoresagainstclusters[clustersamps,k]<-1e-3 # increase the probabilty of clustermembership
        # to get non-uniform starting point
    }

    # this is the weights of each sample for each cluster
    relativeprobabs<-allrelativeprobs(scoresagainstclusters)
    relativeweights<-relativeprobabs # for the starting value

    # main loop
    while((diffy>0) && (county<countlimit)){

        county<-county+1 # update counter

        # first given the current weights we can update tau

        rowtots<-colSums(relativeweights)
        tauvec<-rowtots/sum(rowtots)

        # and the posterior means

        for(kk in 1:kclust){
            weightvec<-relativeweights[,kk]

            if(multiFlag){
              scoresagainstclusters[,kk] <- rep(0,mbig)
              for(jj in 1:nrow(datatocluster)){ #loop over all bulk samples

                if(sum(weightvec)>0){
                    clustermean <- sum(datatocluster[jj,]*weightvec)/sum(weightvec) # the weighted mean
                } else {
                    clustermean <- 0.95 # arbitrary number a long way from the range of 0 to 0.5 which still has a defined sd
                }
                clustersd <- sdfromcoverage(clustermean,coverage) # sd from the coverage through the binomial approximation

                scoresagainstclustersVec <- scoresagainstclusters[,kk]
                # log likelihoods of the observations
                scoresagainstclusters[,kk] <- as.numeric(scoresagainstclustersVec - (datatocluster[jj,]-clustermean)^2/(2*clustersd^2) - log(2*pi*clustersd^2))
              }

            } else {

              if(sum(weightvec)>0){
                clustermean <- sum(datatocluster*weightvec)/sum(weightvec) # the weighted mean
              } else {
                clustermean <- 0.95 # arbitrary number a long way from the range of 0 to 0.5 which still has a defined sd
              }
              clustersd <- sdfromcoverage(clustermean,coverage) # sd from the coverage through the binomial approximation

              # log likelihoods of the observations
              scoresagainstclusters[,kk]<- as.numeric(-(datatocluster-clustermean)^2/(2*clustersd^2) - log(2*pi*clustersd^2))
            }
        }

        # now we can update the weights

        relativeprobabs<-allrelativeprobs(scoresagainstclusters)
        newrelativeweights<-allrelativeweights(relativeprobabs,tauvec)

        # calculate the difference
        diffy3<-sum((newrelativeweights-relativeweights)^2)

        if(diffy3<errortol){
            diffy<-diffy-1
        } else {
            diffy<-diffystart # otherwise reset the counter
        }

        relativeweights<-newrelativeweights # for the next loop
    }

    output<-vector("list",0)
    output$relativeweights<-relativeweights
    output$relativeprobs<-relativeprobabs
    output$scoresagainstclusters<-scoresagainstclusters
    output$tauvec<-tauvec

    return(output)

}
