bestAICsearch <- function(dataVec, minK = 2, maxK, coverage, startseed = 100, nIterations = 40, breakOnIncrease=FALSE, verbose=FALSE) {
    ### Check input parameters
    if(missing(maxK)) stop("Need to input maximum number of clusters k.")
    if(missing(coverage)) stop("Need to input the coverage.")
    if (!identical(maxK - floor(maxK), 0) || maxK < 2) stop("maxK has to be a positive whole number > 1.")
    if (!identical(minK - floor(minK), 0) || minK < 1) stop("minK has to be a positive whole number > 0.")

    output <- list()

    breakFlag <- FALSE

    kk <- minK # starting value

    ### Populate AICs
    while(breakFlag==FALSE){
      if(verbose==TRUE){
        if(kk>1){
          print(paste("Now testing",kk,"clusters."))
        } else {
          print(paste("Now testing",kk,"cluster."))
        }
      }
      bestCluster <- VAFclusterEM::VAFclusterEM(dataVec = dataVec, coverage = coverage, kclust = kk, startseed = startseed, nIterations = nIterations)
      if(verbose==TRUE){
        print(paste("The AIC is",bestCluster$AIC))
      }

      output[[kk - minK + 1]] <- bestCluster

      if(breakOnIncrease==TRUE){
        if(kk==minK){ # keep track of the minimum
          minAIC <- bestCluster$AIC
        } else {
          if(bestCluster$AIC < minAIC){ # if the new number of clusters
            minAIC <- bestCluster$AIC
          } else {
            breakFlag <- TRUE # exit the loop
          }
        }
      }

      kk <- kk + 1 # increase the number of clusters

      if(kk > maxK){
        breakFlag <- TRUE # exit the loop
      }

    }

    return(output)
}
