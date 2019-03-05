\name{bestAICsearch}
\alias{bestAICsearch}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Searches for clustering with the best AIC
}
\description{
Determines the AIC for different numbers of clusters.
}
\usage{
bestAICsearch(dataVec, minK = 2, maxK, coverage, startseed = 100, nIterations = 40, breakOnIncrease=FALSE, verbose=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
    \item{dataVec}{A vector of mutation frequencies, or a matrix or dataframe with each row corresponding to a bulk sample.}
  \item{minK}{Minimal number of clusters}
  \item{maxK}{Maximal number of clusters}
    \item{coverage}{Average coverage of the sequencing data}
  \item{startseed}{Offset for the seed for the random number generator}
  \item{nIterations}{Number of different iterations for each parameter pair}
  \item{breakOnIncrease}{Whether to exit the search over the number of clusters whenever the AIC increases}
  \item{verbose}{Whether to report on results for each number of clusters}
}

\value{
    List of different clustering results, with AIC values, for \code{minK:maxK}.
}
\author{
  Jack Kuipers
}

\seealso{
  \code{\link{VAFclusterEM}}}
\examples{
coverage <- 1e4
kclust <- 3
freqs<-c(0.25,0.32,0.36)
nos <- c(5,3,1)
datas <- c()
for(ii in 1:length(nos)){
  datas<-c(datas,rnorm(nos[ii],mean=freqs[ii],sd=VAFclusterEM:::sdfromcoverage(freqs[ii],coverage)))
}
maxK <- 5
aicMatrix <- bestAICsearch(datas, minK = 2, maxK = maxK, coverage = coverage)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~models }
\keyword{ ~cluster }
\keyword{ ~classif }
\keyword{ ~methods }
