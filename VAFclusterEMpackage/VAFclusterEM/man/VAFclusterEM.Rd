\name{VAFclusterEM}
\alias{VAFclusterEM}
\title{
Gaussian mixture model clustering for VAF data
}
\description{
Parameter estimation and cluster assignment is performed using the EM algorithm, implemented by VAFclusterEM.
}
\usage{
VAFclusterEM(dataVec, coverage, kclust, startseed=100, nIterations=40, verbose=FALSE)
}
\arguments{
    \item{dataVec}{A vector of mutation frequencies, or a matrix or dataframe with each row corresponding to a bulk sample.}
    \item{coverage}{Average coverage of the sequencing data}
    \item{kclust}{Number of clusters}
    \item{startseed}{Offset for the seed for the random number generator}
    \item{nIterations}{Maximum number of iterations}
    \item{verbose}{Whether to report on results for each iteration}
}
\value{
Returns a list with these elements:

\item{AIC}{Minimum AIC over all iterations}
\item{bestseed}{Seed with the minimum AIC}
\item{kclust}{The number of clusters}
\item{relativeweights}{The relative weights of each sample to be in all groups. Maximum value is chosen for group assignment.}
\item{newclustermembership}{Group membership for each sample}
}
\author{
Jack Kuipers
}
\seealso{
\code{\link[flexmix]{flexmix}}
}
\examples{
coverage <- 1e4
kclust <- 3
freqs<-c(0.25,0.32,0.36)
nos <- c(5,3,1)
datas <- c()
for(ii in 1:length(nos)){
  datas<-c(datas,rnorm(nos[ii],mean=freqs[ii],sd=VAFclusterEM:::sdfromcoverage(freqs[ii],coverage)))
}

VAFclusterEM(datas, coverage, kclust)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~models }
\keyword{ ~cluster }
\keyword{ ~classif }
\keyword{ ~methods }
