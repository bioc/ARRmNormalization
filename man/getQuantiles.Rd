\name{getQuantiles}
\alias{getQuantiles}
\title{
Return the percentiles of a betaMatrix for each probe type
}
\description{
It returns the percentiles of a betaMatrix for Type I Green, Type I Red and Type II probes. If no list of probes is provided, all probes are taken into account to compute the percentiles. 
}
\usage{
getQuantiles(betaMatrix,goodProbes=NULL)
}
\arguments{
  \item{betaMatrix}{
\code{matrix} containing the Beta values. Rows are probes, columns are samples.
}
 \item{goodProbes}{
Ids of the probes to be normalized (Id. of the form "cg00000029"). 
}
}
\value{
Returns a \code{list} of three matrices of percentiles. For Type I green and Type I red probes, the corresponding matrices can be accessed by \code{$green} and \code{$red}. For Type II probes, the matrix can be accessed by \code{$II} 
}
	
\author{
Jean-Philippe Fortin <jfortin@jhsph.edu>
}
\examples{
data(greenControlMatrix)
data(redControlMatrix)
data(sampleNames)
data(betaMatrix)
quantiles=getQuantiles(betaMatrix)

}


