\name{positionPlots}
\alias{positionPlots}
\title{
Plots to evalue chip position effects on different percentiles
}
\description{
For each probe type, and for each sample, deviations from the chip mean are computed for a given percentile. These deviations are plotted against on-chip position. 
}
\usage{
positionPlots(quantiles,designInfo,percentiles=c(25,50,75))
}
\arguments{
\item{quantiles}{
A \code{list} containing three matrices. \code{list$green}, \code{list$red} and \code{list$II} must contain respectively the matrices of percentiles obtained from a \code{betaMatrix} for the Type I Green probes, Type I Red probes and Type II probes. See \code{\link{getQuantiles}}.
}
 \item{designInfo}{
\code{designInfo} matrix returned by \code{\link{getDesignInfo}}
}
 \item{percentiles}{
Vector of percentiles to be plotted. By default, the 25th, 50th and 75th percentiles are plotted. (\code{percentiles=c(25,50,75)}).
}
}
\value{
Plots are produced and saved as pdf in the current directory.
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
backgroundInfo=getBackground(greenControlMatrix, redControlMatrix)
designInfo=getDesignInfo(sampleNames)
positionPlots(quantiles, designInfo, percentiles=c(25,50,75))
}


