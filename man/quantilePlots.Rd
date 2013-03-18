\name{quantilePlots}
\alias{quantilePlots}
\title{
Diagnostic plots for evaluation of background effects and dye bias effects on different percentiles
}
\description{
For each probe type, and for each sample, several percentiles are plotted against background intensity, and also against dye bias. 
}
\usage{
quantilePlots(quantiles,backgroundInfo,designInfo,percentilesI=NULL,percentilesII=NULL)
}
\arguments{
\item{quantiles}{
A \code{list} containing three matrices. \code{list$green}, \code{list$red} and \code{list$II} must contain respectively the matrices of percentiles obtained from a \code{betaMatrix} for the Type I Green probes, Type I Red probes and Type II probes. See \code{\link{getQuantiles}}.
}
 \item{designInfo}{
\code{designInfo} matrix returned by \code{\link{getDesignInfo}}
}
 \item{backgroundInfo}{
\code{"backgroundInfo"} matrix returned by \code{\link{getBackground}}
}
 \item{percentilesI}{
List of percentiles to be plotted for Type I probes. Must be a vector of integers from 1 to 100. If set to \code{NULL} (by default), the sequence (5,10,...,95) of percentiles is plotted. 
}
 \item{percentilesII}{List of percentiles to be plotted for Type II probes. Must be a vector of integers from 1 to 100. If set to \code{NULL} (by default), the sequence (10,20,...,90) of percentiles is plotted. 
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
quantilePlots(quantiles, backgroundInfo, designInfo)
}


