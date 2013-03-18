\name{getCoefficients}
\alias{getCoefficients}
\title{
 Return the coefficients from the ARRm linear model
}
\description{
For each probe type, it returns the coefficients of the linear model used in the ARRm normalization. Since the model is applied to each percentile separately, different coefficients are returned for every percentile. Residuals are returned as well. 
}
\usage{
getCoefficients(quantiles,designInfo,backgroundInfo,outliers.perc=0.02)
}
\arguments{
  \item{quantiles}{
A \code{list} containing three matrices. \code{"$green"}, \code{"$red"} and \code{"$II"} must contain respectively the matrices of percentiles obtained from a \code{"betaMatrix"} for the Type I Green probes, Type I Red probes and Type II probes. See \code{\link{getQuantiles}}.
}
 \item{designInfo}{
matrix returned by \code{\link{getDesignInfo}}
}
 \item{backgroundInfo}{
matrix returned by \code{\link{getBackground}}
}
 \item{outliers.perc}{
Percentage of outliers to be removed in the regression. By default, set to 0.02
}
}
\value{
	Returns a \code{list} containing three lists of coefficients for each probe type. (\code{$green} to access coefficients for Type I green probes, \code{$red} to access coefficients for Type I red probes and \code{$II} to access coefficients for Type II probes). Each list of coefficients contains five subfields. \code{res} is a matrix of residuals for the linear model across percentiles (a vector of residuals for each percentile), \code{background.vector} is a vector containing the regression coefficients for background intensity across percentiles; \code{dyebias.vector} is a vector containing the regression coefficients for dye bias across percentiles; \code{chip.variations} is a matrix of chip variations estimated by the linear model; rows correspond to percentiles, columns correspond to chips; \code{position.variations} is a matrix of position deviation from the chip mean estimated by the linear model; rows correspond to percentiles, columns correspond to positions.
}
	
\author{
Jean-Philippe Fortin <jfortin@jhsph.edu>
}
\examples{
data(greenControlMatrix)
data(redControlMatrix)
data(sampleNames)
data(betaMatrix)
backgroundInfo=getBackground(greenControlMatrix,redControlMatrix)
designInfo=getDesignInfo(sampleNames)
quantiles=getQuantiles(betaMatrix)
coefficients=getCoefficients(quantiles,designInfo,backgroundInfo)
}


