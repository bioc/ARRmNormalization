\name{getBackground}
\alias{getBackground}
\title{
 Estimate background intensity from the negative control probes
}
\description{
This function estimates background intensity for the two colors by taking the median of the negative control probes in each color channel. 
}
\usage{
getBackground(greenControlMatrix, redControlMatrix)
}
\arguments{
\item{greenControlMatrix}{
matrix of negative control probes intensities in the green channel. Rows are probes, columns are samples.
}
\item{redControlMatrix}{
matrix of the negative control probes intensities in the red channel. Rows are probes, columns are samples.}
}
\value{
	Returns a \code{data.frame} with two columns; \code{"green"} contains the background intensity in the green channel for each sample and \code{"red"} contains the background intensity in the red channel for each sample 
}
\author{
Jean-Philippe Fortin <jfortin@jhsph.edu>
}
\examples{
data(greenControlMatrix)
data(redControlMatrix)
getBackground(greenControlMatrix,redControlMatrix)
}


