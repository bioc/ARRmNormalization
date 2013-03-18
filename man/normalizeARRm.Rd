\name{normalizeARRm}
\alias{normalizeARRm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Perform ARRm normalization 
}
\description{
This function perform Adaptive Robust Regression method (ARRm) normalization on Beta values. The method corrects for background intensity, dye bias and spatial on-chip position. By default, chip mean correction is also performed. 
}
\usage{
normalizeARRm(betaMatrix, designInfo, backgroundInfo, outliers.perc = 0.02, 
goodProbes = NULL,chipCorrection=TRUE)
}

\arguments{
  \item{betaMatrix}{
\code{matrix} containing the Beta values. Rows are probes, columns are samples.
}
  \item{designInfo}{
A \code{data.frame} containing a column named \code{chipInfo} containing the chip indices and a column named \code{positionInfo} containing the position indices
}
  \item{backgroundInfo}{
A \code{data.frame} containing two columns: \code{green} contains the background intensity in the green channel for each sample and \code{red} contains the background intensity in the red channel for each sample
}
  \item{outliers.perc}{
Proportion (between 0 and 1) of outliers to be removed from the ARRm regression 
}

  \item{goodProbes}{
Ids of the probes to be normalized (Id. of the form "cg00000029")
}
 \item{chipCorrection}{
logical, should normalization correct for chip mean?
}
}
\value{
A \code{matrix} containing the normalized Beta values
}

\author{
Jean-Philippe Fortin <jfortin@jhsph.edu>
}


\seealso{
\code{\link{getBackground}} to see how to obtain background information from control probes, and \code{\link{getDesignInfo}} to see how to obtain position and chip indices
}
\examples{
data(greenControlMatrix)
data(redControlMatrix)
data(sampleNames)
data(betaMatrix)
backgroundInfo=getBackground(greenControlMatrix, redControlMatrix)
designInfo=getDesignInfo(sampleNames)
normMatrix=normalizeARRm(betaMatrix, designInfo, backgroundInfo, outliers.perc = 0.02)
}

