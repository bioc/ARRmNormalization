\name{getDesignInfo}
\alias{getDesignInfo}
\title{
 Build the chip and position indices }
\description{
If a vector of sample names of the form "6793856729_R03C02" is given, the function builds a data frame containing chip and position indices for the samples. If no samples names are provided by the user but explicit postion and chip vectors are provided, the data frame is built with these explicit indices. 
}
\usage{
getDesignInfo(sampleNames = NULL, chipVector = NULL, positionVector = NULL)
}
\arguments{
  \item{sampleNames}{
     Names of the samples of the form "6793856729_R03C02" (Chip ID, Row, Column)
}
  \item{chipVector}{
 Numeric vector of chip indices (one chip contains 12 samples)
}
  \item{positionVector}{
 Numeric vector of on-chip position indices (between 1 and 12) 
}
}

\value{
	
	A \code{data.frame} containing a column named \code{chipInfo} containing the chip indices, a column named \code{positionInfo} containing the position indices, and a column \code{sampleNames} if sample names were provided. 

}
\author{
Jean-Philippe Fortin <jfortin@jhsph.edu>
}

\examples{
data(sampleNames)
getDesignInfo(sampleNames)
}

