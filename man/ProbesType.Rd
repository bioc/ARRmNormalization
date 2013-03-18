\name{ProbesType}
\alias{ProbesType}
\docType{data}
\title{
Probe Design information for the 450k methylation assay
}
\description{
Probe Design information for the Illumina Infinium HumanMethylation 450k array. To each probe is associated the design type, either Infinium I Green, Infinium I Red or Infinium II. Probe names follows Illumina's annotation (names of the form "cg00000029").
}
\usage{data(ProbesType)}
\format{
A data frame containing two columns. \code{$Probe_Name} contains the names of the probes, and \code{$Design_Type} contains the design information ("I Green", "I Red" or "II").
}
\examples{
data(ProbesType)
}

