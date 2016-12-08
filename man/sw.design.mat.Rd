\name{sw.design.mat}
\alias{sw.design.mat}
\title{
Creates a design matrix for a Stepped Wedge Trial
}
\description{
Constructs a basic SWT design matrix
}
\usage{
sw.design.mat(I, J, H = NULL)
}
\arguments{
  \item{I}{
Number of clusters
}
  \item{J}{
Number of time points
}
  \item{H}{
Number of units randomised at each time point
}
}
\value{
Returns a design matrix X}
\references{
Baio, G; Copas, A; Ambler, G; Hargreaves, J; Beard, E; and Omar, RZ Sample size calculation for 
a stepped wedge trial. Trials, 16:354. Aug 2015.
}
\author{
Gianluca Baio
}

