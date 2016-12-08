\name{make.swt}
\alias{make.swt}
\title{
Simulates a 'virtual' Stepped Wedge trial 
}
\description{
Simulates trial data for a SWT with normally distributed outcome
}
\usage{
make.swt(I=NULL,J=NULL,H=NULL,K,design="cross-sec",mu=NULL,b.trt,
b.time=NULL,sigma.y=NULL,sigma.e=NULL,rho,sigma.a=NULL,rho.ind=NULL,
sigma.v=NULL,X=NULL,family="gaussian",natural.scale=TRUE)
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
  \item{K}{
Average size of each cluster
}
  \item{design}{
type of design. Can be \code{'cross-sec'} (default) or \code{'cohort'} (repeated measurements)
}
  \item{mu}{
baseline outcome value
}
  \item{b.trt}{
Treatment effect
}
  \item{b.time}{
Time effect
}
  \item{sigma.y}{
total standard deviation
}
  \item{sigma.e}{
individual standard deviation
}
  \item{rho}{
Intra-class correlation coefficient
}
  \item{sigma.a}{
the sd of the the cluster-level intercept (default at NULL)
}
  \item{rho.ind}{
individual-level ICC (for cohorts)
}
  \item{sigma.v}{
the sd of the cluster-level slope (by intervention, default at NULL)
}
\item{X}{
A design matrix for the SWT. Default at NULL (will be computed automatically)
}
\item{family}{
The model family to be used. Default value is 'gaussian' and other possibile
choices are 'binomial' or 'poisson'
}
\item{natural.scale}{
Indicator for whether the input is passed on the natural scale or on the scale of
the linear predictor. By default is set to TRUE. In the case of family='gaussian'
it does not have any effect, since the link for the linear predictor is the identity.
But for family='binomial' or family='poisson', the user has to specify when the input
is given on the logit or log scale
}
}
\value{
\item{data}{A data frame containing the resulting simulated dataset}
}
\references{
Baio, G; Copas, A; Ambler, G; Hargreaves, J; Beard, E; and Omar, RZ Sample size calculation for 
a stepped wedge trial. Trials, 16:354. Aug 2015.
}
\author{
Gianluca Baio, Rosie Leach
}
\seealso{
See Also \code{\link{sim.power}}
}
\keyword{Stepped wedge design}
\keyword{Continuous outcome}
\keyword{Trial data simulation}
