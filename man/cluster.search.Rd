\name{cluster.search}
\alias{cluster.search}
\title{
Finds the optimum number of clusters or time points
}
\description{
Given some inputs determines the optimal combination of clusters/time points to get
a set level of power.
}
\usage{
cluster.search(target.power=NULL, I=NULL, J=NULL ,H=NULL,K,design="cross-sec",
mu=0,b.trt,b.time=NULL,sigma.y=NULL,sigma.e=NULL,rho=NULL,sigma.a=NULL,rho.ind=NULL,
sigma.v=NULL,n.sims=1000,formula=NULL,family="gaussian",natural.scale=TRUE,
sig.level=0.05,n.cores=NULL,...)
}
\arguments{
  \item{target.power}{
The target power (eg 0.8)
}
  \item{I}{
A vector specifying the range in which to search for the optimal number of clusters,
eg \code{I=c(1,10)}
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
  \item{n.sims}{
Number of simulations to be used (default=1000)
}
  \item{formula}{
Formula describing the model to be used
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
  \item{sig.level}{
Significance level (default=0.05)
}
  \item{n.cores}{
Specifies the number of processors to be used for the computation (default=NULL, which
means that R will try and figure out)
}
\item{...}{Additional arguments}
}
\value{
\item{Optimum_I}{
The value of the optimal number of clusters
}
\item{power}{The estimated power in correspondence of the optimal I}
\item{time2run}{Computational time}
}
\references{
Baio, G; Copas, A; Ambler, G; Hargreaves, J; Beard, E; and Omar, RZ Sample size 
calculation for a stepped wedge trial. Trials, 16:354. Aug 2015.

Hussey M and Hughes J. Design and analysis of stepped wedge cluster randomized trials.
Contemporary Clinical Trials. 28(2):182-91. Epub 2006 Jul 7. Feb 2007
}
\author{
Rosie Leach
}
\examples{
#cluster.search(I=c(4,10),target.power=.8,J=6,K=30,mu=1.5,b.trt=.8,rho=0,
#family="poisson",n.sims=10)
}
\keyword{Stepped wedge design}
