\name{HH.binary}
\alias{HH.binary}
\title{
Power calculation for binary outcome based on analytic formula of Hussey and Hughes
}
\description{
Sample size calculations for binary outcomes based on the formula provided by 
Hussey and Hughes (2007)
}
\usage{
HH.binary(p1, OR, I, J, K, rho = 0, sig.level = 0.05, which.var = "within",X=NULL)
}
\arguments{
  \item{p1}{
Baseline probability of the outcome (for the controls)
}
  \item{OR}{
Value of the expected Odds Ratio (for the intervention vs control)
}
  \item{I}{
Number of clusters
}
  \item{J}{
Number of time points
}
  \item{K}{
Average size of each cluster
}
  \item{rho}{
Intra-class correlation coefficient (default=0)
}
  \item{sig.level}{
Significance level (default=0.05)
}
  \item{which.var}{
String character specifying which variance to be considered (options are the
default value \code{'within'} or \code{'total'}
}
  \item{X}{
A design matrix for the stepped wedge design, indicating the time at which each
of the clusters should switch the active intervention. By default is NULL
and automatically computed, but can be passed as an extra argument as a user-defined
matrix with I rows and (J+1) columns
}
}
\value{
\item{power}{
The resulting power
}
\item{sigma.y}{The estimated total (marginal) sd for the outcome}
\item{sigma.e}{The estimated residual sd}
\item{sigma.a}{The resulting cluster-level sd}
\item{setting}{A list including the following values:
 - n.clusters = The number of clusters
 - n.time.points = The number of 'active' time points
 - avg.cluster.size = The average cluster size
 - design.matrix = The design matrix for the SWT under consideration
}
}
\references{
Baio, G; Copas, A; Ambler, G; Hargreaves, J; Beard, E; and Omar, RZ Sample size 
calculation for a stepped wedge trial. Trials, 16:354. Aug 2015.

Hussey M and Hughes J. Design and analysis of stepped wedge cluster randomized trials.
Contemporary Clinical Trials. 28(2):182-91. Epub 2006 Jul 7. Feb 2007
}
\author{
Gianluca Baio
}
\examples{
HH.binary(p1=.26,OR=.55,I=10,J=5,K=20,rho=.2)}
\keyword{Hussey and Hughes formula}
