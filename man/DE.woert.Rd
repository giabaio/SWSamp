\name{DE.woert}
\alias{DE.woert}
\title{Computes the Design Effect for a Stepped Wedge Trial} 
\description{
Sample size calculations for a SWT using a cross-sectional design. 
This is based on (the correct version) of Woertman et al (2013), as described in 
Baio et al (2015). 
}
\usage{
DE.woert(outcome = "cont", input, K, J, B = 1, T = 1, rho, sig.level = 0.05, power = 0.8)
}
\arguments{
  \item{outcome}{
String. Type of outcome. Options are \code{cont}, \code{bin} or \code{count}
}
  \item{input}{
input = a list containing the arguments. This differs depending on the type of 
outcome, as follows:
 - continuous outcome: 
   1) delta (treatment effect)
   2) sd (standard deviation)
 - binary outcome:
   1) p1 (baseline probability of outcome)
   2) either p2 (treatment probability of outcome), or OR (treatment effect as OR)
 - count outcome: 
   1) r1 (baseline rate of outcome)
   2) either r2 (treatment rate of outcome), or RR (treatment effect as RR)
}
  \item{K}{
average cluster size
}
  \item{J}{
number of time points (excluding baseline)
}
  \item{B}{
number of baseline measurement times
}
  \item{T}{
number of measurement times during each crossover
}
  \item{rho}{
ICC
}
  \item{sig.level}{
significance level (default = 0.05)
}
  \item{power}{
Power (default = 0.8)
}
}
\value{
\item{n.cls.swt}{
Number of clusters required to reach the pre-specified power with the given significance 
level.
}
\item{n.pts}{
The total number of participants required.
}
\item{DE.woert}{
The resulting Design Effect.
}
\item{CF}{
The resulting Correction Factor.
}
\item{n.rct}{
The original individual RCT sample required to reach the pre-specified power with the 
given significance level.
}
}
\references{
Baio, G; Copas, A; Ambler, G; Hargreaves, J; Beard, E; and Omar, RZ Sample size calculation for 
a stepped wedge trial. Trials, 16:354. Aug 2015.
}
\author{
Gianluca Baio
}
\examples{
# Continuous outcome
input <- list(delta=-0.3875,sd=1.55)
K <- 20
J <- 5
rho <- .2
DE.woert(input=input,K=K,J=J,rho=rho)
#
# Binary outcome
input <- list(OR=.53,p1=.26)
DE.woert(outcome="bin",input=input,K=K,J=J,rho=rho)
#
# Count outcome
input <- list(RR=.8,r1=1.5)
DE.woert(outcome="count",input=input,K=K,J=J,rho=rho)
}
\keyword{Design Effect}
