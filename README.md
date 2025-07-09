# SWSamp

SWSamp is a general purpose package to provide a suite of functions for the sample size calculations and power analysis in a Stepped Wedge Trial. Contains functions for closed-form sample size calculation (based on a set of specific models) and simulation-based procedures that can extend the basic framework. A detailed tutorial is available [here](https://gianluca.statistica.it/software/swsamp/tutorial.html), while general instructions and references are available [here](https://gianluca.statistica.it/software/swsamp/) and [here](https://gianluca.statistica.it/research/steppedwedge/).

## Installation
`SWSamp` can be installed running the following code in your `R` terminal.
```R
install.packages(
   "SWSamp", 
   repos = c(
      "https://giabaio.r-universe.dev", 
      "https://cloud.r-project.org",
      "https://inla.r-inla-download.org/R/stable"
   )
)
```
As every other package stored in `GitHub`, you could use `remotes`, but we recommend using the `r-universe` installation, as it is easier to maintain.
