# SWSamp [![Travis-CI Build Status](https://travis-ci.org/giabaio/SWSamp.svg?branch=master)](https://travis-ci.org/giabaio/SWSamp)[![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/giabaio/SWSamp?branch=master&svg=true)](https://ci.appveyor.com/project/giabaio/SWSamp)[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/SWSamp)](https://cran.r-project.org/package=SWSamp)[![CRAN_Download_Badge](http://cranlogs.r-pkg.org/badges/SWSamp)](https://cran.r-project.org/package=SWSamp)
SWSamp is a general purpose package to provide a suite of functions for the sample size calculations and power analysis in a Stepped Wedge Trial. Contains functions for closed-form sample size calculation (based on a set of specific models) and simulation-based procedures that can extend the basic framework.

## Installation
There are two ways of installing `SWSamp`. A "stable" version is packaged and binary files are available for Windows and as source. To install the stable version on a Windows machine, run the following commands
```R
install.packages("SWSamp",
	repos=c("http://www.statistica.it/gianluca/R",
		"https://cran.rstudio.org",
		"https://www.math.ntnu.no/inla/R/stable"),
	dependencies=TRUE
)
```
Note that you need to specify a vector of repositories - the first one hosts `SWSamp`, while the second one should be an official [CRAN mirror](https://cran.r-project.org/index.html). You can select whichever one you like, but a CRAN mirror must be provided, so that `install.packages()` can also install the "dependencies" (e.g. other packages that are required for `SWSamp` to work). The third one is used to install the package [`INLA`](http://www.r-inla.org/), which can be used to perform simulation-based sample size calculations using a Bayesian approach. This process can be quite lengthy, if you miss many of the relevant packages.

To install from source (e.g. on a Linux machine), run
```R
install.packages("SWSamp",
	repos=c("http://www.statistica.it/gianluca/R",
		"https://cran.rstudio.org",
		"https://www.math.ntnu.no/inla/R/stable"),
	type="source",
	dependencies=TRUE

```

The second way involves using the "development" version of `SWSamp` - this will usually be updated more frequently and may be continuously tested. On Windows machines, you need to install a few dependencies, including [Rtools](https://cran.r-project.org/bin/windows/Rtools/) first, e.g. by running
```R
pkgs <- c("foreach", "doParallel", "iterators", "parallel", "Matrix","lme4","INLA","Rtools","devtools")
repos <- c("https://cran.rstudio.com", "https://www.math.ntnu.no/inla/R/stable") 
install.packages(pkgs,repos=repos,dependencies = "Depends")
```
before installing the package using `devtools`:
```R
devtools::install_github("giabaio/SWSamp")
```
Under Linux or MacOS, it is sufficient to install the package via `devtools`:
```R
install.packages("devtools")
devtools:install_github("giabaio/SWSamp")
```
