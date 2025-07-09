# SWSamp 0.3.21

## July 2025

* Refactors the package to prevents failures in compilation for `r-universe.dev`

# SWSamp 0.3.2

* Fixes a minor error in the code for `make.swt`. The previous version would
compute the value for `sigma.a` wrongly, in the case where `(is.null(sigma.e) & is.null(sigma.a))`
The revised version correctly uses `sqrt(rho*sigma.y^2)`. Thanks to George
Papadonatos who has spotted this.
