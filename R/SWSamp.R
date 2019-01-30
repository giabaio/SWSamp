#' Simulates a 'virtual' Stepped Wedge trial
#' 
#' Simulates trial data for a SWT with normally distributed outcome
#' 
#' 
#' @param I Number of clusters
#' @param J Number of time points
#' @param H Number of units randomised at each time point
#' @param K Average size of each cluster
#' @param design type of design. Can be \code{'cross-sec'} (default) or
#' \code{'cohort'} (repeated measurements)
#' @param mu baseline outcome value
#' @param b.trt Treatment effect
#' @param b.time Time effect
#' @param sigma.y total standard deviation
#' @param sigma.e individual standard deviation
#' @param rho Intra-class correlation coefficient
#' @param sigma.a the sd of the the cluster-level intercept (default at NULL)
#' @param rho.ind individual-level ICC (for cohorts)
#' @param sigma.v the sd of the cluster-level slope (by intervention, default
#' at NULL)
#' @param X A design matrix for the SWT. Default at NULL (will be computed
#' automatically)
#' @param family The model family to be used. Default value is 'gaussian' and
#' other possibile choices are 'binomial' or 'poisson'
#' @param natural.scale Indicator for whether the input is passed on the
#' natural scale or on the scale of the linear predictor. By default is set to
#' TRUE. In the case of family='gaussian' it does not have any effect, since
#' the link for the linear predictor is the identity. But for family='binomial'
#' or family='poisson', the user has to specify when the input is given on the
#' logit or log scale
#' @return \item{data}{A data frame containing the resulting simulated dataset}
#' @author Gianluca Baio, Rosie Leach
#' @seealso See Also \code{\link{sim.power}}
#' @references Baio, G; Copas, A; Ambler, G; Hargreaves, J; Beard, E; and Omar,
#' RZ Sample size calculation for a stepped wedge trial. Trials, 16:354. Aug
#' 2015.
#' @keywords Stepped wedge design Continuous outcome Trial data simulation
#' @export make.swt
make.swt <- function(I=NULL,J=NULL,H=NULL,K,design="cross-sec",mu=NULL,b.trt,
                     b.time=NULL,sigma.y=NULL,sigma.e=NULL,rho,sigma.a=NULL,rho.ind=NULL,
                     sigma.v=NULL,X=NULL,family="gaussian",natural.scale=TRUE)  {
  ## Normal outcome
  ## inputs:
  # I = number of clusters
  # J = number of randomisation time points (excluding baseline)
  # K = average sample size in each cluster
  # H = number of units randomised at each time point
  # design = type of design. Can be "cross-sec" (default) or "cohort" (repeated measurements)
  # mu = baseline mean
  # b.trt = intervention effect
  # b.time = (optional) time effect (on linear time trend scale!)
  # sigma.y = *total* individual variability
  # sigma.e = *residual* individual variability
  # rho = ICC at the cluster level
  # sigma.a = the sd of the cluster random effect (default at NULL)
  # rho.ind = ICC at the individual level (for cohort design, default at NULL)
  # sigma.v = the sd of the treatment random effect (default at NULL)
  # X = SW design matrix (default at NULL and will be computed automatically)
  # family = type of outcome to be simulated (options are "gaussian", "binomial" and "poisson")
  # natural.scale = whether the input values for the trt effect etc are passed on the natural scale or not
  #
  # GB + Rosalind Leach (November 2015)
  
  # CHECK ON VALID FAMILY TYPE 
  flag=1
  if(family=="gaussian" | family=="binomial" | family=="poisson"){flag=0}
  if(flag==1){ stop("Available options for the argument 'family' are: gaussian, binomial or poisson")}

  # CREATE DESIGN MATRIX FROM sw.design.mat IF X HAS NOT BEEN PROVIDED
  if(is.null(X)) {X <- sw.design.mat(I=I,J=J,H=NULL)} else {row.names(X) <- sample(1:I,I)}
  
  if(family == "gaussian"){
    # SETS DEFAULT FOR mu AND b.time IF THEY HAVE NOT BEEN PORVIDED
    if(is.null(mu)) {mu=0}                          # Baseline (set to 0 if not specified)
    if(is.null(b.time)) {                           # Time effect (if not specified, set to)
      b.time <- .5*b.trt   }                         #   some pre-defined proportion of the
    # SIMULATES VARIABLES FROM THE INPUTS GIVEN
    treatment <- rep(t(X),each=K,1)                 # Treatment allocation at each time point
    time <- rep(seq(0,J),each=K,I)                  # K measurements during the trial
    person <- seq(1,K)                              # individuals IDs
    cluster.order <- rep(sample(1:I))               # order at which clusters switch
    cluster <- rep(cluster.order, each=(J+1)*K)     # cluster IDs (randomised)
    id <- NULL           
    # Calculate the hyperparameters
    mu.a <- 0					                              # Within cluster (random intercept) mean 
    
    # CROSS-SECTIONAL DESIGN -------------------------------------------------
    if(design=="cross-sec") {
      # CHECK THE RIGHT COMBINATION OF sigma.e, sigma.y, sigma.a, rho HAS BEEN GIVEN
      if((sum(sapply(list(sigma.y,sigma.e),is.null)) !=1) & (sum(sapply(list(sigma.a,rho),is.null)) !=1)) { 	# Need to pass 1 of 'sigma.y' or 'sigma.e' and 1 of 'sigma.a' or 'rho
        stop("Please provide either 'sigma.y' (*total* individual variability) or 'sigma.e' (*residual* variability)
             and either 'sigma.a' (*within cluster* variability) or 'rho' (ICC)")}
      
      # IF STATEMENTS TO CALCULATE ALL RELEVANT PARAMETERS FROM THE INPUTS GIVEN
      if(is.null(sigma.e) & is.null(sigma.a)) {
        sigma.a <- rho*sigma.y
        sigma.e <- sqrt(sigma.y^2-sigma.a^2)
      }
      if(is.null(sigma.e) & is.null(rho)) {
        sigma.e <- sqrt(sigma.y^2-sigma.a^2)
        rho <- sigma.a^2/sigma.y^2
      }
      if(is.null(sigma.y) & is.null(sigma.a)) {
        sigma.a <- sqrt(sigma.e^2*rho/(1-rho))
        sigma.y <- sqrt(sigma.a^2+sigma.e^2)
      }
      if(is.null(sigma.y) & is.null(rho)) {
        sigma.y <- sqrt(sigma.e^2+sigma.a^2)
        rho <- sigma.a^2/sigma.y^2
      }
      
      # SIMULATE REMAINING VARIABLES
      a <- rnorm(I, mu.a, sigma.a)  		            # cluster-level intercept
      
      # CALCULATE LINEAR PREDICTOR
      linpred <- mu + a[cluster] + b.trt*treatment + b.time*time		# Calculates the linear predictor
      } 
    
    # COHORT DESIGN -------------------------------------------------
    if (design=="cohort") { 
      # CHECK THE RIGHT COMBINATION OF sigma.e, sigma.y, sigma.a, rho HAS BEEN GIVEN
      if((sum(sapply(list(sigma.y,sigma.e),is.null)) !=1) || (sum(sapply(list(sigma.a,rho),is.null)) !=1) ||
         (sum(sapply(list(sigma.v,rho.ind),is.null)) !=1)) {
        stop("Please provide either 'sigma.y' (*total* individual variability) or 'sigma.e' (*residual* variability)
             AND either 'sigma.a' (*within cluster* variability) or 'rho' (cluster-level ICC)
             AND either 'sigma.v' (*within cluster-individual* variability) or 'rho.ind (individual-level ICC)") }
      
      # IF STATEMENTS TO CALCULATE ALL RELEVANT PARAMETERS FROM THE INPUTS GIVEN
      if(is.null(sigma.y) & is.null(sigma.a) & is.null(sigma.v)) {
        r1 <- rho.ind/(1-rho.ind)
        r2 <- rho/(1-rho)
        sigma.v <- sqrt((r1*(sigma.e^2/(1-rho)))/(1-r1*r2))
        sigma.a <- sqrt(r2*(sigma.e^2+sigma.v^2))
        sigma.y <- sqrt(sigma.e^2+sigma.a^2+sigma.v^2)
      }
      if(is.null(sigma.y) & is.null(rho) & is.null(sigma.v)) {
        r <- rho.ind/(1-rho.ind)
        sigma.v <- sqrt(r*(sigma.a^2+sigma.e^2))
        rho <- sigma.a^2/(sigma.a^2+sigma.e^2+sigma.v^2)
        sigma.y <- sqrt(sigma.a^2+sigma.e^2+sigma.v^2)
      }
      if(is.null(sigma.y) & is.null(sigma.a) & is.null(rho.ind)) {
        r <- rho/(1-rho)
        sigma.a <- sqrt(r*(sigma.v^2+sigma.e^2))
        sigma.y <- sqrt(sigma.a^2+sigma.e^2+sigma.v^2)
		rho.ind <- sigma.v^2/(sigma.a^2+sigma.e^2+sigma.v^2)
      }
      if(is.null(sigma.y) & is.null(rho) & is.null(rho.ind)) {
        rho <- sigma.a^2/(sigma.a^2+sigma.e^2+sigma.v^2)
        rho.ind <- sigma.v^2/(sigma.a^2+sigma.e^2+sigma.v^2)
        sigma.y <- sqrt(sigma.a^2+sigma.e^2+sigma.v^2)
      }
      if(is.null(sigma.e) & is.null(sigma.a) & is.null(sigma.v)) {
        sigma.a <- sqrt(sigma.y^2*rho)
        sigma.v <- sqrt(sigma.y^2*rho.ind)
        sigma.e <- sqrt(sigma.y^2-sigma.a^2-sigma.v^2)
      }
      if(is.null(sigma.e) & is.null(rho) & is.null(sigma.v)) {
        rho <- sigma.a^2/sigma.y^2
        sigma.v <- rho.ind/sigma.y^2
        sigma.e <- sqrt(sigma.y^2-sigma.a^2-sigma.v^2)
      }
      if(is.null(sigma.e) & is.null(sigma.a) & is.null(rho.ind)) {
        sigma.a <- sqrt(rho*sigma.y^2)
        rho.ind <- sigma.v^2/sigma.y^2
        sigma.e <- sqrt(sigma.y^2-sigma.a^2-sigma.v^2)
      }
      if(is.null(sigma.e) & is.null(rho) & is.null(rho.ind)) {
        rho <- sigma.a^2/sigma.y^2
        rho.ind <- sigma.v^2/sigma.y^2
        sigma.e <- sqrt(sigma.y^2-sigma.a^2-sigma.v^2)
      }
      
      # SIMULATE REMAINING VARIABLES
      start <- seq(1,I*K,by=K)
      id <- numeric                                 # Individual ID (for repeated measurements)
      id <- rep(seq(start[cluster.order[1]],start[cluster.order[1]]+K-1),(J+1))
      for (i in 2:I) {
        id <- c(id,rep(seq(start[cluster.order[i]],start[cluster.order[i]]+K-1),(J+1)))
      }
      mu.v <- 0                                     # within individuals (random intercept) mean
      v <- rnorm(K*I,mu.v,sigma.v)                  # individual-level intercept
      sigma.a <- sqrt(sigma.y^2*rho)                # Within cluster (random intercept) sd 
      a <- rnorm (I, mu.a, sigma.a)                 # cluster-level intercept 
      
      # CALCULATE LINEAR PREDICTOR
      linpred <- mu + a[cluster] + v[id] + b.trt*treatment + b.time*time
      } 
    # SIMULATES DATA 
    n <- I*(J+1)*K
    y <- rnorm(I*(J+1)*K, linpred, sigma.e)
  }
  
  if(family=="binomial"){
    # SETS DEFAULT FOR mu AND b.time IF THEY HAVE NOT BEEN PORVIDED
    if(is.null(mu)) {
      if (natural.scale==TRUE) {
        mu=.5		}		      # Default baseline probability on [0-1]
      else if (natural.scale==FALSE) {
        mu=log(.5/.5)	}		# Default baseline probability on logit scale
    } 
    
    # SIMULATES VARIABLES FROM THE INPUTS GIVEN
    treatment <- rep(t(X),each=K,1)                 # Treatment allocation at each time point
    time <- rep(seq(0,J),each=K,I)                  # K measurements during the trial
    person <- seq(1,K)                              # individuals IDs
    cluster.order <- rep(sample(1:I))               # order at which clusters switch
    cluster <- rep(cluster.order, each=(J+1)*K)   	# cluster IDs (randomised)
    id <- NULL
    
    # If natural.scale = T, then assume that the user is giving values for p0 (mu) and OR directly (b.trt)
    if (natural.scale==TRUE) {
      OR <- b.trt			              # for simplicity defines the OR
      b.trt <- log(b.trt)			      # logOR to be used in the linear predictor
      p0 <- mu				              # for simplicity defines the baseline probability
      mu <- log(mu/(1-mu))          # Converts baseline probability to log odds to be used in the linear predictor
    }  
    # But if natural.scale = F, then the user has passed data on the logit scale for p0 (mu) and on the logOR (b.trt)
    if (natural.scale==FALSE) {
      OR <- exp(b.trt)			        # Defines the OR on the natural scale
      p0 <- exp(mu)/(1+exp(mu))		  # Rescales the baseline to the [0-1] range
    }
    # Now defines p1 and sigma.e
    p1 <- OR*(p0/(1-p0))/(1+OR*(p0/(1-p0)))	      # Estimates the outcome probability for intervention
    sigma.e <- sqrt(((p0*(1-p0))+(p1*(1-p1)))/2)	# Pooled estimate of the within cluster sd
    if(is.null(b.time)) {                         # Time effect (if not specified, set to)
      b.time <- .5*b.trt  }                       #   some pre-defined proportion of the treatment effect
    
    # CALCULATE THE HYPERPARAMETERS
    mu.a <- 0		                              # Within cluster (random intercept) mean 
    
    # CROSS-SECTIONAL DESIGN -------------------------------------------------
    if(design=="cross-sec") {
      # CHECK THE RIGHT COMBINATION OF sigma.a, rho HAS BEEN GIVEN
      if (sum(sapply(list(rho,sigma.a),is.null)) !=1) {			 #cannot provide both sigma.a and rho
        stop("exactly one of 'rho' and 'sigma.a' must be null")  }	
      
      # IF STATEMENTS TO CALCULATE ALL RELEVANT PARAMETERS FROM THE INPUTS GIVEN
      if(is.null(sigma.a)){							# if sigma.a is not given, calculate it as a function of rho
        sigma.a <- sqrt(sigma.e^2*rho/(1-rho))  }       # Within cluster (random intercept) sd
      if(is.null(rho)){							# if sigma.a is not given, calculate it as a function of rho
        rho <-  (sigma.a^2/(sigma.e^2 + sigma.a^2)) }    # if rho is not given, calculate it as a function of sigma.a
      
      # SIMULATE REMAINING VARIABLES
      a <- rnorm(I, mu.a, sigma.a)  		            # cluster-level intercept 
      
      # CALCULATE LINEAR PREDICTOR ON THE LOGIT SCALE!
      linpred <- mu + a[cluster] + b.trt*treatment + b.time*time
    } 

    # COHORT DESIGN -------------------------------------------------
    if (design=="cohort") {  
      # CHECK THE RIGHT COMBINATION OF sigma.a, sigma.v, rho, rho.ind HAS BEEN GIVEN
	  if((sum(sapply(list(rho, sigma.a),is.null)) !=1) || (sum(sapply(list(sigma.a,rho),is.null)) !=1) ||
         (sum(sapply(list(sigma.v,rho.ind),is.null)) !=1)) {
        stop("Please provide either 'sigma.a' (*within cluster* variability) or 'rho' (cluster-level ICC)
             AND either 'sigma.v' (*within cluster-individual* variability) or 'rho.ind (individual-level ICC)") }
      
      # IF STATEMENTS TO CALCULATE ALL RELEVANT PARAMETERS FROM THE INPUTS GIVEN
    if(is.null(sigma.v) & is.null(sigma.a)) {
        r1 <- rho.ind/(1-rho.ind)
        r2 <- rho/(1-rho)
        sigma.v <- sqrt((r1*(sigma.e^2/(1-rho)))/(1-r1*r2))
        sigma.a <- sqrt(r2*(sigma.e^2+sigma.v^2))
      }

	if(is.null(rho.ind) & is.null(sigma.a)) {
         r <- rho/(1-rho)
        sigma.a <- sqrt(r*(sigma.v^2+sigma.e^2))
		#rho.ind <-
      }

	if(is.null(rho.ind) & is.null(rho)) {
         rho <- sigma.a^2/(sigma.a^2+sigma.e^2+sigma.v^2)
        rho.ind <- sigma.v^2/(sigma.a^2+sigma.e^2+sigma.v^2)
      }

	if(is.null(rho) & is.null(sigma.v)) {
         r <- rho.ind/(1-rho.ind)
        sigma.v <- sqrt(r*(sigma.a^2+sigma.e^2))
        rho <- sigma.a^2/(sigma.a^2+sigma.e^2+sigma.v^2)
      }
      
      # SIMULATE REMAINING VARIABLES
      start <- seq(1,I*K,by=K)
      id <- numeric                                 # Individual ID (for repeated measurements)
      id <- rep(seq(start[cluster.order[1]],start[cluster.order[1]]+K-1),(J+1))
      for (i in 2:I) {
        id <- c(id,rep(seq(start[cluster.order[i]],start[cluster.order[i]]+K-1),(J+1)))
      }
      mu.v <- 0                                     # within individuals (random intercept) mean
      v <- rnorm(K*I,mu.v,sigma.v)                  # individual-level intercept
###      sigma.a <- sqrt(sigma.y^2*rho)                # Within cluster (random intercept) sd 
      a <- rnorm (I, mu.a, sigma.a)                 # cluster-level intercept 
      
      # CALCULATE LINEAR PREDICTOR ON THE LOGIT SCALE!
      linpred <- mu + a[cluster] + v[id] + b.trt*treatment + b.time*time
    } 

    # SIMULATES DATA 
    n <- I*(J+1)*K
    # Converts log odds back to probability
    p = (exp(linpred))/(1+exp(linpred))
    y <- rbinom(n, 1, p)
  }
  
  if(family == "poisson"){
    # SETS DEFAULT FOR mu AND b.time IF THEY HAVE NOT BEEN PROVIDED
    if(is.null(mu)) {
      if (natural.scale==TRUE) {
        mu=1}		# Default baseline rate set to 1
      else if (natural.scale==FALSE) {
        mu=log(1)}		# Default baseline rate on log scale
    } 
    
    # SIMULATES VARIABLES FROM THE INPUTS GIVEN
    treatment <- rep(t(X),each=K,1)                 # Treatment allocation at each time point
    time <- rep(seq(0,J),each=K,I)                  # K measurements during the trial
    person <- seq(1,K)                              # individuals IDs
    cluster.order <- rep(sample(1:I))               # order at which clusters switch
    cluster <- rep(cluster.order, each=(J+1)*K)     # cluster IDs (randomised)
    id <- NULL                                      # initialise ID vector
    
    # If natural.scale = T, then assume that the user is giving values for p0 and OR directly
    if (natural.scale==TRUE) {
      RR <- b.trt                 # For convenience defines the RR
      b.trt <- log(b.trt)			    # logRR to be used in the linear predictor
      lambda0 <- mu               # For convenience defines the baseline rate
      mu <- log(mu)               # Converts baseline probability to log rate to be used in the linear predictor
    }  
    # But if natural.scale = F, then the user has passed data on the logit scale for p0 and on the logOR
    if (natural.scale==FALSE) {
      RR <- exp(b.trt)			      # Defines the RR on the natural scale
      lambda0 <- exp(mu)      		# Rescales the baseline rate to the [0-\infty] range
    }
    lambda1 <- lambda0*RR                   # rate for intervention
    sigma.e <- sqrt((lambda0+lambda1)/2)		# Within cluster sd (estimated using pooled value)
    if(is.null(b.time)) {                   # Time effect (if not specified, set to)
      b.time <- .5*b.trt  }                 #  some pre-defined proportion of the

    # CALCULATE THE HYPERPARAMETERS
    mu.a <- 0		                              # Within cluster (random intercept) mean 
    
    # CROSS-SECTIONAL DESIGN -------------------------------------------------
    if(design=="cross-sec") {
      # CHECK THE RIGHT COMBINATION OF sigma.a, rho HAS BEEN GIVEN
	    if (sum(sapply(list(rho,sigma.a),is.null)) !=1) {			 #cannot provide both sigma.a and rho
        stop("exactly one of 'rho' and 'sigma.a' must be null")  }	
		      
      # IF STATEMENTS TO CALCULATE ALL RELEVANT PARAMETERS FROM THE INPUTS GIVEN
	  if(is.null(sigma.a)){							# if sigma.a is not given, calculate it as a function of rho
        sigma.a <- sqrt(sigma.e^2*rho/(1-rho))  }       # Within cluster (random intercept) sd
      if(is.null(rho)){							# if sigma.a is not given, calculate it as a function of rho
        rho <-  (sigma.a^2/(sigma.e^2 + sigma.a^2)) }    # if rho is not given, calculate it as a function of sigma.a
      	  
      # SIMULATE REMAINING VARIABLES
      a <- rnorm(I, mu.a, sigma.a)  		            # cluster-level intercept    
      
      # CALCULATE LINEAR PREDICTOR ON THE LOG SCALE
      linpred <- mu + a[cluster] + b.trt*treatment + b.time*time
    } 

    # COHORT DESIGN -------------------------------------------------
    if (design=="cohort") {  
      # CHECK THE RIGHT COMBINATION OF sigma.a, rho, rho.ind HAS BEEN GIVEN
     if((sum(sapply(list(rho, sigma.a),is.null)) !=1) || (sum(sapply(list(sigma.a,rho),is.null)) !=1) ||
         (sum(sapply(list(sigma.v,rho.ind),is.null)) !=1)) {
        stop("Please provide either 'sigma.a' (*within cluster* variability) or 'rho' (cluster-level ICC)
             AND either 'sigma.v' (*within cluster-individual* variability) or 'rho.ind (individual-level ICC)") }
      
      # IF STATEMENTS TO CALCULATE ALL RELEVANT PARAMETERS FROM THE INPUTS GIVEN
    if(is.null(sigma.v) & is.null(sigma.a)) {
        r1 <- rho.ind/(1-rho.ind)
        r2 <- rho/(1-rho)
        sigma.v <- sqrt((r1*(sigma.e^2/(1-rho)))/(1-r1*r2))
        sigma.a <- sqrt(r2*(sigma.e^2+sigma.v^2))
      }

	if(is.null(rho.ind) & is.null(sigma.a)) {
         r <- rho/(1-rho)
        sigma.a <- sqrt(r*(sigma.v^2+sigma.e^2))
		#rho.ind <-
      }

	if(is.null(rho.ind) & is.null(rho)) {
         rho <- sigma.a^2/(sigma.a^2+sigma.e^2+sigma.v^2)
        rho.ind <- sigma.v^2/(sigma.a^2+sigma.e^2+sigma.v^2)
      }

	if(is.null(rho) & is.null(sigma.v)) {
         r <- rho.ind/(1-rho.ind)
        sigma.v <- sqrt(r*(sigma.a^2+sigma.e^2))
        rho <- sigma.a^2/(sigma.a^2+sigma.e^2+sigma.v^2)
      }	  
	  
      # SIMULATE REMAINING VARIABLES
      start <- seq(1,I*K,by=K)
      id <- numeric                                 # Individual ID (for repeated measurements)
      id <- rep(seq(start[cluster.order[1]],start[cluster.order[1]]+K-1),(J+1))
      for (i in 2:I) {
        id <- c(id,rep(seq(start[cluster.order[i]],start[cluster.order[i]]+K-1),(J+1)))
      }
      mu.v <- 0                                     # within individuals (random intercept) mean
      v <- rnorm(K*I,mu.v,sigma.v)                  # individual-level intercept
###      sigma.a <- sqrt(sigma.y^2*rho)                # Within cluster (random intercept) sd 
      a <- rnorm (I, mu.a, sigma.a)                 # cluster-level intercept 
      
      # CALCULATE LINEAR PREDICTOR
      linpred <- mu + a[cluster] + v[id] + b.trt*treatment + b.time*time
    } 

    # SIMULATES DATA 
    n <- I*(J+1)*K
    y <- rpois(I*(J+1)*K, lambda=exp(linpred))
  }
  # RETURNS DATA FRAME
  data <- data.frame(y, person, time, cluster, treatment, linpred, b.trt)
  if (design=="cohort") {data <- data.frame(data,id)}
  return(data)
}





#' Power calculations based on a simulation approach
#' 
#' Simulation-based power calculations for a SWT with normally distributed
#' outcome
#' 
#' 
#' @param I Number of clusters
#' @param J Number of time points
#' @param H Number of units randomised at each time point
#' @param K Average size of each cluster
#' @param design type of design. Can be \code{'cross-sec'} (default) or
#' \code{'cohort'} (repeated measurements)
#' @param mu baseline outcome value
#' @param b.trt Treatment effect
#' @param b.time Time effect
#' @param sigma.y total standard deviation
#' @param sigma.e individual standard deviation
#' @param rho Intra-class correlation coefficient
#' @param sigma.a the sd of the the cluster-level intercept (default at NULL)
#' @param rho.ind individual-level ICC (for cohorts)
#' @param sigma.v the sd of the cluster-level slope (by intervention, default
#' at NULL)
#' @param n.sims Number of simulations to be used (default=1000)
#' @param formula Formula describing the model to be used
#' @param family The model family to be used. Default value is 'gaussian' and
#' other possibile choices are 'binomial' or 'poisson'
#' @param natural.scale Indicator for whether the input is passed on the
#' natural scale or on the scale of the linear predictor. By default is set to
#' TRUE. In the case of family='gaussian' it does not have any effect, since
#' the link for the linear predictor is the identity. But for family='binomial'
#' or family='poisson', the user has to specify when the input is given on the
#' logit or log scale
#' @param sig.level Significance level (default=0.05)
#' @param n.cores Specifies the number of processors to be used for the
#' computation (default=NULL, which means that R will try and figure out)
#' @param method A string specifying the method to be used for the calculation.
#' The default value is \code{lme}, indicating a standard frequentist analysis,
#' based on (generalised) linear model and including structured (random)
#' effects, when necessary. An alternative version is framed in a Bayesian
#' setting and uses Integrated Nested Laplace Approximation (INLA) to analyse
#' the data and obtain the relevant estimates for the posterior distributions
#' of the model parameters. This can be performed by setting
#' \code{method='inla'}
#' @param plot Shows a plot of the moving average of the resulting power after
#' 10%, 20%, ..., 100% of the simulations. This helps assessing the convergence
#' of the estimate towards some common value. The default is \code{FALSE}, in
#' which case a graph is not shown.
#' @param ...  Additional optional arguments. The user can specify a function
#' named \code{data}, which defines the simulation of the 'virtual trial data'.
#' These can be in any given form, with the only constraint that it should
#' return a data frame containing the relevant variables. The function
#' \code{data} can have no arguments at all, but this can be relaxed and there
#' can be suitable inputs to this function. In this case, the user also needs
#' to specify a list \code{inpts} including all the values for the arguments to
#' be used with the user-defined function \code{data}. When using user-defined
#' data generating processes, the user must be also pass a relevant formula,
#' depending on what the model used (for both generation of the dataset and
#' analysis) is. Another additional argument that can be passed to the call to
#' \code{sim.sw.cont} is \code{treatment}, a string specifying the name of the
#' treatment variable (if not present, \code{SWSamp} assumes that this is
#' exactly 'treatment').
#' @return \item{power}{ The resulting estimated power, for the given
#' configuration. If the model does not include random effects, this is based
#' on the p-value computed by \code{lm}, which is used to analyse the simulated
#' datasets. If the model does include random effects (which is the case for a
#' SWT), then \code{SWSample} assesses whether the 'true' effect is correctly
#' detected by computing the (1-\code{alpha})% confidence interval and checking
#' whether it is entirely above or below 0. This is because it is difficult to
#' assess the correct degrees of freedom of the resulting (linear) mixed model.
#' The p-value could be computed using the Satterthwaite approximation, or by
#' using a rougher Normal approximation, but in line with suggestions by
#' Pinheiro, J. C., and D. M. Bates.  2000. Mixed-effects models in S and
#' S-PLUS. Springer, New York, we sidestep the problem by focussing on
#' estimation, rather than hypothesis testing for this. }
#' \item{time2run}{Running time, in seconds} \item{ci.power}{Estimated 95\%
#' confidence interval for the power - based on normal approximation}
#' \item{theta}{Estimated treatment effect with standard error}
#' \item{rnd.eff.sd}{Estimated variance components} \item{setting}{A list
#' summarising the assumptions in terms of number of clusters, time points,
#' type of model, formula used}
#' @author Gianluca Baio
#' @references Baio, G; Copas, A; Ambler, G; Hargreaves, J; Beard, E; and Omar,
#' RZ Sample size calculation for a stepped wedge trial. Trials, 16:354. Aug
#' 2015.
#' @keywords Stepped wedge design Continuous outcome Simulation-based sample
#' size calculations
#' @examples
#' 
#' mu1=0.3
#' b.trt=-0.3875
#' sigma.e=1.55
#' J=5
#' K=20
#' sig.level=0.05
#' n.sims=10
#' rho=0.5
#' pow.cont <- sim.power(I=14,J=J,H=NULL,K=K,rho=rho,mu=mu1,sigma.e=sigma.e,b.trt=b.trt,
#'                        formula=NULL,n.sims=n.sims,sig.level=sig.level,n.cores=2)
#' pow.cont$power
#' pow.cont$ci
#' pow.cont$time2run.sec
#' pow.cont$rnd.eff.sd^2
#' 
#' @export sim.power
sim.power <- function (I,J,H=NULL,K,design="cross-sec",mu=0,b.trt,b.time=NULL,
                      sigma.y=NULL,sigma.e=NULL,rho=NULL,sigma.a=NULL,
                      rho.ind=NULL,sigma.v=NULL,n.sims=1000,formula=NULL,
                      family="gaussian",natural.scale=TRUE,sig.level=0.05,n.cores=NULL,
                      method="lme",plot=FALSE,...) {
  
  ## Power analysis for repeated cohort design (cross-sectional) - continuous outcome
  ## The optional arguments include the specification of a user-defined function data
  ## which is used to simulate the data, according to any design the user has.
  ## So including data=data, where data=function(){...} specifies the simulation 
  ## of a suitable dataset allows the user to select *any* kind of model for data
  ## generation. If data has inputs, then these must be specified as a list 
  ## inpts, to be passed as an extra argument to sim.sw.cont. 
  ## In addition to this, the user should then select a suitable formula
  ## which is in line with the model used for the simulation of the virtual trial.
  ## Also, the name of the "treatment" variable can be specified as a string in the
  ## extra argument treat.name. If not present then it is assumed that the name of the
  ## treatment variable is, unsurprisingly, "treatment"
  ## method is default at 'lme' indicating that a linear mixed model using lmer will be
  ## used for the analysis. The user can also specify 'INLA' or 'inla' which performs
  ## the analysis in Bayesian framework with INLA
  
  exArgs <- list(...)
  
  requireNamespace("lme4",quietly=TRUE)
  requireNamespace("foreach",quietly=TRUE)
  requireNamespace("doParallel",quietly=TRUE)
  requireNamespace("parallel",quietly=TRUE)
  requireNamespace("iterators",quietly=TRUE)
  
  # If not specified, uses all cores available to the machine
  if(is.null(n.cores)) {n.cores <- parallel::detectCores()}
  # Usually a good idea to leave 1 core free for the OS to use
  doParallel::registerDoParallel(cores=n.cores-1)
  
  # Defines basic formulation for the model
  if(is.null(formula)){
    formula=y~treatment+factor(time)+(1|cluster)
    if(design=="cohort") {formula <- update.formula(formula,.~.+(1|id))}
    
    # Now updates the formula so it's given in INLA notation if method==INLA
    if (method=="INLA" || method=="inla") {
        # Checks the the package INLA is installed
        if (!isTRUE(requireNamespace("INLA", quietly = TRUE))) {
            stop("You need to install the packages 'INLA'. Please run in your R terminal:\n install.packages('INLA', repos='http://www.math.ntnu.no/inla/R/stable')")
        } else {
            if (!is.element("INLA", (.packages()))) {
                attachNamespace("INLA")
            }
        }
        # This substitutes the notation '(1 | var)' to the notation 'f(var)'. Notice that INLA will assume that 
        # this actually means 'f(var,model="iid")' --- if specific prior are required, they need to be passed
        # directly within the formula (which can't be left as NULL), using INLA notation!
        formula <- as.formula(gsub("\\(1 \\| ","f(",deparse(formula)))
      }
  }

  # Defines the name of the treatment variable (for which the main analysis is performed)
  if(exists("treat.name",where=exArgs)) {treatment <- exArgs$treat.name} else {treatment <- "treatment"}

  res <- list()
  tic <- proc.time()
  
  # For user-specified data generating processes
  if(exists("data",where=exArgs)) {
    func <- exArgs$data
    if(exists("inpts",where=exArgs)) {inpts=exArgs$inpts} else {inpts=list()}
    if (method=="INLA" || method=="inla") {
      # INLA-specific options --- can be used to set priors & stuff
      
      # 'control.fixed' can be used to set the priors for the fixed effects in the model
      if(exists("control.fixed",where=exArgs)) {
        control.fixed <- exArgs$control.fixed
      } else {
        control.fixed <- INLA::inla.set.control.fixed.default()
        # prior mean = 0 for *all* fixed effects
        # prior var = 1000 for *all* fixed effects
        # prior mean = 0 for the intercept
        # prior prec -> 0 for the intercept 
      }
      
      # 'offset' can be used with Poisson model. This needs to be a string with the name of the relevant variable
      if(exists("offset",where=exArgs)) {offset <- exArgs$offset} else {offset <- NULL}
      # 'Ntrials' is the Binomial denominator. This needs to be a string with the name of the relevant variable
      if(exists("Ntrials",where=exArgs)) {Ntrials <- exArgs$Ntrials} else {Ntrials <- NULL}

      res <- foreach::foreach(iterators::icount(n.sims), .combine=cbind, .packages="INLA", 
          .export=c("make.swt","sw.design.mat")) %dopar% {
        fake.data <- do.call(what=func,args=inpts)
        inla.args <- list(formula,data=fake.data,family=family,control.fixed=control.fixed)
        if (!is.null(offset)) {inla.args$offset <- fake.data[[offset]]}
        if (!is.null(Ntrials)) {inla.args$Ntrials <- fake.data[[Ntrials]]}
        m <- do.call(INLA::inla,args=inla.args)
        theta <- m$summary.fixed[treatment,1]
        sd.theta <- m$summary.fixed[treatment,2]
        ## This assumes that alpha=0.05 and so can consider the 95% CrI
        ci.fixed <- m$summary.fixed[treatment,c(3,5)]
        lower.lim.check <- ci.fixed[,1]>0
        upper.lim.check <- ci.fixed[,2]<0
        ## Needs to simulate values from the posterior of theta and then get the quantiles if alpha neq 0.05
        signif <- as.numeric(((sign(theta)==1 & lower.lim.check==TRUE) | (sign(theta)==-1 & upper.lim.check==TRUE)))
        tval <- theta / sd.theta
        pvalue <- 2*pnorm(abs(tval), lower.tail = FALSE)
        # Random effects sds
        mat <- INLA::inla.contrib.sd(m)$hyper
        rnd.eff.sd <- as.matrix(mat[,1])
        names(rnd.eff.sd) <- rownames(mat)
        list(power=signif,theta=theta,sd.theta=sd.theta,
          rnd.eff.sd=rnd.eff.sd,method=method,pvalue=pvalue) 
      }
    }
	
	  if (method=="lme") {
	    # Each 10% (default) of the iterations are stored in the results list, labelled res_1, res_2, etc.
	    res <- foreach::foreach(iterators::icount(n.sims), .combine=cbind, .packages="lme4", 
	     .export=c("make.swt","sw.design.mat")) %dopar% {
	    fake.data <- do.call(what=func,args=inpts)
	    # If the formula contains random effects, run lmer
	    check.random.effect <- !is.null(findbars(formula))
	    if(check.random.effect==TRUE) {
	      if(family=="gaussian") {
	        m <- lme4::lmer(formula, data=fake.data)          
	      } else {
          m <- lme4::glmer(formula, data=fake.data,family=family)
	      }
	      which.treat <- which(names(fixef(m))==treatment)
	      theta <- lme4::fixef(m)[which.treat]
	      sd.theta <- summary(m)$coefficients[which.treat,2]
	      Vcov <- vcov(m, useScale = FALSE)
	      se <- sqrt(diag(Vcov))
	      betas <- lme4::fixef(m)
	      tval <- betas / se
	      ### NB: issue with p-values in random effect models (cannot determine df easily)
	      ### An alternative (as recommended by Bates et al) is to use CIs instead of p-values
	      ci.fixed <- confint(m,names(fixef(m)),method="Wald",level=1-sig.level)
	      lower.lim.check <- ci.fixed[,1]>0
	      upper.lim.check <- ci.fixed[,2]<0
	      signif <- as.numeric(((sign(betas)==1 & lower.lim.check==TRUE) | (sign(betas)==-1 & upper.lim.check==TRUE))[2])
	      ### This estimates the p-value based on Normal approximation --- may not be too correct...
	      ### See: http://mindingthebrain.blogspot.co.uk/2014/02/three-ways-to-get-parameter-specific-p.html
	      pval <- 2*pnorm(abs(tval), lower.tail = FALSE)
	      pvalue <- pval[which.treat] 
	      ###signif <- pvalue < alpha
	      VC <- as.data.frame(lme4::VarCorr(m))
	      rnd.eff.sd <- VC[,which(colnames(VC)=="sdcor")]
	      # Finds the rows that needs to be reported
	      to.report <- c(which(is.na(VC[,3]==T), which(VC[,1]=="Residual")))
	      rnd.eff.sd <- VC[to.report,which(colnames(VC)=="sdcor")]
	      names(rnd.eff.sd)[which(is.na(VC[to.report,3])==T)] <- paste(VC[to.report,"grp"],VC[to.report,"var1"])
	      names(rnd.eff.sd) <- gsub(" NA","",names(rnd.eff.sd))
	      if(family=="gaussian") {method <- "lmer"} else {method <- "glmer"}
	    } 
	    
	    # If it doesn't then do lm
	    if (!check.random.effect) {
	       if (family=="gaussian") {
	         m <- lm(formula,data=fake.data)
	       } else {
	         m <- glm(formula, data=fake.data,family=family)
	       }
           which.treat <- which(names(m$coefficients)==treatment)
           theta <- m$coefficients[which.treat]
           sd.theta <- summary(m)$coefficients[which.treat,2]
           se <- summary(m)$coefficients[,2]
           betas <- summary(m)$coefficients[,1]
           tval <- summary(m)$coefficients[,3]
    	   df <- summary(m)$df[2]
    	   pval <- summary(m)$coefficients[,4]
    	   pvalue <- pval[which.treat] 
    	   signif <- pvalue < sig.level
    	   if(family=="gaussian") {method <- "lm"} else {method <- "glm"}
    	   rnd.eff.sd=NULL
	     }
         list(power=signif,theta=theta,sd.theta=sd.theta,
	       rnd.eff.sd=rnd.eff.sd,method=method,pvalue=pvalue) 
	     }
	  } 
    }

  # For standard SWT data generating processes (uses make.swt)
  if(!exists("data",where=exArgs)) {
    if(exists("X",where=exArgs)) {
      X=exArgs$X
      row.names(X) <- sample(1:I,I)
      colnames(X) <- c("Baseline",paste0("Time ",1:J))
    } else {
      X=NULL
    }
    
    if (method=="INLA" || method=="inla") {
      # INLA-specific options --- can be used to set priors & stuff
      
      # 'control.fixed' can be used to set the priors for the fixed effects in the model
      if(exists("control.fixed",where=exArgs)) {
        control.fixed <- exArgs$control.fixed
      } else {
        control.fixed <- INLA::inla.set.control.fixed.default()
        # prior mean = 0 for *all* fixed effects
        # prior var = 1000 for *all* fixed effects
        # prior mean = 0 for the intercept
        # prior prec -> 0 for the intercept 
      }
      
      # 'offset' can be used with Poisson model. This needs to be a string with the name of the relevant variable
      if(exists("offset",where=exArgs)) {offset <- exArgs$offset} else {offset <- NULL}
      # 'Ntrials' is the Binomial denominator. This needs to be a string with the name of the relevant variable
      if(exists("Ntrials",where=exArgs)) {Ntrials <- exArgs$Ntrials} else {Ntrials <- NULL}
      
      res <- foreach::foreach(iterators::icount(n.sims), .combine=cbind, .packages="INLA", 
        .export=c("make.swt","sw.design.mat")) %dopar% {
      fake.data <- fake.data <- make.swt(I=I,J=J,H=H,K=K,design=design,mu=mu,b.trt=b.trt,
                                         b.time=b.time,sigma.y=sigma.y,sigma.e=sigma.e,
                                         rho=rho,sigma.a=sigma.a,rho.ind=rho.ind,sigma.v=sigma.v,
                                         X=X,family=family,natural.scale=natural.scale)
      inla.args <- list(formula,data=fake.data,family=family,control.fixed=control.fixed)
      if (!is.null(offset)) {inla.args$offset <- fake.data[[offset]]}
      if (!is.null(Ntrials)) {inla.args$Ntrials <- fake.data[[Ntrials]]}
      m <- do.call(INLA::inla,args=inla.args)
      theta <- m$summary.fixed[treatment,1]
      sd.theta <- m$summary.fixed[treatment,2]
      ## This assumes that alpha=0.05 and so can consider the 95% CrI
      ci.fixed <- m$summary.fixed[treatment,c(3,5)]
      lower.lim.check <- ci.fixed[,1]>0
      upper.lim.check <- ci.fixed[,2]<0
      ## Needs to simulate values from the posterior of theta and then get the quantiles if alpha neq 0.05
      signif <- as.numeric(((sign(theta)==1 & lower.lim.check==TRUE) | (sign(theta)==-1 & upper.lim.check==TRUE)))
      tval <- theta / sd.theta
      pvalue <- 2*pnorm(abs(tval), lower.tail = FALSE)
      # Random effects sds
      mat <- INLA::inla.contrib.sd(m)$hyper
      rnd.eff.sd <- as.matrix(mat[,1])
      names(rnd.eff.sd) <- rownames(mat)
      list(power=signif,theta=theta,sd.theta=sd.theta,
          rnd.eff.sd=rnd.eff.sd,method=method,pvalue=pvalue) 
      }
    }
    
    if (method=="lme") {

      res <- foreach::foreach(iterators::icount(n.sims), .combine=cbind, .packages="lme4", 
         .export=c("make.swt","sw.design.mat")) %dopar% {
	
        fake.data <- make.swt(I=I,J=J,H=H,K=K,design=design,mu=mu,b.trt=b.trt,
                            b.time=b.time,sigma.y=sigma.y,sigma.e=sigma.e,
                            rho=rho,sigma.a=sigma.a,rho.ind=rho.ind,sigma.v=sigma.v,
                            X=X,family=family,natural.scale=natural.scale)

        # If the formula contains random effects, run lmer
        check.random.effect <- !is.null(findbars(formula))
        if(check.random.effect==TRUE) {
          if(family=="gaussian") {
            m <- lme4::lmer(formula, data=fake.data)          
          } else {
            m <- lme4::glmer(formula, data=fake.data,family=family)
          }
          which.treat <- which(names(fixef(m))==treatment)
          theta <- lme4::fixef(m)[which.treat]
          sd.theta <- summary(m)$coefficients[which.treat,2]
          Vcov <- vcov(m, useScale = FALSE)
          se <- sqrt(diag(Vcov))
          betas <- lme4::fixef(m)
          tval <- betas / se
          ### NB: issue with p-values in random effect models (cannot determine df easily)
          ### An alternative (as recommended by Bates et al) is to use CIs instead of p-values
          ci.fixed <- confint(m,names(fixef(m)),method="Wald",level=1-sig.level)
          lower.lim.check <- ci.fixed[,1]>0
          upper.lim.check <- ci.fixed[,2]<0
          signif <- as.numeric(((sign(betas)==1 & lower.lim.check==TRUE) | (sign(betas)==-1 & upper.lim.check==TRUE))[2])
          ### This estimates the p-value based on Normal approximation --- may not be too correct...
          ### See: http://mindingthebrain.blogspot.co.uk/2014/02/three-ways-to-get-parameter-specific-p.html
          pval <- 2*pnorm(abs(tval), lower.tail = FALSE)
          pvalue <- pval[which.treat] 
          ###signif <- pvalue < alpha
          VC <- as.data.frame(lme4::VarCorr(m))
          rnd.eff.sd <- VC[,which(colnames(VC)=="sdcor")]
          # Finds the rows that needs to be reported
          to.report <- c(which(is.na(VC[,3]==T), which(VC[,1]=="Residual")))
          rnd.eff.sd <- VC[to.report,which(colnames(VC)=="sdcor")]
          names(rnd.eff.sd)[which(is.na(VC[to.report,3])==T)] <- paste(VC[to.report,"grp"],VC[to.report,"var1"])
          names(rnd.eff.sd) <- gsub(" NA","",names(rnd.eff.sd))
          if(family=="gaussian") {method <- "lmer"} else {method <- "glmer"}
        }
        # If it doesn't then do lm
        if (!check.random.effect) {
          if (family=="gaussian") {
            m <- lm(formula,data=fake.data)
          } else {
            m <- glm(formula, data=fake.data,family=family)
          }
          which.treat <- which(names(m$coefficients)==treatment)
          theta <- m$coefficients[which.treat]
          sd.theta <- summary(m)$coefficients[which.treat,2]
          se <- summary(m)$coefficients[,2]
          betas <- summary(m)$coefficients[,1]
          tval <- summary(m)$coefficients[,3]
          df <- summary(m)$df[2]
          pval <- summary(m)$coefficients[,4]
          pvalue <- pval[which.treat] 
          signif <- pvalue < sig.level
          if(family=="gaussian") {method <- "lm"} else {method <- "glm"}
          rnd.eff.sd=NULL
        }
        list(power=signif,theta=theta,sd.theta=sd.theta,
           rnd.eff.sd=rnd.eff.sd,method=method,pvalue=pvalue) 
      }
    }
  }

  toc <- proc.time(); time2run <- (toc-tic)[3]; names(time2run) <- "Time to run (secs)"
  theta=mean(unlist(res[2,]),na.rm=T); names(theta)=NULL
  sd.theta <- mean(unlist(res[3,]),na.rm=T); names(sd.theta)=NULL
  theta=c(theta,sd.theta); names(theta)=c("Estimate","Standard Error")
  pvalue=unlist(res[6,]); names(pvalue)=NULL
  method=unlist(res[5,1]); names(method)=NULL
  power=mean(unlist(res[1,]),na.rm=T)
  ci <- power+c(-qnorm(1-sig.level/2),qnorm(1-sig.level/2))*sqrt(var(pvalue<sig.level,na.rm=T)/n.sims)
  ci <- power+c(-qnorm(1-sig.level/2),qnorm(1-sig.level/2))*sqrt(var(unlist(res[1,]),na.rm=T)/n.sims)
  if(ci[1]<0) {ci[1] <- 0}; if(ci[2]>1) {ci[2] <- 1}    # Constrains ci in [0;1]
  if (!exists("data",where=exArgs)) {
    setting <- list(n.clusters=I,n.time.points=J,avg.cluster.size=K,
                    design=design,formula=formula,method=method,family=family)
  } else {setting=list(formula=formula,method=method,family=family)}
  closeAllConnections()
  
  # If the option 'plot=T' creates a moving average plot of the power
  if (plot==TRUE) {
    # Defines some cutoff points in the simulations to plot the moving average
    percs <- seq(.1,1,by=.1)
    # If the first %tile of the simulations is 0, then sets to 1
    steps <- pmax(round(n.sims*percs),1)
    # Computes the power moving average
    mov_av <- numeric()
    for (s in 1:length(steps)) {
      mov_av[s] <- mean(unlist(res[1,c(1:steps[s])]))
    }
    # plot the moving average of power in a graph
    plot(mov_av, type="o", ylim=c(0,1), ylab="Power", xlab="Percentage of Iterations", xaxt = 'n', main="Moving Average of Power")
    axis(1, at=1:10, labels = c(paste((1:10)*10, "%", sep="")))
    #### NEED TO CHECK THIS
    abline(h=power,col="red",lwd=2)
    abline(h=ci[1],col="red",lwd=2,lty=2)
    abline(h=ci[2],col="red",lwd=2,lty=2)
  }
  if (!is.null(unlist(res[4,]))) {
      rnd.eff.sd=apply(do.call(cbind.data.frame,res[4,]),1,mean,na.rm=T)
  } else {
      rnd.eff.sd=NULL
  }
  list(power=power,time2run=time2run,ci.power=ci,theta=theta,
       rnd.eff.sd=rnd.eff.sd,setting=setting) # pvalue=pvalue,
}

#
## OLD VERSION WHICH DOESN'T WORK ON WINDOWS
# sim.swt <- function (I,J,H=NULL,K,design="cross-sec",mu=0,b.trt,b.time=NULL,
#                          sigma.y=NULL,sigma.e=NULL,rho=NULL,sigma.a=NULL,
#                          rho.ind=NULL,sigma.v=NULL,n.sims=1000,formula=NULL,
#                          alpha=0.05,n.cores=NULL,...) {
#   
#   ## Power analysis for repeated cohort design (cross-sectional) - continuous outcome
#   ## The optional arguments include the specification of a user-defined function data
#   ## which is used to simulate the data, according to any design the user has.
#   ## So including data=data, where data=function(){...} specifies the simulation 
#   ## of a suitable dataset allows the user to select *any* kind of model for data
#   ## generation. If data has inputs, then these must be specified as a list 
#   ## inpts, to be passed as an extra argument to sim.sw.cont. 
#   ## In addition to this, the user should then select a suitable formula
#   ## which is in line with the model used for the simulation of the virtual trial.
#   ## Also, the name of the "treatment" variable can be specified as a string in the
#   ## extra argument treat.name. If not present then it is assumed that the name of the
#   ## treatment variable is, unsurprisingly, "treatment"
#   
#   exArgs <- list(...)
# 
#   ## Uses parallel computation to speed up task
# #   required <- c("lme4","foreach","doParallel","parallel","iterators")
# #   for (i in 1:length(required)) {
# #       txt1 <- paste0("if(!isTRUE(requireNamespace('",required[i],"',quietly=TRUE))) {")
# #       txt2 <- paste0("stop('You need to install the package ",required[i],". Please run in your R terminal:
# # install.packages(",required[i],")')")
# #       txt3 <- paste0("}")
# #       txt <- paste(txt1,txt2,txt3,sep="\n",collapse=" ")
# #       eval(parse(text=txt))
# #   }
#   requireNamespace("lme4",quietly=TRUE)
#   requireNamespace("foreach",quietly=TRUE)
#   requireNamespace("doParallel",quietly=TRUE)
#   requireNamespace("parallel",quietly=TRUE)
#   requireNamespace("iterators",quietly=TRUE)
#   
#   # If not specified, uses all cores available to the machine
#   if(is.null(n.cores)) {n.cores <- parallel::detectCores()}
#   # Usually a good idea to leave 1 core free for the OS to use
#   doParallel::registerDoParallel(cores=n.cores-1)
#   
#   # Defines basic formulation for the model
#   if(is.null(formula)){
#     formula=y~treatment+factor(time)+(1|cluster)
#     if(design=="cohort") {formula <- update.formula(formula,.~.+(1|id))}
#   }
#   # Defines the name of the treatment variable (for which the main analysis is performed)
#   if(exists("treat.name",where=exArgs)) {treatment <- exArgs$treat.name} else {treatment <- "treatment"}
# 
#   # Simulates the datasets & analysis
#   tic <- proc.time()
#   ### NB Needs to use the .export argument to the function foreach, to include
#   ### variables that are passed outside of the scope (eg b.trt) --- need to figure
#   ### this out a bit better. But basically, linux/mac will make stuff available,
#   ### while windows won't and so throw an error because some variable that is input
#   ### the function in foreach is not loaded in memory...
#   res <- foreach::foreach(iterators::icount(n.sims), .combine=cbind, .packages="lme4", 
# 	.export=c("make.swt","sw.design.mat")) %dopar% {
#     if(exists("data",where=exArgs)) {
#       func <- exArgs$data
#       if(exists("inpts",where=exArgs)) {inpts=exArgs$inpts} else {inpts=list()}
#       fake.data <- do.call(what=func,args=inpts)
#     } 
#     if(!exists("data",where=exArgs)) {
#       if(exists("X",where=exArgs)) {
#         X=exArgs$X
#         row.names(X) <- sample(1:I,I)
#         colnames(X) <- c("Baseline",paste0("Time ",1:J))
#       } else {
#         X=NULL
#       }
#       fake.data <- make.swt(I=I,J=J,H=H,K=K,design=design,mu=mu,b.trt=b.trt,
#                                    b.time=b.time,sigma.y=sigma.y,sigma.e=sigma.e,
#                                    rho=rho,sigma.a=sigma.a,rho.ind=rho.ind,
#                                    sigma.v=sigma.v,X=X)
#     }
#     
#     # If the formula contains random effects, run lmer
#     check.random.effect <- !is.null(findbars(formula))
#      if(check.random.effect==TRUE) {
#         m <- lme4::lmer(formula, data=fake.data)
#         which.treat <- which(names(fixef(m))==treatment)
#         theta <- lme4::fixef(m)[which.treat]
#         sd.theta <- summary(m)$coefficients[which.treat,2]
#         Vcov <- vcov(m, useScale = FALSE)
#         se <- sqrt(diag(Vcov))
#         betas <- lme4::fixef(m)
#         tval <- betas / se
#         ### NB: issue with p-values in random effect models (cannot determine df easily)
#         ### An alternative (as recommended by Bates et al) is to use CIs instead of p-values
#         ci.fixed <- confint(m,names(fixef(m)),method="Wald",level=1-alpha)
#         lower.lim.check <- ci.fixed[,1]>0
#         upper.lim.check <- ci.fixed[,2]<0
#         signif <- as.numeric(((sign(betas)==1 & lower.lim.check==TRUE) | (sign(betas)==-1 & upper.lim.check==TRUE))[2])
#         ### This estimates the p-value based on Normal approximation --- may not be too correct...
#         ### See: http://mindingthebrain.blogspot.co.uk/2014/02/three-ways-to-get-parameter-specific-p.html
#         pval <- 2*pnorm(abs(tval), lower.tail = FALSE)
#         pvalue <- pval[which.treat] 
#         ###signif <- pvalue < alpha
#         VC <- as.data.frame(lme4::VarCorr(m))
#         rnd.eff.sd <- VC[,which(colnames(VC)=="sdcor")]
#         # Finds the rows that needs to be reported
#         to.report <- c(which(is.na(VC[,3]==T), which(VC[,1]=="Residual")))
#         rnd.eff.sd <- VC[to.report,which(colnames(VC)=="sdcor")]
#         names(rnd.eff.sd)[which(is.na(VC[to.report,3])==T)] <- paste(VC[to.report,"grp"],VC[to.report,"var1"])
#         names(rnd.eff.sd) <- gsub(" NA","",names(rnd.eff.sd))
#         method <- "lmer"
#      }
#      # If it doesn't then do lm
#      if (!check.random.effect) {
#          m <- lm(formula,data=fake.data)
#          which.treat <- which(names(m$coefficients)==treatment)
#          theta <- m$coefficients[which.treat]
#          sd.theta <- summary(m)$coefficients[which.treat,2]
#          se <- summary(m)$coefficients[,2]
#          betas <- summary(m)$coefficients[,1]
#          tval <- summary(m)$coefficients[,3]
#          df <- summary(m)$df[2]
#  ##        pval <- 2*pt(abs(tval),df=df,lower.tail = FALSE)
#          pval <- summary(m)$coefficients[,4]
#          pvalue <- pval[which.treat] 
#          signif <- pvalue < alpha
#          method <- "lm"
#          rnd.eff.sd=NULL
#      }
#      list(power=signif,theta=theta,sd.theta=sd.theta,rnd.eff.sd=rnd.eff.sd,method=method,pvalue=pvalue) 
#   }
#   toc <- proc.time(); time2run <- (toc-tic)[3]; names(time2run) <- "Time to run (secs)"
#   theta=mean(unlist(res[2,]),na.rm=T); names(theta)=NULL
#   sd.theta <- mean(unlist(res[3,]),na.rm=T); names(sd.theta)=NULL
#   theta=c(theta,sd.theta); names(theta)=c("Estimate","Standard Error")
#   pvalue=unlist(res[6,]); names(pvalue)=NULL
#   power=mean(unlist(res[1,]),na.rm=T)
#   ci <- power+c(-qnorm(1-alpha/2),qnorm(1-alpha/2))*sqrt(var(pvalue<alpha,na.rm=T)/n.sims)
#   if(ci[1]<0) {ci[1] <- 0}; if(ci[2]>1) {ci[2] <- 1}    # Constrains ci in [0;1]
#   if (!exists("data",where=exArgs)) {
#     setting <- list(n.clusters=I,n.time.points=J,avg.cluster.size=K,
#                     design=design,formula=formula)
#   } else {setting=list(formula=formula)}
#   
#   list(power=power,time2run=time2run,ci.power=ci,theta=theta,
#        rnd.eff.sd=unlist(res[4,1]),setting=setting) # pvalue=pvalue,
# }






#' Power calculation for binary outcome based on analytic formula of Hussey and
#' Hughes
#' 
#' Sample size calculations for binary outcomes based on the formula provided
#' by Hussey and Hughes (2007)
#' 
#' 
#' @param p1 Baseline probability of the outcome (for the controls)
#' @param OR Value of the expected Odds Ratio (for the intervention vs control)
#' @param I Number of clusters
#' @param J Number of time points
#' @param K Average size of each cluster
#' @param rho Intra-class correlation coefficient (default=0)
#' @param sig.level Significance level (default=0.05)
#' @param which.var String character specifying which variance to be considered
#' (options are the default value \code{'within'} or \code{'total'}
#' @param X A design matrix for the stepped wedge design, indicating the time
#' at which each of the clusters should switch the active intervention. By
#' default is NULL and automatically computed, but can be passed as an extra
#' argument as a user-defined matrix with I rows and (J+1) columns
#' @return \item{power}{ The resulting power } \item{sigma.y}{The estimated
#' total (marginal) sd for the outcome} \item{sigma.e}{The estimated residual
#' sd} \item{sigma.a}{The resulting cluster-level sd} \item{setting}{A list
#' including the following values: - n.clusters = The number of clusters -
#' n.time.points = The number of 'active' time points - avg.cluster.size = The
#' average cluster size - design.matrix = The design matrix for the SWT under
#' consideration }
#' @author Gianluca Baio
#' @references Baio, G; Copas, A; Ambler, G; Hargreaves, J; Beard, E; and Omar,
#' RZ Sample size calculation for a stepped wedge trial. Trials, 16:354. Aug
#' 2015.
#' 
#' Hussey M and Hughes J. Design and analysis of stepped wedge cluster
#' randomized trials. Contemporary Clinical Trials. 28(2):182-91. Epub 2006 Jul
#' 7. Feb 2007
#' @keywords Hussey and Hughes formula
#' @examples
#' 
#' HH.binary(p1=.26,OR=.55,I=10,J=5,K=20,rho=.2)
#' @export HH.binary
HH.binary <- function(p1,OR,I,J,K,rho=0,sig.level=0.05,which.var="within",X=NULL) {
  # HH sample size calculations for binary outcome
  if(is.null(X)) {
    X <- sw.design.mat(I=I,J=J,H=NULL)
  } else {
    row.names(X) <- sample(1:I,I)
    colnames(X) <- c("Baseline",paste0("Time ",1:J))
  }
  U <- sum(X)
  W <- sum(apply(X,2,sum)^2)
  V <- sum(apply(X,1,sum)^2)
  
  # Data
  p2 <-  OR*(p1/(1-p1))/(1+OR*(p1/(1-p1)))
  theta=abs(p1-p2)
  if (which.var=="within") {
    sigma.e=sqrt((p1*(1-p1)+p2*(1-p2))/2)    
    ####sigma.e=sqrt((p1*(1-p1)))                     # Consistent with Hemming - stata
    sigma.a <- sqrt(rho*sigma.e^2/(1-rho))
    sigma.y <- sqrt(sigma.e^2+sigma.a^2)
  }
  if (which.var=="total") {
    sigma.y=sqrt((p1*(1-p1)+p2*(1-p2))/2)    
    ####sigma.y=sqrt((p1*(1-p1)))
    sigma.a <- sqrt(sigma.y^2*rho)
    sigma.e <- sqrt(sigma.y^2-sigma.a^2)
  }
  sigma <- sqrt(sigma.e^2/K)      
  
  # Power calculations
  v <- (I*sigma^2*(sigma^2+((J+1)*sigma.a^2)))/((I*U-W)*sigma^2+(U^2+I*(J+1)*U-(J+1)*W-I*V)*sigma.a^2)
  power <- pnorm(theta/sqrt(v)-qnorm(1-sig.level/2))
  setting <- list(n.clusters=I,n.time.points=J,avg.cluster.size=K,design.matrix=X)
  list(power=power,p1=p1,p2=p2,sigma.y=sigma.y,sigma.e=sigma.e,sigma.a=sigma.a,setting=setting)
}






#' Power calculation for normal outcome based on analytic formula of Hussey and
#' Hughes
#' 
#' Sample size calculations for normal outcomes based on the formula provided
#' by Hussey and Hughes (2007)
#' 
#' 
#' @param mu Mean value of the outcome for the controls
#' @param b.trt Treatment effect against controls
#' @param sigma Value of the standard deviation (if \code{which.var}='within'
#' then it's assumed to be the residual sd. If \code{which.var}='total', then
#' it's assumed to be the total sd)
#' @param I Number of clusters
#' @param J Number of time points
#' @param K Average size of each cluster
#' @param rho Intra-class correlation coefficient (default=0)
#' @param sig.level Significance level (default=0.05)
#' @param which.var String character specifying which variance to be considered
#' (options are the default value \code{'within'} or \code{'total'}
#' @param X A design matrix for the stepped wedge design, indicating the time
#' at which each of the clusters should switch the active intervention. By
#' default is NULL and automatically computed, but can be passed as an extra
#' argument as a user-defined matrix with I rows and (J+1) columns
#' @return \item{power}{ The resulting power } \item{sigma.y}{The estimated
#' total (marginal) sd for the outcome} \item{sigma.e}{The estimated residual
#' sd} \item{sigma.a}{The resulting cluster-level sd} \item{setting}{A list
#' including the following values: - n.clusters = The number of clusters -
#' n.time.points = The number of 'active' time points - avg.cluster.size = The
#' average cluster size - design.matrix = The design matrix for the SWT under
#' consideration }
#' @author Gianluca Baio
#' @references Baio, G; Copas, A; Ambler, G; Hargreaves, J; Beard, E; and Omar,
#' RZ Sample size calculation for a stepped wedge trial. Trials, 16:354. Aug
#' 2015.
#' 
#' Hussey M and Hughes J. Design and analysis of stepped wedge cluster
#' randomized trials. Contemporary Clinical Trials. 28(2):182-91. Epub 2006 Jul
#' 7. Feb 2007
#' @keywords Hussey and Hughes formula
#' @examples
#' 
#' HH.normal(mu=.3,b.trt=-.3875,I=10,J=5,K=20,rho=.2,sigma=1.55)
#' 
#' @export HH.normal
HH.normal <- function(mu,b.trt,sigma,I,J,K,rho=0,sig.level=0.05,which.var="within",X=NULL) {
  # HH sample size calculations for continuous (normal) outcome
  # Stepped wedge design matrix
  if(is.null(X)) {
    X <- sw.design.mat(I=I,J=J,H=NULL)
  } else {
    row.names(X) <- sample(1:I,I)
    colnames(X) <- c("Baseline",paste0("Time ",1:J))
  }
  U <- sum(X)
  W <- sum(apply(X,2,sum)^2)
  V <- sum(apply(X,1,sum)^2)
  
  # Data
  mu1 <- mu+b.trt
  theta <- abs(mu1-mu)
  # Assumes that the input sigma.y is in fact the within cluster sd
  if (which.var=="within") {
    sigma.e <- sigma
    sigma.a <- sqrt(rho*sigma.e^2/(1-rho))
    sigma.y <- sqrt(sigma.e^2+sigma.a^2)
  }
  # Assumes that the input sigma.y is total sd
  if (which.var=="total") {    
    sigma.a <- sqrt(sigma^2*rho)
    sigma.e <- sqrt(sigma^2-sigma.a^2)
    sigma.y <- sigma
  }
  sigma <- sqrt(sigma.e^2/K)
  
  # Power calculations
  v <- (I*sigma^2*(sigma^2+((J+1)*sigma.a^2)))/((I*U-W)*sigma^2+(U^2+I*(J+1)*U-(J+1)*W-I*V)*sigma.a^2)
  power <- pnorm(theta/sqrt(v)-qnorm(1-sig.level/2))
  setting <- list(n.clusters=I,n.time.points=J,avg.cluster.size=K,design.matrix=X)
  list(power=power,sigma.y=sigma.y,sigma.e=sigma.e,sigma.a=sigma.a,setting=setting)
}






#' Power calculation for count outcome based on analytic formula of Hussey and
#' Hughes
#' 
#' Sample size calculations for count outcomes based on the formula provided by
#' Hussey and Hughes (2007)
#' 
#' 
#' @param lambda1 Baseline value for the rate at which the outcome occurs
#' @param RR Relative risk (of the intervention vs the control)
#' @param I Number of clusters
#' @param J Number of time points
#' @param K Average size of each cluster
#' @param rho Intra-class correlation coefficient (default=0)
#' @param sig.level Significance level (default=0.05)
#' @param which.var String character specifying which variance to be considered
#' (options are the default value \code{'within'} or \code{'total'}
#' @param X A design matrix for the stepped wedge design, indicating the time
#' at which each of the clusters should switch the active intervention. By
#' default is NULL and automatically computed, but can be passed as an extra
#' argument as a user-defined matrix with I rows and (J+1) columns
#' @return \item{power}{ The resulting power } \item{sigma.y}{The estimated
#' total (marginal) sd for the outcome} \item{sigma.e}{The estimated residual
#' sd} \item{sigma.a}{The resulting cluster-level sd} \item{setting}{A list
#' including the following values: - n.clusters = The number of clusters -
#' n.time.points = The number of 'active' time points - avg.cluster.size = The
#' average cluster size - design.matrix = The design matrix for the SWT under
#' consideration }
#' @author Gianluca Baio
#' @references Baio, G; Copas, A; Ambler, G; Hargreaves, J; Beard, E; and Omar,
#' RZ Sample size calculation for a stepped wedge trial. Trials, 16:354. Aug
#' 2015.
#' 
#' Hussey M and Hughes J. Design and analysis of stepped wedge cluster
#' randomized trials. Contemporary Clinical Trials. 28(2):182-91. Epub 2006 Jul
#' 7. Feb 2007
#' @keywords Hussey and Hughes formula
#' @examples
#' 
#' HH.count(lambda1=1.55,RR=.87,I=10,J=5,K=20,rho=.2)
#' 
#' @export HH.count
HH.count <- function(lambda1,RR,I,J,K,rho=0,sig.level=0.05,which.var="within",X=NULL) {
  # HH sample size calculations for continuous (normal) outcome
  # Stepped wedge design matrix
  if(is.null(X)) {
    X <- sw.design.mat(I=I,J=J,H=NULL)
  } else {
    row.names(X) <- sample(1:I,I)
    colnames(X) <- c("Baseline",paste0("Time ",1:J))
  }
  U <- sum(X)
  W <- sum(apply(X,2,sum)^2)
  V <- sum(apply(X,1,sum)^2)
  
  # Data
  lambda2 <- RR*lambda1
  theta <- abs(lambda1-lambda2)
  if (which.var=="within") {
    sigma.e=(sqrt(lambda1)+sqrt(lambda2))/2       # Consistent with Hemming - stata
    sigma.a <- sqrt(rho*sigma.e^2/(1-rho))
    sigma.y <- sqrt(sigma.e^2+sigma.a^2)
  }
  if (which.var=="total") {
    sigma.y <- (sqrt(lambda1)+sqrt(lambda2))/2    # Consisten with Hemming - stata
    sigma.a <- sqrt(sigma.y^2*rho)
    sigma.e <- sqrt(sigma.y^2-sigma.a^2)
  }
  sigma <- sqrt(sigma.e^2/K)  
  
  # Power calculations
  v <- (I*sigma^2*(sigma^2+((J+1)*sigma.a^2)))/((I*U-W)*sigma^2+(U^2+I*(J+1)*U-(J+1)*W-I*V)*sigma.a^2)
  power <- pnorm(theta/sqrt(v)-qnorm(1-sig.level/2))
  setting <- list(n.clusters=I,n.time.points=J,avg.cluster.size=K,design.matrix=X)
  list(power=power,lambda1=lambda1,lambda2=lambda2,sigma.y=sigma.y,sigma.e=sigma.e,sigma.a=sigma.a,setting=setting)
}






#' Creates a design matrix for a Stepped Wedge Trial
#' 
#' Constructs a basic SWT design matrix
#' 
#' 
#' @param I Number of clusters
#' @param J Number of time points
#' @param H Number of units randomised at each time point
#' @return Returns a design matrix X
#' @author Gianluca Baio
#' @references Baio, G; Copas, A; Ambler, G; Hargreaves, J; Beard, E; and Omar,
#' RZ Sample size calculation for a stepped wedge trial. Trials, 16:354. Aug
#' 2015.
#' @export sw.design.mat
sw.design.mat <- function(I,J,H=NULL) {
  ## Creates the design matrix for a stepped wedge design
  ## Checks to see that data are consistent with SW design
  if(sum(sapply(list(I,J,H), is.null)) != 1) {
    warning("exactly one of 'I', 'J' and 'H' must be NULL")
  }
  if (is.null(I)) {
    I <- H*J
  }
  if (is.null(J)) {
    J <- I/H
  }
  if (is.null(H)) {
    H <- I/J
  }
  
  # Stepped wedge design matrix
  X <- matrix(0,I,(J+1))            
  for (i in 2:(J+1)) {
    X[1:((i-1)*H),i] <- 1
  }
  row.names(X) <- sample(1:I,I)
  colnames(X) <- c("Baseline",paste0("Time ",1:J))
  return(X)
}






#' Computes the Design Effect for a Stepped Wedge Trial
#' 
#' Sample size calculations for a SWT using a cross-sectional design.  This is
#' based on (the correct version) of Woertman et al (2013), as described in
#' Baio et al (2015).
#' 
#' 
#' @param outcome String. Type of outcome. Options are \code{cont}, \code{bin}
#' or \code{count}
#' @param input input = a list containing the arguments. This differs depending
#' on the type of outcome, as follows: - continuous outcome: 1) delta
#' (treatment effect) 2) sd (standard deviation) - binary outcome: 1) p1
#' (baseline probability of outcome) 2) either p2 (treatment probability of
#' outcome), or OR (treatment effect as OR) - count outcome: 1) r1 (baseline
#' rate of outcome) 2) either r2 (treatment rate of outcome), or RR (treatment
#' effect as RR)
#' @param K average cluster size
#' @param J number of time points (excluding baseline)
#' @param B number of baseline measurement times
#' @param T number of measurement times during each crossover
#' @param rho ICC
#' @param sig.level significance level (default = 0.05)
#' @param power Power (default = 0.8)
#' @return \item{n.cls.swt}{ Number of clusters required to reach the
#' pre-specified power with the given significance level. } \item{n.pts}{ The
#' total number of participants required. } \item{DE.woert}{ The resulting
#' Design Effect. } \item{CF}{ The resulting Correction Factor. } \item{n.rct}{
#' The original individual RCT sample required to reach the pre-specified power
#' with the given significance level. }
#' @author Gianluca Baio
#' @references Baio, G; Copas, A; Ambler, G; Hargreaves, J; Beard, E; and Omar,
#' RZ Sample size calculation for a stepped wedge trial. Trials, 16:354. Aug
#' 2015.
#' @keywords Design Effect
#' @examples
#' 
#' # Continuous outcome
#' input <- list(delta=-0.3875,sd=1.55)
#' K <- 20
#' J <- 5
#' rho <- .2
#' DE.woert(input=input,K=K,J=J,rho=rho)
#' #
#' # Binary outcome
#' input <- list(OR=.53,p1=.26)
#' DE.woert(outcome="bin",input=input,K=K,J=J,rho=rho)
#' #
#' # Count outcome
#' input <- list(RR=.8,r1=1.5)
#' DE.woert(outcome="count",input=input,K=K,J=J,rho=rho)
#' 
#' @export DE.woert
DE.woert <- function(outcome="cont",input,K,J,B=1,T=1,rho,sig.level=0.05,power=.8) {
  # Computes the power for a SWT using the corrected form of the Woertman Design Effect
  # outcome = a string (default "cont", possible values are "bin" or "count")
  # input = a list containing the arguments
  #   - continuous outcome: 
  #      1) delta (treatment effect)
  #      2) sd (standard deviation)
  #   - binary outcome:
  #      1) p1 (baseline probability of outcome)
  #      2) either p2 (treatment probability of outcome), or OR (treatment effect as OR)
  #   - count outcome: 
  #      1) r1 (baseline rate of outcome)
  #      2) either r2 (treatment rate of outcome), or RR (treatment effect as RR)
  # K = average cluster size
  # J = number of time points (excluding baseline)
  # B = number of baseline measurement times
  # T = number of measurement times during each crossover
  # rho = ICC
  # alpha = significance level (default = 0.05)
  # beta = type 2 error (default = 0.2, implying power = 0.8)
  
  if(outcome=="cont") {
    n.rct <- 2*ceiling(power.t.test(delta=input$delta,sd=input$sd,
                                    sig.level=sig.level,power=power)$n)
  }
  if(outcome=="bin") {
    if(!is.null(input$OR) & is.null(input$p2)) {
        p1 <- input$p1
        p2 <- input$OR*(input$p1/(1-input$p1))/(1+input$OR*(input$p1/(1-input$p1)))
        OR <- input$OR
    }
    n.rct <- 2*ceiling(power.prop.test(p1=p1,p2=p2,sig.level=sig.level,power=power)$n)
  }
  if(outcome=="count") {
    if(!is.null(input$RR) & is.null(input$r2)) {
      r1 <- input$r1
      RR <- input$RR
      r2 <- input$r1*input$RR
    }
    if(is.null(input$RR) & !is.null(input$r2)) {
      r1 <- input$r1
      r2 <- input$r2
      RR <- input$r2/input$r1
    }
    n.rct <- 2*ceiling((((r1*(1+RR))*(qnorm((1-sig.level/2))+(qnorm(power)))^2)/(r1-r1*RR)^2))
  }
  DE.crt <- 1+(K*(J+1)-1)*rho
  n.cls <- ceiling((n.rct*DE.crt)/(K*(J+1)))
  CF <- (1+rho*(J*T*K+B*K-1))/(1+rho*(.5*J*T*K+B*K-1))*(3*(1-rho))/(2*T*(J-(1/J)))
  DE.woert <- (B+J*T)*CF
  
  n.woert <- n.rct*DE.woert
  n.cls.swt <- ceiling(n.woert/(K*(J+1)))
  
  list(n.cls.swt=n.cls.swt,n.pts.swt=n.woert,DE.woert=DE.woert,CF=CF,n.rct=n.rct)
}








#' Finds the optimum number of clusters or time points
#' 
#' Given some inputs determines the optimal combination of clusters/time points
#' to get a set level of power.
#' 
#' 
#' @param target.power The target power (eg 0.8)
#' @param I A vector specifying the range in which to search for the optimal
#' number of clusters, eg \code{I=c(1,10)}
#' @param J Number of time points
#' @param H Number of units randomised at each time point
#' @param K Average size of each cluster
#' @param design type of design. Can be \code{'cross-sec'} (default) or
#' \code{'cohort'} (repeated measurements)
#' @param mu baseline outcome value
#' @param b.trt Treatment effect
#' @param b.time Time effect
#' @param sigma.y total standard deviation
#' @param sigma.e individual standard deviation
#' @param rho Intra-class correlation coefficient
#' @param sigma.a the sd of the the cluster-level intercept (default at NULL)
#' @param rho.ind individual-level ICC (for cohorts)
#' @param sigma.v the sd of the cluster-level slope (by intervention, default
#' at NULL)
#' @param n.sims Number of simulations to be used (default=1000)
#' @param formula Formula describing the model to be used
#' @param family The model family to be used. Default value is 'gaussian' and
#' other possibile choices are 'binomial' or 'poisson'
#' @param natural.scale Indicator for whether the input is passed on the
#' natural scale or on the scale of the linear predictor. By default is set to
#' TRUE. In the case of family='gaussian' it does not have any effect, since
#' the link for the linear predictor is the identity. But for family='binomial'
#' or family='poisson', the user has to specify when the input is given on the
#' logit or log scale
#' @param sig.level Significance level (default=0.05)
#' @param n.cores Specifies the number of processors to be used for the
#' computation (default=NULL, which means that R will try and figure out)
#' @param ... Additional arguments
#' @return \item{Optimum_I}{ The value of the optimal number of clusters }
#' \item{power}{The estimated power in correspondence of the optimal I}
#' \item{time2run}{Computational time}
#' @author Rosie Leach
#' @references Baio, G; Copas, A; Ambler, G; Hargreaves, J; Beard, E; and Omar,
#' RZ Sample size calculation for a stepped wedge trial. Trials, 16:354. Aug
#' 2015.
#' 
#' Hussey M and Hughes J. Design and analysis of stepped wedge cluster
#' randomized trials. Contemporary Clinical Trials. 28(2):182-91. Epub 2006 Jul
#' 7. Feb 2007
#' @keywords Stepped wedge design
#' @examples
#' 
#' #cluster.search(I=c(4,10),target.power=.8,J=6,K=30,mu=1.5,b.trt=.8,rho=0,
#' #family="poisson",n.sims=10)
#' 
#' @export cluster.search
cluster.search <- function(target.power=NULL, I=NULL, J=NULL ,H=NULL,K,design="cross-sec",mu=0,b.trt,b.time=NULL,
                      sigma.y=NULL,sigma.e=NULL,rho=NULL,sigma.a=NULL,
                      rho.ind=NULL,sigma.v=NULL,n.sims=1000,formula=NULL,
                      family="gaussian",natural.scale=TRUE,sig.level=0.05,n.cores=NULL,...){
	
# Find the optimum I or J from a given range
# Enter a range, e.g. I=c(1,10), to optimise on that parameter
	
#start clock
tic <- proc.time()
	
# Determine which parameter needs to be optimised	
if(length(I)==2 & length(J)==1){print("Optimising on I ...")
	
	I.lower <- I[1]
	I.upper <- I[2]
	
	# Check they haven't entered 0 as a possible cluster size	
	if(I.upper==0 | I.lower==0){stop("Error: I must be greater than 0.")}
		
	power.storage <- matrix(c(rep(x=0, times=(I.upper - I.lower)+1)))
	rownames(power.storage) <- c(I.lower:I.upper)
	
	# Initially check upper and lower limit 	
	test.lower <- sim.power(I=I.lower,J=J,H=H,K=K,design=design,mu=mu,b.trt=b.trt,b.time=b.time,
                    sigma.y=sigma.y,sigma.e=sigma.e,rho=rho,sigma.a=sigma.a,
                    rho.ind=rho.ind,sigma.v=sigma.v,n.sims=n.sims,formula=formula,
                    family=family,natural.scale=natural.scale,sig.level=sig.level,n.cores=n.cores)
					
	power.lower <- test.lower$power
	# Store the lower power limit in the storeage list
	power.storage[I.lower - (I.lower-1)] <- power.lower
					
	if(power.lower >= target.power){print("Lower bound of range already above power threshold. Try a smaller lower bound.")
										return(I.lower)
									   break}			
					
	test.upper <- sim.power(I=I.upper,J=J,H=H,K=K,design=design,mu=mu,b.trt=b.trt,b.time=b.time,
                    sigma.y=sigma.y,sigma.e=sigma.e,rho=rho,sigma.a=sigma.a,
                    rho.ind=rho.ind,sigma.v=sigma.v,n.sims=n.sims,formula=formula,
                    family=family,natural.scale=natural.scale,sig.level=sig.level,n.cores=n.cores)
					
	power.upper <- test.upper$power
	# Store the upper power limit in the storeage list
	power.storage[I.upper - (I.lower-1)] <- power.upper
					
	if(power.upper <= target.power){print("Range does not include optimum. Try a larger upper bound.")
										return(I.upper)
										break}
										
										
	# Loop that successively narrows range by half based on power of the midpoint		
	while(TRUE){


		range <- c(I.lower:I.upper)			# Range of I to be checked
		midpoint <- ceiling(median(range))	# Middle value - must be an integer	

	if(length(range) != 2){	
			
	# If power has already been calculated (i.e. is stored in power.storage), use that
	if(power.storage[midpoint - (I[1]-1)]!=0){midpoint.power <- power.storage[midpoint+I[1]]
										  cat("used store value \n")
										  cat("midpoint is:", midpoint, "range is:", range, "power is:", midpoint.power, "\n")
										  cat("powers are", power.storage, "\n")}
	
	# Otherwise calculate the power using sim.power
	else if(power.storage[midpoint - (I[1]-1)]==0){	

		# If the range has more than 2 elements, use the midpoint of the range to discard one half of the range based on the power
		
		
			# Calculate power for SWT based on midpoint
			midpoint.power <- sim.power(I=midpoint,J=J,H=H,K=K,design=design,mu=mu,b.trt=b.trt,b.time=b.time,
                    sigma.y=sigma.y,sigma.e=sigma.e,rho=rho,sigma.a=sigma.a,
                    rho.ind=rho.ind,sigma.v=sigma.v,n.sims=n.sims,formula=formula,
                    family=family,natural.scale=natural.scale,sig.level=sig.level,n.cores=n.cores)$power
				
		power.storage[midpoint - (I[1]-1)] <- midpoint.power
		
		}			
			# If the power of the midpoint is too low, discard the lower half of the range
			if (midpoint.power < target.power) {I.upper <- I.upper
												I.lower <- midpoint}
			# If the power of the midpoint is high enough, discard the upper half of the range
			if (midpoint.power >= target.power) {I.upper <- midpoint
											     I.lower < I.lower}
		}	

		# If the range has only 2 elements, compare them directly
		if(length(range) == 2){
		
			# Calculate power for SWT based on the smaller value			
			if(power.storage[I.lower - (I[1]-1)]!=0){lower.power <- power.storage[I.lower - (I[1]-1)]}	# Utilise stored value if appropriate
			
			else{lower.power <- sim.power(I=I.lower,J=J,H=H,K=K,design=design,mu=mu,b.trt=b.trt,b.time=b.time,
                    sigma.y=sigma.y,sigma.e=sigma.e,rho=rho,sigma.a=sigma.a,
                    rho.ind=rho.ind,sigma.v=sigma.v,n.sims=n.sims,formula=formula,
                    family=family,natural.scale=natural.scale,sig.level=sig.level,n.cores=n.cores)$power
					}
			
			# Calculate power for SWT based on the larger value			
			if(power.storage[I.upper - (I[1]-1)]!=0){upper.power <- power.storage[I.upper - (I[1]-1)]}	# Utilise stored value if appropriate
			
			else{upper.power <- sim.power(I=I.upper,J=J,H=H,K=K,design=design,mu=mu,b.trt=b.trt,b.time=b.time,
                    sigma.y=sigma.y,sigma.e=sigma.e,rho=rho,sigma.a=sigma.a,
                    rho.ind=rho.ind,sigma.v=sigma.v,n.sims=n.sims,formula=formula,
                    family=family,natural.scale=natural.scale,sig.level=sig.level,n.cores=n.cores)$power
					}
					
			#end time counter
			toc <- proc.time(); time2run <- (toc-tic)[3]; names(time2run) <- "Time to run (secs)"
					
			# If either has power < target.power, discard 
			if(lower.power < target.power && upper.power < target.power){ return(0)}
			# If upper limit has enough power, but lower limit does not, return the upper limit
			if(lower.power < target.power && upper.power >= target.power){ return(list(Optimum_I=I.upper, power=upper.power, time2run=time2run))}
			# If both lmits have sufficient power, return the lower limit
			if(lower.power >= target.power && upper.power >= target.power){ return(list(Optimum_I=I.lower, power=lower.power, time2run=time2run))}
			# Other situations should not be possible
			else{return(0)
			# Optimum has been found, so break the loop and terminate the function
			break}
			}						
	}
}
		
# Determine which paramteter needs to be optimised	
if(length(J)==2 & length(I)==1){print("Optimising on J ...")}	

	J.lower <- J[1]
	J.upper <- J[2]
	
	# Check they haven't entered 0 as a possible number of time points	
	if(J.upper==0 | J.lower==0){stop("Error: J must be greater than 0.")}
		
		
	power.storage <- matrix(c(rep(x=0, times=(J.upper - J.lower)+1)))
	rownames(power.storage) <- c(J.lower:J.upper)
	
	# Initially check upper and lower limit 	
	test.lower <- sim.power(J=J.lower,I=I,H=H,K=K,design=design,mu=mu,b.trt=b.trt,b.time=b.time,
                    sigma.y=sigma.y,sigma.e=sigma.e,rho=rho,sigma.a=sigma.a,
                    rho.ind=rho.ind,sigma.v=sigma.v,n.sims=n.sims,formula=formula,
                    family=family,natural.scale=natural.scale,sig.level=sig.level,n.cores=n.cores)
					
	power.lower <- test.lower$power
	# Store the lower power limit in the storeage list
	power.storage[J.lower - (J.lower-1)] <- power.lower
					
	if(power.lower >= target.power){print("Lower bound of range already above power threshold. Try a smaller lower bound.")		
										return(J.lower)
										break}			
					
	test.upper <- sim.power(J=J.upper,I=I,H=H,K=K,design=design,mu=mu,b.trt=b.trt,b.time=b.time,
                    sigma.y=sigma.y,sigma.e=sigma.e,rho=rho,sigma.a=sigma.a,
                    rho.ind=rho.ind,sigma.v=sigma.v,n.sims=n.sims,formula=formula,
                    family=family,natural.scale=natural.scale,sig.level=sig.level,n.cores=n.cores)
					
	power.upper <- test.upper$power
	# Store the upper power limit in the storeage list
	power.storage[J.upper - (J.lower-1)] <- power.upper
					
	if(power.upper <= target.power){print("Range does not include optimum. Try a larger upper bound.")
										return(J.upper)
										break}
										
										
	# Loop that successively narrows range by half based on power of the midpoint		
	while(TRUE){


		range <- c(J.lower:J.upper)			# Range of I to be checked
		midpoint <- ceiling(median(range))	# Middle value - must be an integer	

	if(length(range) != 2){	
			
	# If power has already been calculated (i.e. is stored in power.storage), use that
	if(power.storage[midpoint - (J[1]-1)]!=0){midpoint.power <- power.storage[midpoint+J[1]]
										  cat("used store value \n")
										  cat("midpoint is:", midpoint, "range is:", range, "power is:", midpoint.power, "\n")
										  cat("powers are", power.storage, "\n")}
	
	# Otherwise calculate the power using sim.power
	else if(power.storage[midpoint - (J[1]-1)]==0){	

		# If the range has more than 2 elements, use the midpoint of the range to discard one half of the range based on the power
		
		
			# Calculate power for SWT based on midpoint
			midpoint.power <- sim.power(J=midpoint,I=I,H=H,K=K,design=design,mu=mu,b.trt=b.trt,b.time=b.time,
                    sigma.y=sigma.y,sigma.e=sigma.e,rho=rho,sigma.a=sigma.a,
                    rho.ind=rho.ind,sigma.v=sigma.v,n.sims=n.sims,formula=formula,
                    family=family,natural.scale=natural.scale,sig.level=sig.level,n.cores=n.cores)$power
				
		power.storage[midpoint - (J[1]-1)] <- midpoint.power
		
		
		}			
			# If the power of the midpoint is too low, discard the lower half of the range
			if (midpoint.power < target.power) {J.upper <- J.upper
												J.lower <- midpoint}
			# If the power of the midpoint is high enough, discard the upper half of the range
			if (midpoint.power >= target.power) {J.upper <- midpoint
											     J.lower < J.lower}
		}	

		# If the range has only 2 elements, compare them directly
		if(length(range) == 2){
		
			# Calculate power for SWT based on the smaller value			
			if(power.storage[J.lower - (J[1]-1)]!=0){lower.power <- power.storage[J.lower - (J[1]-1)]}  # Utilise stored value if appropriate
			
			else{lower.power <- sim.power(J=J.lower,I=I,H=H,K=K,design=design,mu=mu,b.trt=b.trt,b.time=b.time,
                    sigma.y=sigma.y,sigma.e=sigma.e,rho=rho,sigma.a=sigma.a,
                    rho.ind=rho.ind,sigma.v=sigma.v,n.sims=n.sims,formula=formula,
                    family=family,natural.scale=natural.scale,sig.level=sig.level,n.cores=n.cores)$power
					}
			
			# Calculate power for SWT based on the larger value			
			if(power.storage[J.upper - (J[1]-1)]!=0){upper.power <- power.storage[J.upper - (J[1]-1)]}	# Utilise stored value if appropriate
			
			else{upper.power <- sim.power(J=J.upper,I=I,H=H,K=K,design=design,mu=mu,b.trt=b.trt,b.time=b.time,
                    sigma.y=sigma.y,sigma.e=sigma.e,rho=rho,sigma.a=sigma.a,
                    rho.ind=rho.ind,sigma.v=sigma.v,n.sims=n.sims,formula=formula,
                    family=family,natural.scale=natural.scale,sig.level=sig.level,n.cores=n.cores)$power
					}
					
			#end timer
			toc <- proc.time(); time2run <- (toc-tic)[3]; names(time2run) <- "Time to run (secs)"			
					
			# If either has power < target.power, discard 
			if(lower.power < target.power && upper.power < target.power){ return(0)}
			# If upper limit has enough power, but lower limit does not, return the upper limit
			if(lower.power < target.power && upper.power >= target.power){ return(list(Optimum_J=J.upper, power=upper.power, time2run=time2run))}
			# If both lmits have sufficient power, return the lower limit
			if(lower.power >= target.power && upper.power >= target.power){ return(list(Optimum_J=J.lower, power=lower.power, time2run=time2run))}
			# Other situations should not be possible
			else{return(0)
			# Optimum has been found, so break the loop and terminate the function
			break}
			}						
	}
	
if(length(I) != 2 & length(J) != 2){stop("Error: exactly one of I or J must be a vector of length 2.")}	
	
		#stop clock
		proc.time() - tic
}
