# --------------------------------------------------------------------------
# power (or alpha) of 2 TOST via simulations
#
# Author D.L.
# --------------------------------------------------------------------------

power.2TOST.sim <- function(alpha=c(0.05,0.05), logscale=TRUE, theta0, 
                            theta1, theta2,  CV, n, rho, design="2x2", 
                            robust=FALSE, nsims, setseed=TRUE, details=FALSE)
{
  # Args:
  #   alpha: Vector of one-sided alpha levels for each procedure
  #   logscale: logical; if TRUE, log-transformed data is assumed
  #   theta0: Vector of 'true' assumed ratio of geometric means
  #   theta1: Vector of lower (bio-)equivalence limits
  #   theta2: Vector of upper (bio-)equivalence limits
  #   CV: Vector of coefficient of variations (use e.g. 0.3 for 30%)
  #       for the two metrics 
  #   n: For balanced allocation to treatment/sequence groups n is the total 
  #      sample size. For unbalanced allocation n is a vector with sample size
  #      per treatment/sequence group (e.g. n = c(8, 10) for an unbalanced
  #      2x2 design with total sample size 18)
  #   rho: Correlation between the two parameters under consideration
  #   design: Character string describing the study design
  #   robust: logical; if TRUE, use robust degrees of freedom (n - #sequences)
  #   setseed: pmvt() is based on randomized quasi Monte Carlo methods.
  #            Set seed value for (pseudo) random number generator?
  #   t1e: logical; if TRUE, Type I Error will be calculated
  #   details: logical; if TRUE, return P(Type I error) for all intersection
  #            null sets, max P(Type I error) otherwise
  if(length(alpha)==1) {
    alpha <- rep(alpha, 2)
    message("Scalar alpha used for both metrics.")
  }
  # is length(alpha=2) an real option
  # never heard that two parallel TOST's were done with different alpha
  if(length(alpha) != 2) stop("alpha must have two elements.")
  
  if (!missing(theta0) && length(theta0) != 2)
    stop("Two theta0 values must be given!")
  if (!missing(theta1) && length(theta1) != 2)
    stop("Two theta1 values must be given!")
  if (!missing(theta2) && length(theta2) != 2)
    stop("Two theta2 values must be given!")
  
  if (logscale) {
    if (missing(theta0))  theta0 <- c(0.95, 0.95)
    if (any(theta0 <= 0)) stop("theta0 must be > 0.")
    if (missing(theta1))  theta1 <- c(0.8, 0.8)
    if (missing(theta2))  theta2 <- 1/theta1
    if (any(theta1 <= 0) || any(theta1 > theta2))
      stop("theta1 and/or theta2 not correctly specified.")
  } else {
    if (missing(theta0)) theta0 <- c(0.05, 0.05)
    if (missing(theta1)) theta1 <- c(-0.2, -0.2)
    if (missing(theta2)) theta2 <- -theta1
    if (any(theta1 > theta2))
      stop("theta1 and/or theta2 not correctly specified.")
  }
  theta <- cbind(theta1, theta2)
    
  if(missing(CV))  stop("CVs must be given.")
  else {
    if(length(CV)==1) CV <- rep(CV,2)
    if(length(CV)!=2) stop("CV must have 2 elements.")
    if(any(CV<=0))   stop("CVs must be >0.")
  }
  
  if(missing(n)) stop("Number of subjects must be given.")

  # correlation
  if (missing(rho))
    stop("Correlation between the two endpoints must be given!")
  if (length(rho) != 1)
    stop("One rho must be given!")
  if (rho <= -1 || rho > 1)
    stop("Correlation must be > -1 and < 1.")

    # design characteristics
  # check if design is implemented
  d.no <- .design.no(design)
  if (is.na(d.no)) stop("Design ",design, " unknown!", call.=FALSE)
  
  # design characteristics
  ades <- .design.props(d.no)
  # degrees of freedom as expression
  dfe  <- .design.df(ades, robust=robust)
  # we use always bkni
  # step size = no of sequences
  steps <- ades$steps
  
  # handle n = ntotal if scalar else n's of the sequence groups
  if (length(n) == 1) {
    # total n given    
    # for unbalanced designs we divide the ns by ourself
    # to have only small imbalance (function nvec() from Helper_dp.R)
    if(is.finite(n)) n <- nvec(n=n, grps=ades$steps) else n <- rep(Inf, times=steps)
    if (n[1]!=n[length(n)]){
      message("Unbalanced design. n(i)=", paste(n, collapse="/"), " assumed.")
    } 
  } else {
    if (length(n) != ades$steps) {
      stop("Length of n vector must be ", ades$steps, "!")
    }
    if (any(n<1)) stop("All n(i) have to be >0.")
  }
  nc <- sum(1/n)
  n  <- sum(n)
  # Cfact*mse gives the variance of the mean difference
  Cfact <- ades$bkni * nc
  
  if(missing(nsims)){
    nsims <- 1E5
    if(any(theta0<=theta[, 1]) | any(theta0>=theta[,2])) nsims <- 1E6
  }
  
  # ---------------------------------------------------------------------------
  # now the calculations are starting
  # start timer
  ptm  <- proc.time()
  
  if (setseed) set.seed(1234567)
  
  sigma <- matrix(0, nrow=2, ncol=2)
  if(logscale) {
    ltheta0 <- log(theta0)
    ltheta  <- log(theta)
    diag(sigma) <- CV2mse(CV)
    sigma[1,2]  <- sigma[2,1] <- rho*sqrt(sigma[1,1]*sigma[2,2])
  } else {
    ltheta0 <- theta0
    ltheta  <- theta
    diag(sigma) <- CV^2 # ??? is this correct?
    sigma[1,2]  <- sigma[2,1] <- rho*sqrt(sigma[1,1]*sigma[2,2])
  }

  df    <- eval(dfe) 
  tvals <- qt(1-alpha, df)
  # cave sqrt() in case of rho<0
  s2m   <- sigma*Cfact
  
  # to avoid memory problems
  chunksize <- 1e7
  chunks    <- 1
  nsi       <- nsims
  if(nsims > chunksize) {
     chunks <- ceiling(nsims/chunksize)
     nsi    <- chunksize
  } 
  BE <- 0
  # over chunks of 1E7 if nsims > 1E7 to avoid memory problems
  for (iter in 1:chunks){
    # if chunks*1E7 >nsims correct nsi to given nsims
    if(iter==chunks & chunks>1) nsi <- nsims-(chunks-1)*nsi
    pes  <- matrix(0, nrow=nsi, ncol=2)
    mses <- matrix(0, nrow=nsi, ncol=2)
    # next: if directive with rho==0 was only for comparative purposes
    # no need to have it. multivariate NV and Wishart can handle this case
    # but retained by reason of an eventual speed gain
    if (rho==0){
      # simulate point est. via 2 independen normal distributions
      pes[ , 1] <- rnorm(n=nsi, mean=ltheta0[1], sd=sqrt(s2m[1,1])) # metric 1, f.i. AUC
      pes[ , 2] <- rnorm(n=nsi, mean=ltheta0[2], sd=sqrt(s2m[2,2])) # metric 2, f.i. Cmax
      # simulate mse via 2 independent chi-squared distributions
      mses[ , 1] <- sigma[1,1]*rchisq(n=nsi, df=df)/df
      mses[ , 2] <- sigma[2,2]*rchisq(n=nsi, df=df)/df
    } else {
      # pes multivariate normal with rho
      pes <- mvtnorm::rmvnorm(nsi, mean=ltheta0, sigma=s2m)
      # simulate covariance matrices via Wishart distribution
      covm      <- stats::rWishart(n=nsi, df=df, Sigma=sigma)/df 
      mses[, 1] <- covm[1, 1, ]
      mses[, 2] <- covm[2, 2, ]
      # take care of memory
      rm(covm)
    } 

    BE_test <- function(nu)
    {
      # nu = number of the metric
      hw <- tvals[nu]*sqrt(Cfact*mses[, nu])
      loCL <- pes[, nu] - hw
      upCL <- pes[, nu] + hw
      BE <- loCL >= ltheta[nu, 1] & upCL <= ltheta[nu, 2]
      BE
    }
    # make BE decision for both metrics
    BE_m1 <- BE_test(nu=1)
    BE_m2 <- BE_test(nu=2)
    
    # overall BE 
    BE <- BE + sum(BE_m1 & BE_m2)
  } # end over chunks
  
  if(details){
    cat(nsims, "simulations. Time consumed (secs)\n")
    print(proc.time()-ptm)
    cat("\n")
  }
  return(BE/nsims)
  
} #end function
