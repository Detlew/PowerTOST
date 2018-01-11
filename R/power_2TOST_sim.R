# --------------------------------------------------------------------------
# power (or alpha) of 2 TOST via simulations
#
# Author D.L.
# --------------------------------------------------------------------------

power.2TOST.sim <- function(alpha=c(0.05,0.05), logscale=TRUE, CV, rho=0, n, 
                            theta0, theta, design="2x2", robust=FALSE,
                            nsims, setseed=TRUE, details=FALSE)
{
  if(missing(CV))  stop("CVs must be given.")
  else {
    if(length(CV)==1) CV <- rep(CV,2)
    if(length(CV)!=2) stop("CV must have 2 elements.")
    if(any(CV<=0))   stop("CVs must be >0.")
  }
  if(missing(n)) stop("Number of subjects must be given.")

  if(length(alpha)==1) {
    alpha <- rep(alpha, 2)
    message("Scalar alpha used for both metrics.")
  }
  if(length(alpha) != 2) stop("alpha must have two elements.")

  # 'true' GMR or difference  
  if(missing(theta0)){
    if(logscale){
      theta0 <- rep(0.95, 2)
    } else {
      theta0 <- rep(0.05, 2)
    }
  } else {
    if(length(theta0)==1) theta0 <- rep(theta0, 2)
    if(length(theta0)!=2) stop("theta0 must have 2 elements.")
    # more to check???
    
  }
  
  # BE acceptance limits, row 1 for 1. metric, row 2 for second
  if( missing(theta)) {
    if(logscale){
      theta <- matrix(c(0.8, 0.8, 1.25, 1.25), nrow=2)
    } else {
      theta <- matrix(c(-0.2, -0.2, 0.2, 0.2), nrow=2)
    }
    # more to check?
  }

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
  #Cfact <- sqrt(ades$bkni * nc)
  Cfact <- ades$bkni * nc
  
  #if(rho!=0) warning("rho != 0 is only experimental.", call. = FALSE)
  stopifnot(length(rho)==1, rho >= -1, rho <= 1)
  
  if(missing(nsims)){
    nsims <- 1E5
    if(any(theta0<=theta[, 1]) | any(theta0>=theta[,2])) nsims <- 1E6
  }
  
  # start timer
  ptm  <- proc.time()
  
  if (setseed) set.seed(1234567)
  
  sigma <- matrix(0,2,2)
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
  
  pes  <- matrix(0, nrow=nsims, ncol=2)
  mses <- matrix(0, nrow=nsims, ncol=2)
  if (rho==0){
    # simulate point est. via normal distribution
    pes[ , 1] <- rnorm(n=nsims, mean=ltheta0[1], sd=sqrt(s2m[1,1])) # metric 1, f.i. AUC
    pes[ , 2] <- rnorm(n=nsims, mean=ltheta0[2], sd=sqrt(s2m[2,2])) # metric 2, f.i. Cmax
    # simulate mse via chi-squared distribution
    mses[ , 1] <- sigma[1,1]*rchisq(n=nsims, df=df)/df
    mses[ , 2] <- sigma[2,2]*rchisq(n=nsims, df=df)/df
  } else {
    # multivariate normal with rho
    # TODO: check this
    pes <- mvtnorm::rmvnorm(nsims, mean=ltheta0, sigma=s2m)
    # simulate covariance matrices via Wishart distribution
    #or do we have here df=n-1?
    covm      <- rWishart(n=nsims, df=df, Sigma=sigma)/df 
    mses[, 1] <- covm[1, 1, ]
    mses[, 2] <- covm[2, 2, ]
    # take care of memory
    rm(covm)
  }
  #browser()
  
  BE <- function(nu)
  {
    # nu = number of the metric
    hw <- tvals[nu]*sqrt(Cfact*mses[, nu])
    loCL <- pes[, nu] - hw
    upCL <- pes[, nu] + hw
    BE <- loCL >= ltheta[nu, 1] & upCL <= ltheta[nu, 2]
    BE
  }
  # make BE decision for both metrics
  BE_m1 <- BE(nu=1)
  BE_m2 <- BE(nu=2)
  
  # overall BE 
  BE <- BE_m1 & BE_m2

  if(details){
    cat(nsims, "simulations. Time consumed (secs)\n")
    print(proc.time()-ptm)
    cat("\n")
  }
  return(sum(BE)/nsims)
  
} #end function
