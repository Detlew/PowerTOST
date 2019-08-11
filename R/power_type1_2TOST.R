# --------------------------------------------------------------------------
# power (or alpha aka type 1 error aka TIE) of 2 TOST via simulations
#
# Author(s) D. Labes, B. Lang
# --------------------------------------------------------------------------
power.2TOST <- function(alpha = c(0.05, 0.05), logscale = TRUE, theta0, theta1,
                        theta2, CV, n, rho, design = "2x2", robust = FALSE, 
                        nsims, setseed = TRUE, details = FALSE) {
  prob.2TOST(alpha = alpha, logscale = logscale, theta0 = theta0, 
             theta1 = theta1, theta2 = theta2, CV = CV, n = n, rho = rho, 
             design = design, robust = robust, nsims = nsims, setseed = setseed,
             details = details)
}
# type1error.2TOST without theta0
type1error.2TOST <- function(alpha = c(0.05, 0.05), logscale = TRUE, 
                             theta1, theta2, CV, n, rho, design = "2x2", 
                             robust = FALSE, nsims, setseed = TRUE, 
                             details = FALSE) {
  .Defunct(msg=paste0("This function is no longer available since it suffers\n",
                      "from insufficient precision to obtain the type 1 error (TIE)\n",
                      "via simulations."))
  
  prob.2TOST(alpha = alpha, logscale = logscale,  
             theta1 = theta1, theta2 = theta2, CV = CV, n = n, rho = rho, 
             design = design, robust = robust, nsims = nsims, setseed = setseed,
             t1e = TRUE, details = details)
}

prob.2TOST <- function(alpha=c(0.05,0.05), logscale=TRUE, theta0, theta1,
                       theta2, CV, n, rho, design="2x2", robust=FALSE, nsims, 
                       setseed=TRUE, t1e = FALSE, details=FALSE)
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
  #   nsims: number of simulations
  #   setseed: Set seed value for (pseudo) random number generator?
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

  if(missing(CV))  stop("CVs must be given.")
  else {
#    if(length(CV)==1) CV <- rep(CV,2)
    if(length(CV)!=2) stop("CV must have 2 elements.")
    if(any(CV<=0))   stop("CVs must be >0.")
  }
  
  if(missing(n)) stop("Number of subjects must be given.")

  # correlation
  if (missing(rho))
    stop("Correlation between the two endpoints must be given!")
  if (length(rho) != 1){
    warning("Only one rho must be given. First entry used.")
    rho <- rho[1]
  }
  if (rho < -1 || rho > 1)
    stop("Correlation must be within {-1, +1}.")
  # allow -1, +1 within machine epsilon [HS]
  if (rho == -1) rho <- -1 + .Machine $double.eps
  if (rho == +1) rho <- +1 - .Machine $double.eps

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
  
  if (missing(nsims)) {
    nsims <- 1E5
    if (t1e) nsims <- 1E6
  }
  
  # 'true' variance-covariance matrix
  sigma <- matrix(0, nrow=2, ncol=2)
  if(logscale) {
    ltheta0 <- log(theta0)
    ltheta1 <- log(theta1)
    ltheta2 <- log(theta2)
    diag(sigma) <- CV2mse(CV)
    sigma[1,2]  <- sigma[2,1] <- rho*sqrt(sigma[1,1]*sigma[2,2])
  } else {
    ltheta0 <- theta0
    ltheta1  <- theta1
    ltheta2  <- theta2
    diag(sigma) <- CV^2 # ??? is this correct?
    sigma[1,2]  <- sigma[2,1] <- rho*sqrt(sigma[1,1]*sigma[2,2])
  }
  
  # now the calculations are starting
  # start timer
  ptm  <- proc.time()
  df  <- eval(dfe)
  if (t1e) {
    # Calculate Type I Error
    
    twodim = FALSE # via 2 dimensional optim() or via one dimensional optimize()
    
    if (setseed) set.seed(1234567)
    nullsets <- vector("list", 8)
    lim <- 100  # 100 instead of Inf; suffices and avoids potential
    # optimization problems when using Inf
    #lim <- Inf
    H_A01 <- c(-lim, ltheta1[1])
    H_A02 <- c(ltheta2[1], lim)
    H_C01 <- c(-lim, ltheta1[2])
    H_C02 <- c(ltheta2[2], lim)
    H_Aa <- c(ltheta1[1], ltheta2[1])
    H_Ca <- c(ltheta1[2], ltheta2[2])
    nullsets[[1]] <- rbind(H_A01, H_Ca)
    nullsets[[2]] <- rbind(H_A02, H_Ca)
    nullsets[[3]] <- rbind(H_Aa, H_C01)
    nullsets[[4]] <- rbind(H_Aa, H_C02)
    nullsets[[5]] <- rbind(H_A01, H_C01)
    nullsets[[6]] <- rbind(H_A01, H_C02)
    nullsets[[7]] <- rbind(H_A02, H_C01)
    nullsets[[8]] <- rbind(H_A02, H_C02)
    onedim <- function(H) 
      {
      # function to optimize
      fun <- function(x, ix, ltheta0,
                      alpha = alpha, df = df, Cfact = Cfact, 
                      ltheta1 = ltheta1, ltheta2 = ltheta2, 
                      sigma = sigma, rho = rho, nsims = nsims)
        {
        ltheta0[ix] <- x
        pwr <- .prob.2TOST(ltheta0=ltheta0, alpha = alpha, df = df,
        Cfact = Cfact, ltheta1 = ltheta1, ltheta2 = ltheta2, 
        sigma = sigma, rho = rho, nsims = nsims)
        pwr
      } 
      lower <- H[, 1]
      upper <- H[, 2]
      ltheta0 <- c(.getStart(H[1, ], lim), .getStart(H[2, ], lim))
      if(lower[1] == -lim | upper[1] == lim){
        ix <- 2
      } else ix <- 1
      res <- optimize(f=fun, interval=H[ix, ], ix=ix, ltheta0=ltheta0,
                      alpha = alpha, df = df, Cfact = Cfact, 
                      ltheta1 = ltheta1, ltheta2 = ltheta2, 
                      sigma = sigma, rho = rho, nsims = nsims,
                      maximum=TRUE, tol=1e-5)
      #browser()
      ltheta0[ix] <- res$maximum
      argmax <- if (logscale) exp(ltheta0) else ltheta0
      c(res$objective, argmax)
    } # function onedim
    
    size.H <- function(H) {  # size of test (supremum over H)
      # Determine starting values (theta0[1], theta0[2])
      starting <- c(.getStart(H[1, ], lim), .getStart(H[2, ], lim))
      # Perform maximization over intersection nullset H via optim()
      res <- optim(par = starting, fn = .prob.2TOST, alpha = alpha, df = df,
                   Cfact = Cfact, ltheta1 = ltheta1, ltheta2 = ltheta2, 
                   sigma = sigma, rho = rho, nsims = nsims, 
                   method = "L-BFGS-B", lower = H[, 1], upper = H[, 2], 
                   control = list(fnscale = -1, ndeps=rep(1e-4,2)))
      if (!(res$convergence %in% c(0, 1))) {
        nsname <- paste(row.names(H), collapse=" n ")
        warning("Result of maximization over nullset '", nsname, 
                "' may not be reliable.\n optim() message: ", 
                res$message, call. = FALSE) 
      }
      # Select maximum
      res.argmax <- if (logscale) exp(res$par) else res$par
      res.max <- res$value
      c(res.max, res.argmax)
    }  # End of size.H
    #browser()
    # Combine results
    # original:
    #probs <- as.data.frame(t(vapply(nullsets, size.H, numeric(3))))
    nullsets14 <- nullsets[1:4]
    if( twodim) probs <- as.data.frame(t(vapply(nullsets14, size.H, numeric(3))))
      else probs <- as.data.frame(t(vapply(nullsets14, onedim, numeric(3))))
    # recalculate with setseed to avoid discepancies to call of power.2TOST()
    if(setseed) {
      for(i in 1:4) {
        lth0 <- as.numeric(probs[i, 2:3])
        if(logscale) lth0 <- log(lth0)
        set.seed(123456)
        pwr <- .prob.2TOST(ltheta0=lth0, alpha = alpha, df = df,
                           Cfact = Cfact, ltheta1 = ltheta1, ltheta2 = ltheta2, 
                           sigma = sigma, rho = rho, nsims = nsims)
        probs$V1[i] <- pwr
      }
    }
    # nullsets [5:8] via direct call of .prob.2TOST
    nullsets58 <- nullsets[5:8]
    for(i in 1:4) {
      H <- nullsets58[[i]]
      starting <- c(.getStart(H[1, ], lim), .getStart(H[2, ], lim))
      if(setseed) set.seed(123456)
      pwr <- .prob.2TOST(ltheta0=starting, alpha = alpha, df = df,
                         Cfact = Cfact, ltheta1 = ltheta1, ltheta2 = ltheta2, 
                         sigma = sigma, rho = rho, nsims = nsims)
      if(logscale) starting <- exp(starting)
      prow <- c(pwr, starting)
      probs <- rbind(probs, prow)
    }
    nullsets.label <- c("H_A01 n H_Ca", "H_A02 n H_Ca", "H_Aa n H_C01",
                        "H_Aa n H_C02", "H_A01 n H_C01", "H_A01 n H_C02",
                        "H_A02 n H_C01", "H_A02 n H_C02")
    probs <- cbind(nullsets.label, probs)
    colnames(probs) <- c("Intersection null", "P(Type I Error)", "theta0 #1", 
                         "theta0 #2")
    prob <- if (details) probs else max(probs[["P(Type I Error)"]])
  } else {
    # Calculate Power
    if (setseed) set.seed(123456)
    prob <- .prob.2TOST(ltheta0 = ltheta0, alpha = alpha, df = df, Cfact = Cfact, 
                        ltheta1 = ltheta1, ltheta2 = ltheta2, sigma = sigma, 
                        rho = rho, nsims = nsims)
  }
  if (details) {
    cat(nsims, "simulations. Time consumed (secs)\n")
    print(proc.time()-ptm)
    cat("\n")
  }
  prob
}

# ---------------------------------------------------------------------------
# working horse. to be used also in sampleN.2TOST()
# shifting setseed into this function doubles the run-time in
# type1error.2TOST()
.prob.2TOST <- function(ltheta0, alpha, df, Cfact, ltheta1, ltheta2, sigma, 
                        rho, nsims)
{
  tvals <- qt(1-alpha, df)
  # cave sqrt() in case of rho<0
  s2m   <- sigma*Cfact
  
  # over chunks of 1E7 if nsims > 1E7 to avoid memory problems
  chunksize <- 1e7
  chunks    <- 1
  nsi       <- nsims
  if(nsims > chunksize) {
     chunks <- ceiling(nsims/chunksize)
     nsi    <- chunksize
  } 
  BE <- 0
  for (iter in 1:chunks){
    # if chunks*1E7 > nsims then correct nsi to given rest of nsims
    if(iter==chunks & chunks>1) nsi <- nsims-(chunks-1)*nsi
    # reserve memory
    pes  <- matrix(logical(0), nrow=nsi, ncol=2)
    mses <- matrix(logical(0), nrow=nsi, ncol=2)
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
      #pes <- mvnfast::rmvn(nsi, mu=ltheta0, sigma=s2m)
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
      BE <- loCL >= ltheta1[nu] & upCL <= ltheta2[nu]
      BE
    }
    # make BE decision for both metrics
    BE_m1 <- BE_test(nu=1)
    BE_m2 <- BE_test(nu=2)
    
    # overall BE 
    BE <- BE + sum(BE_m1 & BE_m2)
  } # end over chunks
  
  return(BE/nsims)
  
} #end function

.getStart <- function(H, lim) {
  if (H[1] <= -lim) {
    return(H[2])
  } else if (H[2] >= lim) {
    return(H[1])
  } else {
    return((H[1] + H[2]) / 2)
  }
}
