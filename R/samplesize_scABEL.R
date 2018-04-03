#---------------------------------------------------------------------------
# unified function
# chooses the sample size function according to regulator$est_method
# former is now sampleN.scABEL1
#
# author dlabes
#---------------------------------------------------------------------------
sampleN.scABEL <- function(alpha=0.05, targetpower=0.8, theta0, theta1, 
                           theta2, CV, design=c("2x3x3", "2x2x4", "2x2x3"), 
                           regulator, nsims=1E5, nstart, imax=100, print=TRUE, 
                           details=TRUE, setseed=TRUE)
{
  # design must be checked outside
  desi <- match.arg(design)
  # check regulator
  if (missing(regulator)) regulator <- "EMA"
  reg  <- reg_check(regulator)
  ssfun <- "sampleN.scABEL1"
  if (reg$est_method=="ISC") ssfun <- "sampleN.scABEL2"
  # print(sys.call())
  # next doesn't function if arguments are missing
  # r <- do.call(ssfun,
  #              list(alpha, targetpower, theta0, theta1, theta2, CV, 
  #                   design=desi, regulator=reg, nsims, nstart, imax,
  #                   print, details, setseed))
  if (reg$est_method!="ISC"){
    r <- sampleN.scABEL1(alpha, targetpower, theta0, theta1, theta2, CV, 
                         design=desi, regulator=reg, nsims, nstart, imax,
                         print, details, setseed)
  } else {
    r <- sampleN.scABEL2(alpha, targetpower, theta0, theta1, theta2, CV, 
                         design=desi, regulator=reg, nsims, nstart, imax,
                         print, details, setseed)
  } 
  
}  

#-----------------------------------------------------------------------------
# Sample size for partial and full replicate design and scaled ABE 
# via simulated (empirical) power
# estimation method ANOVA and BE decision via ABEL (Average bioequivalence with 
# expanding limits)
# 
# Author: dlabes
#-----------------------------------------------------------------------------

# helper function: sample size for pe in a range?
# definition is more or less empirical (i.e. not understood by me)
.sampleN0.2 <- function(targetpower, ltheta2, diffm, se, bk, steps)
{
  n <- qnorm(targetpower)^2*se^2*bk/(abs(diffm)-ltheta2)^2
  #make an even multiple of step (=2 in case of 2x2 cross-over)
  n <- steps*trunc(n/steps)
  n
}  

sampleN.scABEL1 <- function(alpha=0.05, targetpower=0.8, theta0, theta1, 
                            theta2, CV, design=c("2x3x3", "2x2x4", "2x2x3"), 
                            regulator, nsims=1E5, nstart, imax=100, print=TRUE, 
                            details=TRUE, setseed=TRUE)
{
  if (missing(theta1) & missing(theta2)) theta1 <- 0.8
  if (missing(theta2)) theta2 <- 1/theta1
  if (missing(theta1)) theta1 <- 1/theta2
  # the two Laszlo's recommend theta0=0.9 for HVD's
  if (missing(theta0)) theta0 <- 0.9
  if ( (theta0<=theta1) | (theta0>=theta2) ) {
    stop("True ratio ",theta0," not between margins ",theta1," / ",
         theta2,"!", call.=FALSE)
  }
  if (missing(CV)) stop("CV must be given!", call.=FALSE)
  
  # subject-by-formulation interaction can't play a role here I think
  # since the EMA model doesn't allow such term
  CVwT <- CV[1]
  # should we allow different variabilities in the EMA method at all?
  if (length(CV)==2) CVwR <- CV[2] else CVwR <- CVwT
  s2wT <- log(1.0 + CVwT^2)
  s2wR <- log(1.0 + CVwR^2)
  
  if(missing(regulator)) regulator <- "EMA"
  # check regulator and get 
  # constants acc. to regulatory bodies (function in scABEL.R)
  reg <- reg_check(regulator)
  CVcap     <- reg$CVcap
  CVswitch  <- reg$CVswitch
  r_const   <- reg$r_const
  pe_constr <- reg$pe_constr
  
  # check design
  design <- match.arg(design)
  # we are treating only balanced designs
  # thus we use here bk - design constant for ntotal
  # expressions for the df's
  if (design=="2x3x3") {
    desi <- "2x3x3 (partial replicate)"
    bk <- 1.5; seqs <- 3
    dfe   <- parse(text="2*n-3", srcfile=NULL)
    dfRRe <- parse(text="n-2", srcfile=NULL)
    #sd2  <- (s2wT + s2wR)/2 # used in v1.1-00 - v1.1-02, wrong
    # simulations with s2D=0 show:
    Emse  <- (s2wT + 2.0*s2wR)/3
    cvec  <- c(1, 2) # for sim of mses from s2wT and s2wR
  }
  if (design=="2x2x4") {
    desi <- "2x2x4 (full replicate)"
    bk <- 1.0; seqs <- 2
    # only EMA settings
    dfe   <- parse(text="3*n-4", srcfile=NULL)
    dfRRe <- parse(text="n-2", srcfile=NULL)
    # sd^2 (variance) of the differences T-R from their components
    Emse  <- (s2wT + s2wR)/2
    cvec  <- c(1, 1)
  }
  if (design=="2x2x3") {
    desi <- "2x2x3 (TRT|RTR)"
    bk <- 1.5; seqs <- 2
    # only EMA settings
    dfe   <- parse(text="2*n-3", srcfile=NULL)
    dfRRe <- parse(text="n/2-1", srcfile=NULL)
    # sd^2 (variance) of the differences T-R from their components
    Emse  <- (s2wT + s2wR)/2 # for balanced designs we use here
    cvec  <- c(1, 1) # dummy
  }
  
  mlog <- log(theta0)
  
  if (print){
    cat("\n+++++++++++ scaled (widened) ABEL +++++++++++\n")
    cat("            Sample size estimation\n")
    cat("   (simulation based on ANOVA evaluation)\n")
    cat("---------------------------------------------\n")
    cat("Study design: ",desi,"\n")
    cat("log-transformed data (multiplicative model)\n")
    cat(nsims,"studies for each step simulated.\n\n")
    cat("alpha  = ", alpha,", target power = ", targetpower,"\n", sep="")
    cat("CVw(T) = ", CVwT,"; CVw(R) = ", CVwR,"\n", sep="")
    cat("True ratio = ",theta0,"\n", sep="")
    cat("ABE limits / PE constraint =", theta1,"...", theta2,"\n")
    if (details | reg$name=="USER") { 
      print(reg)
      cat("\n")
    } else {
      cat("Regulatory settings:", reg$name,"\n")
    }
  }
  
  # -----------------------------------------------------------------
  # nstart? from sampleN0 with widened limits
  # does'nt fit really good if theta0>=1.2! ways out? see sampleN0.2
  ltheta1 <- -sqrt(s2wR)*r_const
  ltheta2 <- -ltheta1
  if (CVwR <= CVswitch){
    ltheta1 <- log(theta1)
    ltheta2 <- log(theta2)
  }
  if (CVwR > CVcap){
    ltheta1 <- -sqrt(log(1.0 + CVcap^2))*r_const
    ltheta2 <- -ltheta1
  }
  if (missing(nstart)){
    # start from ABE start with widened limits
    n01 <- .sampleN0_2(alpha=alpha, targetpower, ltheta1, ltheta2, diffm=mlog, 
                     se=sqrt(Emse), steps=seqs, bk=bk)
    # empirical correction in the vicinity of CV=0.3 for ratios 
    # outside 0.86 ... 1/0.86
    # does this fit also for CVswitch in case of ANVISA = 0.4?
    if(Emse < CV2mse(CVswitch+0.005) & Emse > CV2mse(CVswitch-0.005) 
        & abs(mlog)>log(1/0.865)) {
      if (reg$name=="EMA")     n01 <- 0.9*n01
      if (reg$name=="FDA")     n01 <- 0.65*n01
      if (reg$name=="ANVISA")  n01 <- 0.6*n01
      n01 <- seqs*trunc(n01/seqs)
    }
    
    # start from PE constraint sample size
    n02 <- .sampleN0.2(targetpower, ltheta2=log(theta2), diffm=mlog, 
                       se=sqrt(Emse), steps=seqs, bk=bk)
    # debug print
    # cat(n01,n02,"\n")
    n <- max(c(n01,n02))
  } else n <- seqs*round(nstart/seqs)           
  nmin <- 6
  if (n<nmin) n <- nmin
  # we are simulating for balanced designs
  C2 <- bk/n
  # sd of the sample mean T-R (point estimator)
  sdm  <- sqrt(Emse*C2)
  df   <- eval(dfe)
  dfRR <- eval(dfRRe)
  dfTT <- dfRR

  if(setseed) set.seed(123456)
  p <- .pwr.ABEL.ANOVA(mlog, sdm, C2, Emse, cvec, df, s2wR, dfRR, s2wT, dfTT,
                       design, nsims, CVswitch, r_const, CVcap, pe_constr,
                       ln_lBEL=log(theta1),ln_uBEL=log(theta2), alpha=alpha)
  pwr <- as.numeric(p["BE"]);
  pd <- max(4,round(log10(nsims),0)-1)  # digits for power
  if (details) {
    cat("\nSample size search\n")
    cat(" n     power\n")
    # do not print first too high
    # this is for cases with only one step-down and than step up
    if (pwr<=targetpower) cat( n," ", formatC(pwr, digits=pd, format="f"),"\n")
  }
  iter <- 0
  # iter>100 is emergency brake
  # --- loop until power <= target power, step-down
  down <- FALSE
  up   <- FALSE
  while (pwr>targetpower) {
    down <- TRUE
    if (n<=nmin) { 
      if (details & iter==0) cat( n," ", formatC(pwr, digits=pd, format="f"),"\n")
      break
    }
    n  <- n-seqs     # step down if start power is to high
    iter <- iter + 1
    C2 <- bk/n
    # sd of the sample mean T-R (point estimator)
    sdm  <- sqrt(Emse*C2)
    df   <- eval(dfe)
    dfRR <- eval(dfRRe)
    dfTT <- dfRR
    if(setseed) set.seed(123456)
    p <- .pwr.ABEL.ANOVA(mlog, sdm, C2, Emse, cvec, df, s2wR, dfRR, s2wT, dfTT,
                         design, nsims, CVswitch, r_const, CVcap, pe_constr,
                         ln_lBEL=log(theta1),ln_uBEL=log(theta2), alpha=alpha)
    pwr <- as.numeric(p["BE"]);
    
    if (details) cat( n," ", formatC(pwr, digits=pd, format="f"),"\n")
    if (iter>imax) break  
    # loop results in n with power too low
    # must step up again one step. is done in the next loop
  }
  while (pwr<targetpower) {
    up   <- TRUE; down <- FALSE
    n    <- n+seqs   # step-up
    iter <- iter+1
    C2   <- bk/n
    sdm  <- sqrt(Emse*C2)
    df   <- eval(dfe)
    dfRR <- eval(dfRRe)
    dfTT <- dfRR
    if(setseed) set.seed(123456)
    p <- .pwr.ABEL.ANOVA(mlog, sdm, C2, Emse, cvec, df, s2wR, dfRR, s2wT, dfTT,
                         design, nsims, CVswitch, r_const, CVcap, pe_constr,
                         ln_lBEL=log(theta1),ln_uBEL=log(theta2), alpha=alpha)
    pwr <- as.numeric(p["BE"]);
    
    if (details) cat( n," ", formatC(pwr, digits=pd, format="f"),"\n")
    if (iter>imax) break 
  }

  nlast <- n
  if (up & pwr<targetpower) {
    n <- NA
    if (details) cat("Sample size search failed!\n")
  }
  if (down & pwr>targetpower & n != nmin) {
    n <- NA
    if (details) cat("Sample size search failed!\n")
  }
  
  if (print && !details) {
    cat("\nSample size\n")
    cat(" n     power\n")
    cat( n," ", formatC(pwr, digits=pd, format="f"),"\n")
    if (is.na(n)) cat("Sample size search failed!\n")
  }
  if (print) cat("\n")
  
  #return results as data.frame
  res <- data.frame(design=design, alpha=alpha, CVwT=CVwT, CVwR=CVwR,
                    theta0=theta0, theta1=theta1, theta2=theta2, n=n, power=pwr, 
                    targetpower=targetpower,nlast=nlast)
  names(res) <-c("Design", "alpha", "CVwT", "CVwR", "theta0", "theta1", "theta2",
                 "Sample size", "Achieved power", "Target power", "nlast")

  #cat("iter=",iter,"\n")
  
  if (print | details) return(invisible(res)) else return(res)
  
} # end function
