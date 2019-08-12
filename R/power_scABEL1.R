# sims based on EMA key statistics
# small CV correction of the ABE contribution (CV<=0.2)
power.scABEL1 <- function(alpha=0.05, theta1, theta2, theta0, CV, n,   
                          design=c("2x3x3", "2x2x4", "2x2x3"), regulator,
                          nsims=1E5, details=FALSE, setseed=TRUE)
{
  if (missing(CV)) stop("CV must be given!")
  if (missing(n))  stop("Number of subjects n must be given!")
  
  if (missing(theta0)) theta0 <- 0.90
  if (length(theta0)>1) {
    theta0 <- theta0[2]
    warning(paste0("theta0 has to be scalar. theta0 = ",
                   theta0, " used."), call. = FALSE)
  }
  if (missing(theta1) & missing(theta2)) theta1 <- 0.8
  if (missing(theta2)) theta2 <- 1/theta1
  if (missing(theta1)) theta1 <- 1/theta2
  
  ptm <- proc.time()
  
  # subject-by-formulation interaction can't play a role here
  # since the model doesn't allow such term
  CVwT <- CV[1]
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
  # paranoia
  if(is.null(pe_constr)) pe_constr <- TRUE
  
  # check design argument
  design <- match.arg(design)
  if (design=="2x3x3") {
    seqs <- 3
    bkni <- 1/6
    bk   <- 1.5
    dfe   <- parse(text="2*n-3", srcfile=NULL)
    dfRRe <- parse(text="n-2", srcfile=NULL)
    # sd^2 (variance) of the differences T-R from their components
    # According to McNally (sD=0):
    #sd2  <- s2wT + s2wR/2
    # but this gives too low power compared to the tables of the 2 Laszlos
    #sd2  <- (s2wT + s2wR)/2 # used in v1.1-00 - v1.1-02, wrong
    # simulations with s2D=0 show:
    Emse  <- (s2wT + 2.0*s2wR)/3
    # warning in case of CVwR != CVwT
    if (abs(s2wT-s2wR)>1e-4){
      warning(paste("Heteroscedasticity in the 2x3x3 design may lead", 
                    "to poor accuracy of power;\nconsider the function",
                    "power.scABEL.sdsims() instead."), call.=FALSE)
    }
  }
  if (design=="2x2x4") {
    seqs  <- 2
    bkni  <- 1/4
    bk    <- 1
    dfe   <- parse(text="3*n-4", srcfile=NULL)
    dfRRe <- parse(text="n-2", srcfile=NULL)
    Emse  <- (s2wT + s2wR)/2
  }
  if (design=="2x2x3") { # TRT|RTR design or TRR|RTT
    seqs  <- 2
    bkni  <- 3/8
    bk    <- 1.5
    dfe   <- parse(text="2*n-3", srcfile=NULL)
    dfRRe <- parse(text="n/2-1", srcfile=NULL) # balanced design
    Emse  <- (s2wT + s2wR)/2
  }
  if (length(n)==1){
    # for unbalanced designs we divide the ns by ourself
    # to have 'smallest' imbalance
    nv <- nvec(n=n, grps=seqs)
    if (nv[1]!=nv[length(nv)]){
      message("Unbalanced design. n(i)=", paste(nv, collapse="/"), " assumed.")
    }
    C2 <- sum(1/nv)*bkni
    n <- sum(nv)
  } else {
    # check length
    if (length(n)!=seqs) stop("n must be a vector of length=",seqs,"!")
    C2 <- sum(1/n)*bkni
    nv <- n
    n <- sum(n)
  }
  
  if (design=="2x2x3"){
    dfTT <- nv[1]-1
    dfRR <- nv[2]-1
    Emse <- (nv[1]*(2*s2wT+s2wR)/3+nv[2]*(s2wT+2*s2wR)/3)/n
  }
  
  # point est. in log domain
  mlog <- log(theta0)
  # sd of the sample mean T-R (point estimate)
  sdm  <- sqrt(Emse*C2)
  
  df   <- eval(dfe)
  if (design!="2x2x3"){
    dfRR <- eval(dfRRe)
    # next is purely empirical for 2x3x3
    dfTT <- dfRR
  } 
  
  if(setseed) set.seed(123456)
  p <- .pwr.ABEL.ANOVA(mlog, sdm, Ccon=C2, Emse, df, s2wR, dfRR, s2wT, dfTT, 
                       design, nsims, CVswitch, r_const, CVcap, pe_constr,
                       ln_lBEL=log(theta1),ln_uBEL=log(theta2), alpha=alpha)
  
  if (details) {
    ptm <- summary(proc.time()-ptm)
    message(nsims," sims. Time elapsed (sec): ", 
            formatC(ptm["elapsed"], digits=2), "\n")
    #print(ptm)
    # return the vector of all counts
    names(p) <- c("p(BE)", "p(BE-wABEL)", "p(BE-pe)", "p(BE-ABE)")
    if (!pe_constr) p <- p[-3] # without pe constraint
    p
  } else {
    # return only the 'power'
    as.numeric(p["BE"])
  }
}

# --- working horse for power calculation
# sims based on EMA ANOVA
.pwr.ABEL.ANOVA <- function(mlog, sdm, Ccon, Emse, df, s2wR, dfRR, s2wT, dfTT, 
                            design, nsims, CVswitch, r_const, CVcap, pe_constr, 
                            ln_lBEL, ln_uBEL, alpha)
{
  tcrit    <- qt(1-alpha,df)
  s2switch <- log(1.0 + CVswitch^2)
  s2cap    <- log(1.0 + CVcap^2)
  capABEL  <- sqrt(s2cap)*r_const
  
  # result vector
  counts        <- rep.int(0, times=4)
  names(counts) <- c("BE", "BEwl", "BEpe", "BEabe")
  # to avoid memory problems
  chunks <- 1
  nsi    <- nsims
  if (nsims>1E7) {
    chunks <- ceiling(nsims/1E7)
    nsi    <- 1E7
  } 
  for (iter in 1:chunks){
    # if chunks*1E7 >nsims correct nsi to given nsims
    if(iter==chunks & chunks>1) nsi <- nsims-(chunks-1)*nsi
    # simulate sample mean via its normal distribution
    means  <- rnorm(nsi, mean=mlog, sd=sdm)
    # simulate sample value s2wT/s2wR via chi-square distri
    # shifting this after the sample mse sims will give
    # changes in the order of max. 4E-4 compared to V1.1-02!
    s2wRs <- s2wR*rchisq(nsi, dfRR)/dfRR
    s2wTs <- s2wT*rchisq(nsi, dfTT)/dfTT
    # now simulate sample mse
    if (design=="2x3x3"){
      mses_ABE  <- Emse*rchisq(nsi, df)/df
      # simulate sample mse not for this design
      # s2wT is empirical because dfTT is not defined and  
      # artificially set to equal dfRR
      mses  <- (s2wTs + 2*s2wRs)/3
    } 
    if (design=="2x2x4"){
      # simulate sample mse 
      mses_ABE  <- Emse*rchisq(nsi, df)/df
      #--- 'mean' of both attempts V1.1-02/V1.1-03 in case of 2x2x4
      mses  <- (mses_ABE + (s2wTs + s2wRs)/2)/2
    }
    if (design=="2x2x3"){
      # simulate sample mse 
      mses_ABE  <- Emse*rchisq(nsi, df)/df
      # --- 'mean' of mses and mses calculated as sum from components
      mses  <- 0.5*mses_ABE + 
        0.5*((dfTT+1)*(2*s2wTs + s2wRs)/3
             +(dfRR+1)*(s2wTs + 2*s2wRs)/3)/(dfTT+dfRR+2)
    }
    #--- EMA widened limits in log-domain
    uABEL   <- +sqrt(s2wRs)*r_const
    lABEL   <- -uABEL
    # use ABE limits if CV < CVswitch
    # lABEL[s2wRs<=s2switch] <- ln_lBEL
    # uABEL[s2wRs<=s2switch] <- ln_uBEL
    # cap limits if CVcap not Inf
    if (is.finite(CVcap)){
      lABEL[s2wRs>s2cap] <- -capABEL
      uABEL[s2wRs>s2cap] <- +capABEL
    }
    #--- 90% CIs for T-R
    hw  <- tcrit*sqrt(Ccon*mses)
    lCL <- means - hw 
    uCL <- means + hw
    #--- 90% CI in 'widened' limits?
    BE   <- (lABEL<=lCL) & (uCL<=uABEL)
    #--- conventional ABE for low CV (0.2 is purely empiric)
    if (s2wR<=CV2mse(0.2)){
      hw  <- tcrit*sqrt(Ccon*mses_ABE)
      lCL <- means - hw 
      uCL <- means + hw
    }
    BEABE <- (ln_lBEL<=lCL) & (uCL<=ln_uBEL)
    # if CV < CV switch use ABE, else scABEL
    BE    <- ifelse(s2wRs>s2switch, BE, BEABE)
    
    #--- point est. constraint true?
    BEpe <- (means>=ln_lBEL) & (means<=ln_uBEL)
    
    rm(hw)
    
    counts["BEabe"] <- counts["BEabe"] + sum(BEABE)
    counts["BEpe"]  <- counts["BEpe"]  + sum(BEpe)
    counts["BEwl"]  <- counts["BEwl"]  + sum(BE)
    if (pe_constr) {
      counts["BE"] <- counts["BE"] + sum(BE & BEpe)
    } else {
      counts["BE"] <- counts["BE"] + sum(BE)
    }
  } # end over chunks
  
  # return the counts
  counts/nsims
}
