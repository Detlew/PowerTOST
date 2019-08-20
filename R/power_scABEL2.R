#---------------------------------------------------------------------------
# Simulate partial and full replicate design and scaled ABE power
# using 'widened' limits (EMA)
# But the estimation method via intra-subject contrasts (FDA, HC)
#
# Author: dlabes
#---------------------------------------------------------------------------

# degrees of freedom for the TR/RR  analysis: 
# Using the intrasubject contrasts T-R and R-R and analyze them  
# by sequence groups the df's = n-seq (robust df's).
# 2x3x3  dfRR = n-3
# 2x2x4  dfRR = n-2
# 2x2x3  dfRR = n/2 - 2

power.scABEL2 <- function(alpha=0.05, theta1, theta2, theta0, CV, n,   
                          design=c("2x3x3", "2x2x4", "2x2x3"), regulator,
                          nsims=1E5, details=FALSE, setseed=TRUE)
{
  # function is no longer visible via NAMESPACE
  # .Deprecated(new="power.scABEL", msg=paste0("'power.scABEL2' is deprecated.",
  #             " Use 'power.scABEL' instead."))
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
  
  ptm <- proc.time()
  
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
  if(is.null(pe_constr)) pe_constr <- TRUE
  
  design <- match.arg(design)
  if (design=="2x3x3") {
    seqs  <- 3
    dfe   <- parse(text="n-3", srcfile=NULL)
    # this is only for testing purposes
    dfRRe <- parse(text="n-3", srcfile=NULL)
    #sd2  <- s2D + (s2wT + s2wR)/2 # used in v1.1-00 - v1.1-02
    # according to McNally et al.
    # verified via simulations:
    Emse  <- s2wT + s2wR/2
  }
  if (design=="2x2x4") {
    seqs  <- 2
    dfe   <- parse(text="n-2", srcfile=NULL)
    dfRRe <- parse(text="n-2", srcfile=NULL)
    # sd^2 of the differences T-R from their components
    Emse  <- (s2wT + s2wR)/2 
  }
  if (design=="2x2x3") {
    seqs  <- 2
    dfe   <- parse(text="n-2", srcfile=NULL)
    dfRRe <- parse(text="n/2-1", srcfile=NULL) # for balanced designs
    # sd^2 of the differences T-R from their components
    Emse  <- 1.5*(s2wT + s2wR)/2               # for balanced design 
  }
  
  if (length(n)==1){
    # then we assume n=ntotal
    # for unbalanced designs we divide the ns by ourself
    # to have only small imbalance
    nv <- nvec(n=n, grps=seqs)
    if (nv[1]!=nv[length(nv)]){
      message("Unbalanced design. n(i)=", paste(nv, collapse="/"), " assumed.")
    } 
  } else {
    # then we assume n = vector of n's in sequences
    # check length
    if (length(n)!=seqs) stop("n must be a vector of length=",seqs,"!")
    
    nv <- n
  }
  C3 <- sum(1/nv)/seqs^2
  n  <- sum(nv)
  
  df   <- eval(dfe)

  if (design=="2x2x3"){
    dfTT <- nv[1]-1
    dfRR <- nv[2]-1
    # where does this next came from?
    Emse <- (dfRR*(s2wT + s2wR/2)+dfTT*(s2wT/2 + s2wR))/(dfRR+dfTT)
    # warning in case of unbalanced design and heteroscdasticity
    # here not? TODO: check this against subject sims
    # if (abs(s2wT - s2wR)>1e-5 & abs(dfRR-dfTT)>2){
    #   warning(paste("Heteroscedasticity in unbalanced 2x2x3 design may led",
    #           "to poor accuracy of power!"), call.=FALSE)
    # }
  } else {
    dfRR <- eval(dfRRe)
  }
  
  # sd of the mean T-R (point estimate)
  sdm  <- sqrt(Emse*C3)
  # pe in the log domain
  mlog <- log(theta0)
  
  if(setseed) set.seed(123456)
  p <- .pwr.ABEL.ISC(mlog, sdm, C3, Emse, df, s2wR, dfRR, nsims, 
                     CVswitch=CVswitch, r_const=r_const, CVcap=CVcap, 
                     pe_constr=pe_constr, ln_lBEL=log(theta1),ln_uBEL=log(theta2), 
                     alpha=alpha)
    
  if (details) {
    ptm <- summary(proc.time()-ptm)
    message(nsims," sims. Time elapsed (sec): ", 
            formatC(ptm["elapsed"], digits=2), "\n")
    #print(ptm)
    # return also the components
    names(p) <- c("p(BE)", "p(BE-ABEL)", "p(BE-pe)", "p(BE-ABE)")
    if (!pe_constr) p <- p[-3] # without pe constraint
    p
  } else {
    # return only the 'power'
    as.numeric(p["BE"])
  }
}
