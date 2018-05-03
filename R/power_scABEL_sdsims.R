# ----------------------------------------------------------------------------
# simulate replicate design subject data and evaluate via EMA ABEL method
# ANOVA & average BE with expanding limits
#
# Author D. Labes with suggestions by Ben
# This is the version trying to optimize more things
# using Ben's function to simulate multiple right-hand sides
# Trying to made more self-contained, df and C2 from a first, single lm() call
# design definition now alternatively via data.frame 'design_dta'
# a data.frame with columns subject, sequence, period and tmt and optional logval
# sequence and tmt in R/T notation
# ----------------------------------------------------------------------------
power.scABEL.sdsims <- function(alpha=0.05, theta1, theta2, theta0, CV, n,   
                                design=c("2x3x3", "2x2x4", "2x2x3"), 
                                design_dta=NULL, regulator,
                                nsims=1E5, details=FALSE, setseed=TRUE, progress)
{
  # check theta0 and ABE limits
  if (missing(theta0)) theta0 <- 0.90
  if (length(theta0)>1) {
    theta0 <- theta0[2]
    warning(paste0("theta0 has to be scalar. theta0 = ",
                   theta0, " used."), call. = FALSE)
  }
  if (missing(theta1) & missing(theta2)) theta1 <- 0.8
  if (missing(theta2)) theta2 <- 1/theta1
  if (missing(theta1)) theta1 <- 1/theta2
  
  # check regulator
  if (missing(regulator)) regulator <- "EMA"
  reg  <- reg_check(regulator)
  if (reg$est_method=="ISC") stop("ISC evaluation not allowed in this function.")
  # Check CV
  if (missing(CV)) stop("CV must be given!")
  # CV scalar or vector
  CVwT <- CV[1]
  if (length(CV)==2) CVwR <- CV[2] else CVwR <- CVwT
  # intra-subject variabilities from CV
  s2wT <- CV2mse(CVwT)
  s2wR <- CV2mse(CVwR)
  
  if (is.null(design_dta)){
    # check design
    desi <- match.arg(design)
    # degrees of freedom no longer set here
    # will be determined in the working horse via a first, single call of lm
    if(desi=="2x3x3"){
      bk    <- 1.5
      bkni  <- 1/6
      seqs  <- c("TRR", "RTR", "RRT")
    }
    if(desi=="2x2x3"){
      bk    <- 1.5
      bkni  <- 3/8
      seqs  <- c("TRT", "RTR")
    }
    if(desi=="2x2x4"){
      bk    <- 1
      bkni  <- 1/4
      seqs  <- c("TRTR", "RTRT")
    }
    seqn <- length(seqs)
    # check n
    if (missing(n))  stop("Number of subjects n must be given!")
    # check n: vector or scalar
    if (length(n)==1){
      # for unbalanced designs we divide the ns by ourself
      # to have 'smallest' imbalance
      nv <- nvec(n=n, grps=seqn)
      if (nv[1]!=nv[length(nv)]){
        message("Unbalanced design. n(i)=", paste(nv, collapse="/"), " assumed.")
      }
      C2 <- sum(1/nv)*bkni
      n <- sum(nv)
    } else {
      # check length
      if (length(n)!=seqn) stop("n must be a vector of length=",seqn,"!")
      C2 <- sum(1/n)*bkni
      nv <- n
      n <- sum(n)
    }
    design_dta <- prep_data2(seqs, nseq=nv, muR=log(10), ldiff=log(theta0), 
                             s2wT=s2wT, s2wR=s2wR)
  } else {
    # check data.frame design_dta 
    # TODO
    
    # delete NA in logval if logval is defined in data.frame
    if ("logval" %in% names(design_dta)){
      design_dta <- design_dta[!is.na(design_dta$logval),]
    }
    seqs <- unique(design_dta$sequence)
    seqn <- length(seqs)
    nv   <- as.numeric(table(unique(design_dta[,c("subject", "sequence")])[,"sequence"]))
    n    <- sum(nv)
    # C2 is calculated in the working horse
    # above for definition via design and n it is given to the working horse 
    # for comparative purposes only
    C2   <- NULL
  }
  # progressbar or not
  if(missing(progress)) {
    progress <- FALSE
    if(nsims>=5E5) progress <- TRUE
    if(nsims>=1E5 & n>72) progress <- TRUE
  }
  
  # after introducing multiple right-hand sides .lm.fit seems to have no longer
  # an advantage in run-time, fitmethod removed May 2018
  
  # call the working horse
  pwr <- .pwr.ABEL.sds(design_dta=design_dta, ldiff=log(theta0), 
                       s2WR=s2wR, s2WT=s2wT, C2=C2, nsims=nsims, regulator=reg, 
                       ln_lBEL=log(theta1), ln_uBEL=log(theta2), alpha=alpha, 
                       setseed=setseed, details=details, progress=progress)
  pwr
}
