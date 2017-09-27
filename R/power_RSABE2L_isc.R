#---------------------------------------------------------------------------
# Simulate partial or full replicate designs and scaled ABE power
# estimation method via intra-subject contrasts of T-R and R-R
# BE decision via "exact" method of the 2L, ABEL or hyslop's method
#  
# no cap if regulator = "FDA", else CVcap is taken into account if defined
# in regulator argument of class 'regSet'
# 
# Author: dlabes
#---------------------------------------------------------------------------

# degrees of freedom for the TR/RR  analysis: 
# Using the intrasubject contrasts T-R and R-R and analyze them  
# by sequence groups the df's are = n-seq.
# 2x3x3  dfRR = n-3
# 2x2x4  dfRR = n-2
# 2x2x3  dfRR = n/2 - 2

power.RSABE2L.isc <- function(alpha=0.05, theta1, theta2, theta0, CV, n,   
                              design=c("2x3x3", "2x2x4", "2x2x3"), regulator,
                              SABE_test="exact",
                              nsims=1E5, details=FALSE, setseed=TRUE)
{
  if (missing(CV)) stop("CV must be given!")
  if (missing(n))  stop("Number of subjects n must be given!")

  if (missing(theta0)) theta0 <- 0.90
  if (missing(theta1) & missing(theta2)) theta1 <- 0.8
  if (missing(theta2)) theta2 <- 1/theta1
  if (missing(theta1)) theta1 <- 1/theta2
  
  ptm <- proc.time()
  
  # for later enhancement taking into account the 
  # subject-by-formulation interaction
  s2D  <- 0 
  CVwT <- CV[1]
  if (length(CV)==2) CVwR <- CV[2] else CVwR <- CVwT
  s2wT <- CV2mse(CVwT)
  s2wR <- CV2mse(CVwR)

  # regulator here only FDA, EMA
  # other regulatory bodies ("HC", "ANVISA") use all the EMA regulatory constant
  if (missing(regulator)) regulator <- "FDA"
  rc <- reg_check(regulator, choices=c("FDA", "EMA"))
  CVswitch  <- rc$CVswitch
  r_const   <- rc$r_const
  pe_constr <- rc$pe_constr
  # CVcap doesn't apply to the FDA recommended method
  # but in Munoz et al. method = Howe-EMA
  CVcap     <- rc$CVcap

  # check SABE_test
  SABE_test <- tolower(SABE_test)
  SABE_test <- match.arg(SABE_test, choices=c("exact", "abel", "hyslop"))
  
  # check design and give the design characteristics
  design <- match.arg(design)
  if (design=="2x3x3") {
    seqs  <- 3
    dfe   <- parse(text="n-3", srcfile=NULL)
    dfRRe <- parse(text="n-3", srcfile=NULL)
    #sd2  <- s2D + (s2wT + s2wR)/2 # used in v1.1-00 - v1.1-02
    # according to McNally et al.
    # verified via simulations:
    Emse  <- s2D + s2wT + s2wR/2
  }
  if (design=="2x2x4") {
    seqs  <- 2
    dfe   <- parse(text="n-2", srcfile=NULL)
    dfRRe <- parse(text="n-2", srcfile=NULL)
    # sd^2 of the differences T-R from their components
    Emse  <- (s2D + (s2wT + s2wR)/2) 
  }
  if (design=="2x2x3") {
    seqs  <- 2
    dfe   <- parse(text="n-2", srcfile=NULL)
    # next was pre-V1.2-08
#     dfRRe <- parse(text="n/2-2", srcfile=NULL) # for balanced designs
#     dfTTe <- parse(text="n/2-2", srcfile=NULL) # for balanced designs
    # correct should be (only 1 sequence for each, f.i. RR from RTR):
    dfRRe <- parse(text="n/2-1", srcfile=NULL) # for balanced designs
    dfTTe <- parse(text="n/2-1", srcfile=NULL) # for balanced designs
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
    C3 <- sum(1/nv)/seqs^2
    n  <- sum(nv)
  } else {
    # then we assume n = vector of n's in sequences
    # check length
    if (length(n)!=seqs) stop("n must be a vector of length=",seqs,"!")
    
    C3 <- sum(1/n)/seqs^2
    nv <- n
    n  <- sum(n)
  }
  
  df   <- eval(dfe)
  if (design=="2x2x3"){
    dfTT <- nv[1]-1
    dfRR <- nv[2]-1
    # where does this next came from?
    Emse <- (dfRR*(s2wT + s2wR/2)+dfTT*(s2wT/2 + s2wR))/(dfRR+dfTT)
    # warning in case of unbalanced design and heteroscdasticity
    if (abs(s2wT - s2wR)>1e-5 & abs(dfRR-dfTT)>2){
      warning(paste("Heteroscedasticity in unbalanced 2x2x3 design may lead", 
              "to poor accuracy of power!"), call.=FALSE)
    }
  } else {
    dfRR <- eval(dfRRe)
  }
  #cat("dfRR=", dfRR," dfTT=",dfTT," E(s2I)=", Emse, "\n")
  # sd of the mean T-R (point estimator)
  sdm  <- sqrt(Emse*C3)
  mlog <- log(theta0)
  
  if(setseed) set.seed(123456)
  p <- .power.RSABE.isc(mlog, sdm, C3, Emse, df, s2wR, dfRR, nsims, 
                        CVswitch, r_const, pe_constr, CVcap, SABE_test=SABE_test,
                        ln_lBEL=log(theta1),ln_uBEL=log(theta2), alpha=alpha)
    
  if (details) {
    ptm <- summary(proc.time()-ptm)
    message(nsims," sims. Time elapsed (sec): ", 
            formatC(ptm["elapsed"], digits=2), "\n")
    #print(ptm)
    # return also the components
    names(p) <- c("p(BE)", "p(BE-RSABE)", "p(BE-pe)", "p(BE-ABE)")
    if (SABE_test=="abel") names(p) <- c("p(BE)", "p(BE-ABEL)", "p(BE-pe)", "p(BE-ABE)")
    if (!pe_constr) p <- p[-3] # without pe constraint
    p
  } else {
    # return only the 'power'
    as.numeric(p["BE"])
  }
}

# working horse of RSABE
.power.RSABE.isc <- function(mlog, sdm, C3, Emse, df, s2wR, dfRR, nsims, 
                             CVswitch, r_const, pe_constr, CVcap, SABE_test="exact",
                             ln_lBEL, ln_uBEL, alpha=0.05)
{
  tval     <- qt(1-alpha,df)
  chisqval <- qchisq(1-alpha, dfRR)
  r2const  <- r_const^2
  s2switch <- log(CVswitch^2+1) 
  s2cap    <- log(CVcap^2+1)
  if(is.finite(CVcap)) wABEL_cap <- r_const*CV2se(CVcap) else wABEL_cap <- Inf
  
  counts <- rep.int(0, times=4)
  names(counts) <- c("BE", "BEul", "BEpe", "BEabe")
  # to avoid memory problems for high number of sims we work in chunks
  chunks <- 1
  nsi    <- nsims
  if (nsims>1E7) {
    chunks <- ceiling(nsims/1E7)
    nsi    <- 1E7
  } 
  for (iter in 1:chunks) {
    # if chunks*1E7 >nsims correct nsi to given nsims
    if(iter==chunks) nsi <- nsims-(chunks-1)*nsi
    # debug code
    # cat("nsi=", nsi, "\n")
    # simulate sample mean via its normal distribution
    means  <- rnorm(nsi, mean=mlog, sd=sdm)
    # simulate sample sd2s via chi-square distri
    sd2s   <- Emse*C3*rchisq(nsi, df)/df
    # simulate sample value s2wRs via chi-square distri
    s2wRs  <- s2wR*rchisq(nsi, dfRR)/dfRR
    
    SEs <- sqrt(sd2s)
    # conventional (1-2*alpha) CI's for T-R
    hw  <- tval*SEs
    lCL <- means - hw 
    uCL <- means + hw
    # conventional ABE
    BE_ABE <- ((ln_lBEL<=lCL) & (uCL<=ln_uBEL))
    
    if (SABE_test=="exact"){
      # step 1: compute k
      k <- SEs/sqrt(s2wRs)
      # Hedges correction
      Hf <- 1-3/(4*dfRR-1)
      # step 2: compute L/U using eqn. (26)
      # warnings wrt to full precision in pnt
      op <- options()
      options(warn=-1)
        Ltheta <- qt(1-alpha, df, -(Hf/k)*r_const)
        Utheta <- qt(alpha, df, +(Hf/k)*r_const)
      options(op)
      # effect size
      es <- means/sqrt(s2wRs)/k
      # RSABE ("exact") decision
      BE_RSABE  <- (Ltheta < es) & (es < Utheta)
      
    } else if (SABE_test=="hyslop" | SABE_test=="howe") {
      # linearized SABE criterion + 95%CI
      # with -SEs^2 the 'unknown' x from the progesterone guidance
      Em <- means^2 #- SEs^2  
      Es <- r2const*s2wRs
      Cm <- (abs(means) + hw)^2
      Cs <- Es*dfRR/chisqval    
      SABEc95 <- Em - Es + sqrt((Cm-Em)^2 + (Cs-Es)^2)
      # save memory
      rm(SEs, hw, Em, Es, Cm, Cs)
      BE_RSABE <- SABEc95<=0
    } else {
      # ABEL method
      # 'widened' acceptance limits or conventional limits below s2switch
      wABEL  <- ifelse(s2wRs<=s2switch, ln_uBEL, r_const*sqrt(s2wRs))
      # cap on widening
      wABEL  <- ifelse(s2wRs>s2cap, wABEL_cap, wABEL)
      # scaled ABE (ABEL) decision
      BE_RSABE  <- (lCL >= -wABEL) & (uCL <= wABEL)
    }
    # mixed procedure
    BE <- ifelse(s2wRs>s2switch, BE_RSABE, BE_ABE)
    # use capped acceptance limits if CVwR > CVcap
    # already incorporated with "abel"
    if (is.finite(CVcap) & SABE_test!="abel"){
      s2cap <- CV2mse(CVcap)
      # calculate the capped widened acceptance limits in log domain
      uprABEL <- r_const*sqrt(s2cap)
      lwrABEL <- -uprABEL
      BE <- ifelse(s2wRs>=s2cap, ((lwrABEL<=lCL) & (uCL<=uprABEL)), BE)
    }
    
    # point est. constraint true?
    BEpe  <- ( means>=ln_lBEL & means<=ln_uBEL )

    counts["BEabe"] <- counts["BEabe"] + sum(BE_ABE)
    counts["BEpe"]  <- counts["BEpe"]  + sum(BEpe)
    counts["BEul"]  <- counts["BEul"]  + sum(BE)
    if(pe_constr) {
      counts["BE"]    <- counts["BE"]    + sum(BE & BEpe)
    } else {
      counts["BE"]    <- counts["BE"]    + sum(BE) # no pe constraint
    }
    
  } # end over chunks
  # return the pBEs
  counts/nsims
}