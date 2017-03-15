# ----------------------------------------------------------------------------
# simulate replicate design subject data and evaluate via EMA ABEL method
# ANOVA & average BE with expanding limits
#
# Author D. Labes with suggestions by Ben
# This is the version trying to optimize more things
# using Ben's function to simulate multiple right-hand sides
# ----------------------------------------------------------------------------
power.scABEL.sdsims <- function(alpha=0.05, theta1, theta2, theta0, CV, n,   
                                design=c("2x3x3", "2x2x4", "2x2x3"), regulator,
                                nsims=1E5, details=FALSE, setseed=TRUE, progress)
{
  if(missing(progress)) {
    progress <- FALSE
    if(nsims>=5E5) progress <- TRUE
  }

  # check design
  desi <- match.arg(design)
  # check regulator
  if (missing(regulator)) regulator <- "EMA"
  reg  <- reg_check(regulator)
  if (reg$est_method=="ISC") stop("ISC evaluation not allowed in this function.")

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
  if(desi=="2x3x3"){
    bk    <- 1.5
    bkni  <- 1/6
    dfc   <- "2*n-3"
    dfcRR <- "n-2"
    seqs  <- c("TRR", "RTR", "RRT")
  }
  if(desi=="2x2x3"){
    bk    <- 1.5
    bkni  <- 3/8
    dfc   <- "2*n-3"
    dfcRR <- "n/2-1"        #balanced design only, is set explicitely later
    seqs  <- c("TRT", "RTR")
  }
  if(desi=="2x2x4"){
    bk    <- 1
    bkni  <- 1/4
    dfc   <- "3*n-4"
    dfcRR <- "n-2"
    seqs  <- c("TRTR", "RTRT")
  }
  seqn <- length(seqs)
  # degrees of freedom as expression
  dfe   <- parse(text=dfc, srcfile=NULL)
  dfRRe <- parse(text=dfcRR, srcfile=NULL) 
  # CV scalar or vector
  CVwT <- CV[1]
  if (length(CV)==2) CVwR <- CV[2] else CVwR <- CVwT
  s2wT <- log(1.0 + CVwT^2)
  s2wR <- log(1.0 + CVwR^2)
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
  df <- eval(dfe)
  dfRR <- eval(dfRRe)
  if(desi=="2x2x3") dfRR <- nv[2]-1 
  
  # auto fit method: set this to .lm.fit if n<=18
  fitm <- "qr"
  # after introducing multiple right-hand sides .lm.fit seems to have no longer
  # an advantage in run-time
  # if(n<=18) fitm <- ".lm.fit"
  
  # call the working horse
  pwr <- .pwr.ABEL.sdsims(seqs=seqs, nseq=nv, ldiff=log(theta0), s2WR=s2wR, 
                          s2WT=s2wT, C2=C2, df=df, dfRR=dfRR, nsims=nsims, 
                          regulator=reg, ln_lBEL=log(theta1), ln_uBEL=log(theta2), 
                          alpha=alpha, fitmethod=fitm, setseed=setseed, 
                          details=details, progress=progress)
  pwr
}
  

# working horse
.pwr.ABEL.sdsims <- function(seqs, nseq, muR=log(10), ldiff, s2WR, s2WT, C2,
                             df, dfRR, nsims, regulator,
                             ln_lBEL=log(0.8), ln_uBEL=log(1.25), 
                             alpha=0.05, fitmethod="lm.fit", setseed=TRUE, 
                             details=FALSE, progress=FALSE)
{
  # start time measurement
  ptm <- proc.time()
  
  CVcap     <- regulator$CVcap
  CVswitch  <- regulator$CVswitch
  r_const   <- regulator$r_const
  pe_constr <- regulator$pe_constr
  # paranoia
  if(is.null(pe_constr)) pe_constr <- TRUE
  
  if(progress) pb <- txtProgressBar( min = 0, max = 1, style = 3)

  if(setseed) set.seed(123456)
#  set.seed(146389) # seed for the scripts in directory /workspace/replicate_simul
  
  tcrit <- qt(1-alpha, df)
  
  # prep_data2 gives all logvals = 0 back
  data <- prep_data2(seqs, nseq, muR=log(10), ldiff=ldiff, s2wT=s2WT, s2wR=s2WR)
  data$tmt     <- as.factor(data$tmt)
  data$period  <- as.factor(data$period)
  data$subject <- as.factor(data$subject)
  # change coding to 1=R, 2=T
  data_tmt <- as.numeric(data$tmt)
  # measurements under treatments, values to simulate
  nT <- length(data_tmt[data_tmt==2])
  nR <- length(data_tmt[data_tmt==1])

  oc <- options(contrasts=c("contr.treatment","contr.poly"))
  # save the model matrices for reuse in the simulation loop
  # the inclusion of sequence doesn't change the residual ms
  # model.matrix full
  mm <- model.matrix(~tmt+period+subject, data=data)
  # model.matrix for R data only
  dataR   <- data[data$tmt=="R",]
  mmR     <- model.matrix(~period+subject, data=dataR)
  
  if(is.finite(CVcap)) wABEL_cap <- r_const*CV2se(CVcap) else wABEL_cap <- Inf
  s2switch <- CV2mse(CVswitch)
  s2cap    <- CV2mse(CVcap)
  
  # allocate memory space for the results of the simulation loop
  pes   <- vector(mode="numeric", length=nsims)
  mses  <- vector(mode="numeric", length=nsims)
  s2wRs <- vector(mode="numeric", length=nsims)
  
  # working with multiple right-hand side logvals 
  # attention! the code breaks if no_rhs = 1. then the returns of .lm.fit
  # or qr.coef are no longer matrices
  no_rhs <- 500
  # at least nsims sims
  nsi <- ceiling(nsims/no_rhs)
  nsims <- nsi*no_rhs
  j1  <- 1
  j2  <- no_rhs
  # only ".lm.fit" and "qr" retained. if n<=18 ".lm.fit" is used else "qr"
  # programming each fitmethod in own loop to avoid 1 Mio if's. but this give 
  # no notable difference in run-time
  if(fitmethod==".lm.fit") {
    # using .lm.fit()
    for(j in 1:nsi){
      logval <- sim_mrhs(data_tmt=data_tmt, nT=nT, nR=nR, ldiff=ldiff, 
                         s2wT=s2WT, s2wR=s2WR, no=no_rhs)
      logvalR <- logval[data_tmt==1,]
      
      # determine pe, mse of all data
      model   <- .lm.fit(x=mm, y=logval)
      # lm.fit does'nt have anova() and confint() methods
      # thus we have to calculate them explicite
      mses[j1:j2] <- colSums(model$residuals^2)/df
      pes[j1:j2]  <- model$coefficients[2,]
      
      # determine s2wR from R data only
      modelR       <- .lm.fit(x=mmR, y=logvalR)
      s2wRs[j1:j2] <- colSums(modelR$residuals^2)/dfRR
      
      # show progress
      if(progress){
        jsim <- j*no_rhs
        if(100*trunc(jsim/100)==jsim) setTxtProgressBar(pb, jsim/nsims)
      }
      j1 <- j1+no_rhs
      j2 <- j2+no_rhs
      #browser()
    }  
  } else {
    # qr decomp saved for re-use
    qr_all <- qr(mm)
    qr_R   <- qr(mmR)
    for(j in 1:nsi){
      logval <- sim_mrhs(data_tmt=data_tmt, nT=nT, nR=nR, ldiff=ldiff, 
                         s2wT=s2WT, s2wR=s2WR, no=no_rhs)
      logvalR <- logval[data_tmt==1,]
      
      # original attempt
      # pes[j]   <- qr.coef(qr_all, logval)["tmtT"]
      # mses[j]  <- sum((qr.resid(qr_all, logval))^2)/df
      
      # Instead of qr.resid() we use faster approach via y - X * coefs
      coefs <- qr.coef(qr_all, logval)
      pes[j1:j2]  <- coefs["tmtT", ]
      mses[j1:j2] <- colSums((logval - mm %*% coefs)^2)/df
      
      # For reference: do not use this alternative approach as some
      # coefficients may be NA due to non-full rank
      s2wRs[j1:j2] <- colSums((qr.resid(qr_R, logvalR))^2)/dfRR
      
      # doing all by "hand"
      # coeffs  <- qr.coef(qr_all, logval)
      # resid   <- as.numeric(mm %*% coeffs - logval)
      # pes[j]  <- coeffs["tmtT"]
      # mses[j] <- sum(resid^2)/df
      # 
      # coeffsR  <- qr.coef(qr_R, logvalR)
      # residR   <- as.numeric(mmR %*% coeffsR - logvalR)
      # s2wRs[j] <- sum(residR^2)/dfRR
      
      # show progress
      if(progress){
        jsim <- j*no_rhs
        if(100*trunc(jsim/100)==jsim) setTxtProgressBar(pb, jsim/nsims)
      }  
      j1 <- j1+no_rhs
      j2 <- j2+no_rhs
    } 
  }
  
  # reset options
  options(oc)
  # vectorized calculations of CI's and widened acceptance limits 
  # (moved outside of loops)
  hw <- tcrit*sqrt(C2*mses)
  loCL <- pes - hw
  upCL <- pes + hw
  # conventional ABE decision
  BE_ABE <- (loCL >=  ln_lBEL) & (upCL <= ln_uBEL)
  # widened acceptance limits
  wABEL  <- ifelse(s2wRs<=s2switch, ln_uBEL, r_const*sqrt(s2wRs))
  # cap on widening
  wABEL  <- ifelse(s2wRs>s2cap, wABEL_cap, wABEL)
  # scaled ABE (ABEL) decision
  BE_sc  <- (loCL >= -wABEL) & (upCL <= wABEL)
  # pe constraint
  BE_pe <- (pes >= ln_lBEL & pes <= ln_uBEL)
  if(pe_constr) BE <- BE_sc & BE_pe else BE <- BE_sc
  
  # done with the progressbar
  if(progress) close(pb)
  
  # same output as power.scABEL
  p <- vector("numeric", length=4)
  names(p) <- c("p(BE)", "p(BE-wABEL)", "p(BE-pe)", "p(BE-ABE)")
  p[1] <- sum(BE)/nsims
  p[2] <- sum(BE_sc)/nsims
  p[3] <- sum(BE_pe)/nsims
  p[4] <- sum(BE_ABE)/nsims
  
  if (details){
    ptm <- summary(proc.time()-ptm)
    tunit <- "sec"
    if(ptm["elapsed"]>60){
      ptm <- ptm/60; tunit <- "min"
    }
    if (fitmethod!="qr") message("Using ", fitmethod, "\n")
    message(nsims," sims. Time elapsed (",tunit,"): ", 
            formatC(ptm["elapsed"], digits=3), "\n")
    
    if (!pe_constr) p <- p[-3] # without pe constraint
    p
  } else {
    as.numeric(p["p(BE)"])
  }
}
