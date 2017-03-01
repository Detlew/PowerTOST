# ----------------------------------------------------------------------------
# simulate replicate design subject data and evaluate via EMA ABEL method
# ANOVA & average BE with expanding limits
#
# Author D. Labes with suggestions by Ben
# ----------------------------------------------------------------------------
power.scABEL.sdsims <- function(alpha=0.05, theta1, theta2, theta0, CV, n,   
                                design=c("2x3x3", "2x2x4", "2x2x3"), regulator,
                                nsims=1E5, details=FALSE, setseed=TRUE, progress,
                                fitmethod=c(".lm.fit", "lm.fit", "fastlm", "qr"))
{
  fitm <- match.arg(fitmethod)
  
  if(missing(progress)) {
    progress <- FALSE
    if(nsims>1E4) progress <- TRUE
  }
  # check if on windows, if not don't show progressbar TODO
  
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
  # check if RcppEigen is installed
  # call the working horse
  pwr <- .pwr.ABEL.sdsims(seqs=seqs, nseq=nv, ldiff=log(theta0), s2WR=s2wR, 
                          s2WT=s2wT, C2=C2, df=df, dfRR=dfRR, nsims=nsims, 
                          regulator=reg, 
                          ln_lBEL=log(theta1), ln_uBEL=log(theta2), 
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
  pt <- proc.time()
  
  CVcap     <- regulator$CVcap
  CVswitch  <- regulator$CVswitch
  r_const   <- regulator$r_const
  pe_constr <- regulator$pe_constr
  # paranoia
  if(is.null(pe_constr)) pe_constr <- TRUE
  
  if(progress) pb <- winProgressBar(title = "sims progress", min = 0, max = nsims, width = 400)

  if(setseed) set.seed(123456)
#  set.seed(146389) # seed for the scripts in directory /workspace/replicate_simul
  
  tcrit <- qt(1-alpha, df)
  
  data <- prep_data2(seqs, nseq, muR=log(10), ldiff=ldiff, s2wT=s2WT, s2wR=s2WR)
  data$tmt     <- as.factor(data$tmt)
  data$period  <- as.factor(data$period)
  data$subject <- as.factor(data$subject)
  
  data_tmt <- data$tmt
  logval   <- data$logval

  oc <- options(contrasts=c("contr.treatment","contr.poly"))
  # save the model matrices for reuse in the simulation loop
  # the inclusion of sequence doesn't change the residual ms
  # model.matrix full
  mm <- model.matrix(logval~tmt+period+subject, data=data)
  # model matrix for R data only
  dataR   <- data[data$tmt=="R",]
  logvalR <- dataR$logval
  mmR     <- model.matrix(logval~period+subject, data=dataR)
  
  if(is.finite(CVcap)) wABEL_cap <- r_const*CV2se(CVcap) else wABEL_cap <- Inf
  s2switch <- CV2mse(CVswitch)
  s2cap    <- CV2mse(CVcap)
  
  # reserve space for the results of the simulation loop
  pes   <- vector(mode="numeric", length=nsims)
  mses  <- vector(mode="numeric", length=nsims)
  s2wRs <- vector(mode="numeric", length=nsims)
  # programming each fitmethod in own loop to avoid 1 Mio if's. but this give 
  # no notable difference in run-time
  if(fitmethod=="fastlm"){
    # or better work in chunks?
    for(j in 1:nsims){
      # use RcppEigen
      # Ben's method=2L gives different results to lm.fit() for
      # design="2x2x4" and "2x2x3"
      # but method=0L has no big advantage compared to lm.fit()
      model   <- RcppEigen::fastLmPure(X=mm, y=logval, method=2L)
      mses[j] <- (model$s)^2
      pes[j]  <- model$coefficients["tmtT"]
      
      # R data only, rank deficient in case of design="2x2x4" and "2x2x3"?
      modelR   <- RcppEigen::fastLmPure(X=mmR, y=logvalR, method=0L)
      s2wRs[j] <- (modelR$s)^2
      
      # show progress
      if(progress){
        if(100*trunc(j/100)==j) 
          setWinProgressBar(pb, j, title=paste(j, "sims done (",
                                               trunc(1000*j/nsims)/10,"%)"))
      }  
      # simulate for next step
      # this gives one step more as necessary but avoiding here an if may have some
      # run-time advantage
      logval  <- sim_data2_y(data_tmt=data_tmt, ldiff=ldiff, s2wT=s2WT, s2wR=s2WR)
      logvalR <- logval[data_tmt == "R"]
    }  
  } else if(fitmethod=="lm.fit") {
    # using lm.fit()
    for(j in 1:nsims){
      # determine pe, mse of all data
      model   <- lm.fit(x=mm, y=logval)
      # lm.fit does'nt have anova() and confint() methods
      # thus we have to calculate them explicite
      mses[j] <- sum(model$residuals^2)/model$df.residual
      pes[j]  <- coef(model)["tmtT"]
      
      # determine s2wR from R data only
      modelR   <- lm.fit(x=mmR, y=logvalR)
      s2wRs[j] <- sum(modelR$residuals^2)/modelR$df.residual
      
      # show progress
      if(progress){
        if(100*trunc(j/100)==j) 
          setWinProgressBar(pb, j, title=paste(j, "sims done (",
                                               trunc(1000*j/nsims)/10,"%)"))
      }  
      # simulate for next step
      # this gives one step more as necessary but avoiding here an if may have some
      # run-time advantage
      logval  <- sim_data2_y(data_tmt=data_tmt, ldiff=ldiff, s2wT=s2WT, s2wR=s2WR)
      logvalR <- logval[data_tmt == "R"]
    }  
  } else if(fitmethod==".lm.fit") {
    # using lm.fit()
    for(j in 1:nsims){
      # determine pe, mse of all data
      model   <- .lm.fit(x=mm, y=logval)
      # lm.fit does'nt have anova() and confint() methods
      # thus we have to calculate them explicite
      mses[j] <- sum(model$residuals^2)/df
      pes[j]  <- model$coefficients[2]
      
      # determine s2wR from R data only
      modelR   <- .lm.fit(x=mmR, y=logvalR)
      s2wRs[j] <- sum(modelR$residuals^2)/dfRR
      
      # show progress
      if(progress){
        if(100*trunc(j/100)==j) 
          setWinProgressBar(pb, j, title=paste(j, "sims done (",
                                               trunc(1000*j/nsims)/10,"%)"))
      }  
      # simulate for next step
      # this gives one step more as necessary but avoiding here an if may have some
      # run-time advantage
      logval  <- sim_data2_y(data_tmt=data_tmt, ldiff=ldiff, s2wT=s2WT, s2wR=s2WR)
      logvalR <- logval[data_tmt == "R"]
    }  
  } else {
    # qr decomp saved for re-use
    qr_all <- qr(mm)
    qr_R   <- qr(mmR)
    for(j in 1:nsims){
      pes[j]   <- qr.coef(qr_all, logval)["tmtT"]
      mses[j]  <- sum((qr.resid(qr_all, logval))^2)/df
      
      s2wRs[j] <- sum((qr.resid(qr_R, logvalR))^2)/dfRR
      
      # show progress
      if(progress){
        if(100*trunc(j/100)==j) 
          setWinProgressBar(pb, j, title=paste(j, "sims done (",
                                               trunc(1000*j/nsims)/10,"%)"))
      }  
      # simulate for next step
      # this gives one step more as necessary but avoiding here an if may have some
      # run-time advantage
      logval  <- sim_data2_y(data_tmt=data_tmt, ldiff=ldiff, s2wT=s2WT, s2wR=s2WR)
      logvalR <- logval[data_tmt == "R"]
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
  
  if (details){
    cat("Using", fitmethod, "\n")
    cat("Time consumed [min]:\n")
    print(round((proc.time()-pt)/60,2));cat("\n")
      
    cat("pBE(sABE)      =", sum(BE_sc)/nsims, "\n")
    cat("pBE(ABE)       =", sum(BE_ABE)/nsims, "\n")
    cat("pBE(pe constr.)=", sum(BE_pe)/nsims, "\n")
    cat("\n")
  }
  # result: probability of study passed:
  # CI within wABEL and pe in acceptance range
  sum(BE)/nsims
}
