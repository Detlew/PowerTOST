# ----------------------------------------------------------------------------
# simulate subject data & evaluate via model with group effect
# gmodel==1 is full FDA model but only for testing grptreatment interaction
# gmodel==2 is full FDA model but without group by treatment interaction
# gmodel==3 is model with pooled groups, i.e. without any group term
#
# Author D. Labes based on power.scABEL.sdsims()
# ----------------------------------------------------------------------------
power.TOST.sds <- function(alpha=0.05, theta1, theta2, theta0, CV, n,   
                           design=c("2x2", "2x3x3", "2x2x4", "2x2x3"), 
                           design_dta=NULL, grps=2, ngrp = NULL, gmodel=2,
                           nsims=1E5, details=FALSE, setseed=TRUE, progress)
{
  # Check CV
  if (missing(CV)) stop("CV must be given!")
  # check theta0 and ABE limits
  if (missing(theta0)) theta0 <- 0.95
  if (length(theta0)>1) {
    theta0 <- theta0[2]
    warning(paste0("theta0 has to be scalar. theta0 = ",
                   theta0, " used."), call. = FALSE)
  }
  if (missing(theta1) & missing(theta2)) theta1 <- 0.8
  if (missing(theta2)) theta2 <- 1/theta1
  if (missing(theta1)) theta1 <- 1/theta2
  
  # CV scalar or vector with length=2
  # is this here reasonable?
  CVwT <- CV[1]
  if (length(CV)==2) CVwR <- CV[2] else CVwR <- CVwT
  # intra-subject variabilities from CV
  s2wT <- CV2mse(CVwT)
  s2wR <- CV2mse(CVwR)
  
  
  if (is.null(design_dta)){
    # check design
    desi <- match.arg(design)
    # degrees of freedom no longer set here
    # will be determined in the working horse via a first, single call of lm()
    if(desi=="2x2"){
      bk <- 2
      bkni <- 1/2
      seqs <- c("TR", "RT")
    }
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
    design_dta <- prep_data3(seqs, nseq=nv, muR=log(10), ldiff=log(theta0), 
                             s2wT=s2wT, s2wR=s2wR)
  } else {
    # check the data.frame design_dta 
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
    # above for definition via design and n it is feed to the working horse 
    # for comparative purposes only
    C2   <- NULL
  }
  # groups, ngrp= number if subjects in groups
  if (is.null(ngrp)){
    ngrp <- nvec(n=n, grps=grps)
  } else {
    # check the giving in ngrp
    # TODO
  }
  
  # grp no of subjects
  subs <- vector("numeric")
  for (i in 1:grps){
    subs <- c(subs, rep.int(i, ngrp[i]))
  }
  design_dta$grp <- 1
  for(i in seq_along(design_dta$subject)){
    design_dta$grp[i] <- subs[design_dta$subject[i]]
  }

  # progressbar or not
  if(missing(progress)) {
    progress <- FALSE
    if(nsims>=5E5) progress <- TRUE
    if(nsims>=1E5 & n>72) progress <- TRUE #???
    if(gmodel==1) progress <- TRUE
  }
  
  # call the working horse
  pwr <- .pwr.ABE.sds(muR=log(10), design_dta=design_dta, ldiff=log(theta0), 
                      s2WR=s2wR, s2WT=s2wT, C2=C2, nsims=nsims, gmodel = gmodel,  
                      ln_lBEL=log(theta1), ln_uBEL=log(theta2), alpha=alpha, 
                      setseed=setseed, details=details, progress=progress)
  pwr
}

# alias
power.ABE.sds <- power.TOST.sds

# ------ working horse ----------------------------------------------------
.pwr.ABE.sds <- function(muR=log(10), design_dta, ldiff, s2WR, s2WT, C2, 
                        nsims, ln_lBEL, ln_uBEL, alpha=0.05, gmodel, p.level=0.1,
                        setseed=TRUE, details=FALSE, progress=FALSE)
{
  # start time measurement
  ptm <- proc.time()
  
  if(progress) pb <- txtProgressBar( min = 0, max = 1, style = 3)

  if(setseed) set.seed(123456)

  dta <- design_dta
  # make a first simulation of logval
  dta$logval <- sim_data2_y(data_tmt=dta$tmt, ldiff=ldiff, s2wT=s2WT, s2wR=s2WR)
  dta$tmt      <- as.factor(dta$tmt)
  dta$period   <- as.factor(dta$period)
  dta$subject  <- as.factor(dta$subject)
  dta$grp      <- as.factor(dta$grp)
  dta$sequence <- as.factor(dta$sequence)
  # change coding to 1=R, 2=T
  dta_tmt <- as.numeric(dta$tmt)
  # measurements under treatments, values to simulate
  nT <- length(dta_tmt[dta_tmt==2])
  nR <- length(dta_tmt[dta_tmt==1])

  oc <- options(contrasts=c("contr.treatment","contr.poly"))

  # make a first evaluation via lm() to obtain degrees of freedom
  if(gmodel==2 || gmodel==1) {
    mudel <- lm(logval ~ tmt + grp + sequence + subject%in%(grp:sequence) 
                       + period%in%grp + grp:sequence, data=dta)
  }
  if(gmodel==3) {
    #mudel <- lm(logval ~ tmt + period + sequence + subject%in%sequence, data=dta)
    # without decomposiotion of the subject effect into sequence & subjet in sequence
    mudel <- lm(logval ~ tmt + period + subject, data=dta)
  }
  df     <- mudel$df.residual
  # save the model matrices for reuse in the simulation loop
  mm <- model.matrix(mudel)
  
  # calculate C2 for the equation CI half width hw = tcrit*sqrt(C2*mse) 
  # may be different to argument C2 in case of missing data
  # in case of alpha=0 C21 becomes NaN (tcrit=Inf, hw1=Inf)
  # therefore we calculate C21 with alpha=0.05
  tcrit <- qt(1-0.05, df)
  hw1   <- as.numeric(confint(mudel, level=1-2*0.05)["tmtT", 2]-coef(mudel)["tmtT"])
  mse1  <- summary(mudel)$sigma^2
  C21   <- (hw1/tcrit)^2/mse1
  C2    <- C21
  # now the correct tcrit to be used, df from the ANOVA of T-R
  tcrit <- qt(1-alpha, df)
  
  # allocate memory space for the results of the simulation loop
  pes   <- vector(mode="numeric", length=nsims)
  mses  <- vector(mode="numeric", length=nsims)
  BE_ABE <- vector(mode="logical", length=nsims)

  # yet another set.seed because the first sim is now done in prep_data2 above
  if(setseed) set.seed(123456)
  
  # working with multiple right-hand side (no_rhs) logvals 
  # attention! the code breaks if no_rhs = 1. then the returns of .lm.fit
  # or qr.coef are no longer matrices
  no_rhs <- 1
  if(gmodel!=1) no_rhs <- 500
  # at least nsims sims
  nsi <- ceiling(nsims/no_rhs)
  nsims <- nsi*no_rhs
  j1  <- 1
  j2  <- no_rhs
  # qr decomposition saved for re-use in the simulation loop
  qr_all <- qr(mm)
  
  p.GxT <- vector(mode="numeric", length=nsims) # p-vals of grp by tmt interaction
  gmod  <- vector(mode="numeric", length=nsims) # which model after check if significance
  #browser()
  if (gmodel==1) {
    # determine largest group(s?)
    largest <- as.numeric(which(summary(dta$grp) == max(summary(dta$grp))))
    if (length(largest)>1) {
      warning("More than 1 max. group not implemented yet.")
      largest <- largest[1]
    }
    dta3G <- dta[dta$grp==largest, ]
    m3G <- lm(logval ~ tmt + period + subject, data=dta3G)
    df3G <- m3G$df.residual
    t3G <- qt(1-alpha, df3G)
    # save model matrix and QR decomposition
    mm3G <- model.matrix(m3G)
    qr3G <- qr(mm3G)
  }
  # ---------------------------------------------------------------------------
  # loop of simulations
  for(j in 1:nsi){
    logval <- sim_mrhs(data_tmt=dta_tmt, nT=nT, nR=nR, ldiff=ldiff, 
                       s2wT=s2WT, s2wR=s2WR, no=no_rhs)
    if(gmodel==1){
      # we are working with only 1 right-hand-side!
      dta$logval <- logval[,1]
      # test the group by tmt interaction
      mud1 <- lm(logval ~ grp + sequence + tmt +
                          subject%in%(grp*sequence) + period%in%grp +
                          grp:sequence + grp:tmt, data=dta)
      p.GxT[j] <- anova(mud1)[["grp:tmt", "Pr(>F)"]]
      if(p.GxT[j] >= p.level){
        # interaction not significant
        # use model 2
        gmod[j1:j2] <- 2
        coefs <- qr.coef(qr_all, logval)
        pes[j1:j2]  <- coefs["tmtT", ]
        mses[j1:j2] <- colSums((qr.resid(qr_all, logval))^2)/df  
        # standard error of the difference T-R
        seD  <- sqrt(C2*mses[j1:j2])
        # ABE test = 1-2*alpha CI, df are the df of the ANOVA
        hw   <- tcrit*seD
        loCL <- pes[j1:j2] - hw
        upCL <- pes[j1:j2] + hw
        # conventional ABE decision
        BE_ABE[j1:j2] <- (loCL >=  ln_lBEL) & (upCL <= ln_uBEL)
      } else {
        # interaction significant, not allowed to pool
        # use model 3 with data of group with max. group size
        #browser()
        # if logval isn't a matrix qr.coeff only returns a named numeric vector
        logval <- as.matrix(logval[dta$grp==largest, 1])
        gmod[j1:j2] <- 3
        coefs <- qr.coef(qr3G, logval)
        pes[j1:j2]  <- coefs["tmtT", ]
        mses[j1:j2] <- colSums((qr.resid(qr3G, logval))^2)/df  
        # standard error of the difference T-R
        seD  <- sqrt(C2*mses[j1:j2])
        # ABE test = 1-2*alpha CI, df are the df of the ANOVA
        hw   <- t3G*seD
        loCL <- pes[j1:j2] - hw
        upCL <- pes[j1:j2] + hw
        # conventional ABE decision
        BE_ABE[j1:j2] <- (loCL >=  ln_lBEL) & (upCL <= ln_uBEL)

      }
    } else {
      # model 2 and 3
      coefs <- qr.coef(qr_all, logval)
      pes[j1:j2]  <- coefs["tmtT", ]
      # astonishing enough the next line gives NA only in case of gmodel==2
      # due to not full rank of the model matrix?
      #mses[j1:j2] <- colSums((logval - mm %*% coefs)^2)/df
      # thus we have to use qr.resid() for obtaining mse
      mses[j1:j2] <- colSums((qr.resid(qr_all, logval))^2)/df  
      # standard error of the difference T-R
      seD  <- sqrt(C2*mses[j1:j2])
      # ABE test = 1-2*alpha CI, df are the df of the ANOVA
      hw   <- tcrit*seD
      loCL <- pes[j1:j2] - hw
      upCL <- pes[j1:j2] + hw
      # conventional ABE decision
      BE_ABE[j1:j2] <- (loCL >=  ln_lBEL) & (upCL <= ln_uBEL)
    }

    # show progress
    if(progress){
      jsim <- j*no_rhs
      if(100*trunc(jsim/100)==jsim) setTxtProgressBar(pb, jsim/nsims)
    }  
    j1 <- j1 + no_rhs
    j2 <- j2 + no_rhs
  } 
  # reset options
  options(oc)
  
  # done with the progressbar
  if(progress) close(pb)
  
  pwr <- sum(BE_ABE)/nsims

  if (details){
    # return run-time
    ptm <- summary(proc.time()-ptm)
    tunit <- "sec"
    if(ptm["elapsed"]>60){
      ptm <- ptm/60; tunit <- "min"
    }
    message(nsims," sims. Time elapsed (",tunit,"): ", 
            formatC(ptm["elapsed"], digits=3), "\n")
    if (gmodel!=1) {
      # return only pwr
      pwr
    } else {
      #return pwr + some summary givings
      res <- list(pBE=pwr, 'p.GxT > p.level'=sum(p.GxT>p.level)/nsims)
      df <- data.frame(gmodels=gmod, BE=BE_ABE)
      res$xtab <- table(df)
      res
    }
  } else {
    pwr
  }
}
