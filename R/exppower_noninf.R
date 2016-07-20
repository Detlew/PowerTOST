#------------------------------------------------------------------------------
# Author: dlabes, Benjamin Lang
#------------------------------------------------------------------------------

# --- Approximate "expected" power according to Julious book
# taking into account the uncertainty of an estimated se with 
# dfse degrees of freedom
# Only for log-transformed data.
# Raw function: see the call in exppower.noninf()
.approx.exppower.noninf <- function(alpha=0.025, lmargin, diffm, sem, dfse, df)
{
  tval <- qt(1-alpha, df, lower.tail = TRUE)
  tau  <- sqrt((diffm-lmargin)^2/sem^2)
  # in case of diffm=lmargin and se=0
  # tau has the value NaN
  tau[is.nan(tau)] <- 0
  pow  <- pt(tau, dfse, tval)
  # values <0 not possible?
  pow
}

# --- Exact implementation of expected power
# Author(s) B. Lang & D. Labes
.exact.exppower.noninf <- function(alpha=0.025, lmargin, ldiff, se.fac, se, 
                                   dfCV, df) 
{
  # infinite df of CV should give expected power identical to power   
  if (!is.finite(dfCV)) {
    return(.power.noninf(alpha, lmargin, ldiff, se.fac*se, df))
  }
  # Define assurance function (expected power)
  f <- function(v) {
    .power.noninf(alpha, lmargin, ldiff, se.fac*sqrt(v), df) * 
         dinvgamma(x = v, shape = dfCV/2, scale = dfCV/2 * se^2)
  }
  
  # Detlew's attempt to stay with one integrate() call
  # modified by Ben
  
  # mean of inverse gamma with alpha=dfCV/2, beta=se^2*dfCV/2
  minvg <- (dfCV/2 * se^2)/(dfCV/2 - 1)
  # or should we use the mode?
  minvg <- (dfCV/2 * se^2)/(dfCV/2 + 1)
  # variance of inverse gamma
  vinvg <- (dfCV/2 * se^2)^2/(dfCV/2 - 1)^2/(dfCV/2 - 2)
  # Chebyshev's inequality with k = 10
  k <- 10
  lwr  <- max(0, minvg - k*sqrt(vinvg))
  # heavier tail to the right. may be an other value than 2 fits better TODO
  upr  <- minvg + 2*k*sqrt(vinvg)       
  pwr <- integrate(f, lwr, upr, rel.tol = 1e-5, stop.on.error = FALSE)
  if(pwr$message!="OK") warning(pwr$message)
  # return expected power
  pwr$value
  
}

# --- working horse for exppower.noninf() and expsampleN.noninf()
.exppower.noninf <- function(alpha=0.025, lmargin, ldiff, se.fac, se, dfCV, 
                             df, method)
{
  if (method == "exact") {
    return(.exact.exppower.noninf(alpha, lmargin, ldiff, se.fac, se, dfCV, df))
  } else if (method == "approx") {
    return(.approx.exppower.noninf(alpha, lmargin, ldiff, se.fac*se, dfCV, df))
  } else {
    # Paranoia. Is tested outside.
    stop("Method '", method, "' unknown!\n", call. = TRUE)
  }
}
# Main function for expected power
exppower.noninf <- function(alpha=0.025, logscale=TRUE, theta0, margin, 
                            CV, dfCV, n, design="2x2",  robust=FALSE,
                            method=c("exact", "approx")) 
{
  # check if design is implemented
  d.no <- .design.no(design)
  if (is.na(d.no)) stop("Design ",design, " unknown!", call.=FALSE)
  
  # design characteristics
  ades <- .design.props(d.no)
  #df as expression
  dfe  <- .design.df(ades, robust=robust)
  #design const.
  #bk   <- ades$bk
  
  if (missing(CV) | missing(dfCV)) stop("CV and df must be given!", call.=FALSE)
  if (any(CV<=0)) stop("CV has to be >0", call.=FALSE)
  if (length(CV)>1) {
    CV <- CV[1]
    warning("CV has to be a scalar. Only CV[1] used.", call.=FALSE)
  }
  if (length(dfCV)>1) {
    dfCV <- dfCV[1]
    warning("dfCV has to be a scalar here. Only dfCV[1] used.", call.=FALSE)
  }
  if (dfCV<=4) stop("dfCV has to be >4", call.=FALSE)
  
  if (missing(n)) stop("Number of subjects must be given!", call.=FALSE)
  
  if (logscale){
    if (missing(theta0)) theta0 <- 0.95
    if (missing(margin)) margin <- 0.8
    lmargin <- log(margin)
    ldiff   <- log(theta0)
    se      <- CV2se(CV)
  } else {
    if (missing(theta0)) theta0 <- -0.05
    if (missing(margin)) margin <- -0.2
    lmargin <- margin
    ldiff   <- theta0
    se      <- CV
#    message("It is assumed that CV is the standard deviation.")
  }
  if (length(n) == 1) {
    # total n given    
    # for unbalanced designs we divide the ns by ourself
    # to have only small imbalance (function nvec() from Helper_dp.R)
    n <- nvec(n=n, grps=ades$steps)
    if (n[1] != n[length(n)]) {
      message("Unbalanced design. n(i)=", paste(n, collapse="/"), " assumed.")
    } 
  } else {
    if (length(n) != ades$steps) {
      stop("Length of n vector must be ", ades$steps, "!")
    }
  }
  nc <- sum(1/n)
  n  <- sum(n)
  se.fac <- sqrt(ades$bkni * nc)
  df     <- eval(dfe)
  if (any(df < 1)) stop("n too low. Degrees of freedom <1!", call.=FALSE)

  # check method
  method <- tolower(method)
  method <- match.arg(method)
  
  pow <- .exppower.noninf(alpha, lmargin, ldiff, se.fac, se, dfCV, df, method)
  
  pow
}