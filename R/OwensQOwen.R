#-----------------------------------------------------------------------------
# Calculation of Owens Q-function by the algo given by Owen himself:
# repeated integration by parts
# 
# Author: dlabes Mar 2012
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
# Owen's T-function: 
# Calculates integral 0 to a of exp(-0.5*h^2*(1+x^2))/(1+x^2)/(2*pi)
# most simple implementation via integrate()
# TODO: consider eventually a numerical algo with high precision
# according to M. PATEFIELD, D. TANDY

OT_integrand <- function(x, h) exp(-0.5*h^2*(1+x^2))/(1+x^2)/(2*pi)

#This function has quirks if h is large. see below
OwensT_old <- function(h, a)
{ 
  int <- integrate(OT_integrand, lower=0, upper=abs(a), h=h)$value
  # in case of a=Inf or -Inf the condition T(h,-a)=-T(h,a) is not maintained!
  # thus do it by ourself
  int <- ifelse(a<0, -int, int)
  return(int) 
}
# ----------------------------------------------------------------------------
# first attempt of reworked Owen's T function to handle numeric issues if h is large
# example OwensT_old(0, 1e5) 
# gives an error if stop.on.error=TRUE or gives erroneously
# -1.289739e-06, 
# correct is 0.2499984, nearly the same as
# OwensT_old(0, Inf) = 0.25
# ----------------------------------------------------------------------------
# using formulas (2) found in 
# Patefield M, Tandy D
# "Fast and Accurate Calculation of Owen's T-Function"
# https://www.jstatsoft.org/article/view/v005i05/t.pdf
#
# but using still integrate()
OwensT2 <- function(h, a){
  # T(h, -a) = -T(h, a) thus we takes abs and care later for the sign
  aa <- abs(a)
  # T(-h, a) = T(h, a) thus we can take abs
  hh <- abs(h)
  if (aa<=1){
    #T(h, a)
    int <- integrate(OT_integrand, lower=0, upper=aa, h=hh, rel.tol=1e-5)$value
  } else {
    ah <- aa*hh
    #T(a*h, 1/a)
    int <- integrate(OT_integrand, lower=0, upper=1/aa, h=ah, rel.tol=1e-5)$value
    #equation 2.3, Owen 1956
    phh <- pnorm(hh)
    pah <- pnorm(ah)
    int <- 0.5*(phh + pah) - phh*pah - int
  }
  # take care of sign
  # T(h, -a) = -T(h, a)
  int <- ifelse(a<0, -int, int)
  int
}
# OwensT according AS76 (see below) needs ~25-50% longer run time as OwensT2 
# if h,a are not among the special cases and implementation of Gauss quadraturde
# is via a loop

# ----------------------------------------------------------------------------
# Owen's T-function according to algorithm AS76 and remarks AS R65, AS R80
# R port of the FORTRAN code in the References and matlab code given on
# https://people.sc.fsu.edu/~jburkardt/m_src/asa076/asa076.html
# by J. Burkhardt, license GNU LGPL
#
# no trouble with integrate()
# arguments must be scalars!
OwensT <- function(h, a)
{
  eps <- .Machine$double.eps # 2.220446e-16 on my machine
  if (abs(a)<eps | !is.finite(h) | abs(1-abs(a))<eps | abs(h)<eps | 
      !is.finite(abs(a))) {
    # special cases
    # D.B. Owen
    # "Tables for computing bivariate normal Probabilities"
    # The Annals of Mathematical Statistics, Vol. 27 (4) Dec. 1956, pp. 1075-1090 
    if (abs(a)<eps)        return(0)                                 # a==0
    if (!is.finite(h))     return(0)                                 # h==Inf, -Inf
    if (abs(1-abs(a))<eps) return(sign(a)*0.5*pnorm(h)*(1-pnorm(h))) # a==1
    if (abs(h)<eps)        return(atan(a)/2/pi)                      # h==0
    # a==Inf, -Inf
    if (!is.finite(abs(a))){
      if(h<0) tha <- pnorm(h)/2. else tha <- (1-pnorm(h))/2
      return(sign(a)*tha)
    }
  }
  
  # Boys AS R80
  # AS R89 recommends *not* to use this
  # if (abs(h) < 0.3 && 7.0 < abs(a)){
  #   lam <- abs(a*h)
  #   ex  <- exp(-lam*lam/2.0)
  #   g   <- pnorm(lam)
  #   c1  <- (ex/lam + sqrt(2.0*pi)*(g - 0.5) )/(2.0*pi)
  #   c2  <- ((lam*lam + 2.0)*ex/lam/lam/lam + sqrt(2.0*pi)*(g - 0.5))/(12.0*pi)
  #   ah  <- abs(h)
  #   tha <- abs(0.25 - c1*ah + c2*ah^3)
  #   # next is paranoia?
  #   if (a<0.0) tha <- -tha
  #   return(tha)
  # } 
  aa <- abs(a)
  if (aa <= 1.0){
    # sign is controlled in tfn()
    tha <- tfn(h, a)
    return(tha)
  } else {
    ah  <- aa * h
    gh  <- pnorm(h)
    gah <- pnorm (ah);
    tha <- 0.5*(gh + gah) - gh*gah - tfn(ah, 1.0/aa);
  } 
  if (a < 0.0) tha <- -tha
  return(tha)
  
} #end function

# Auxillary function 
# R port of the matlab code given on
# https://people.sc.fsu.edu/~jburkardt/m_src/owens/owens.html
# by J. Burkhardt license GNU LGPL
# the home page with the MATLAB code is gone meanwhile!
#
# is called as tfn(h, a) if a<=1
# otherwise as tfn(a*h, 1/a)
tfn <- function(x, fx)
{
  # constants
  ng  <- 5
  #r   <- c(0.1477621, 0.1346334, 0.1095432, 0.0747257, 0.0333357)
  #u   <- c(0.0744372, 0.2166977, 0.3397048, 0.4325317, 0.4869533)
  
  # coefficients with more decimals (15) for the 10-point Gauss quadrature 
  # THX to PharmCat who took the coefficients from Julia
  # but does it make sense?
  r <- c(0.14776211235737646, 0.13463335965499826, 0.10954318125799103, 
         0.07472567457529028, 0.033335672154344104)
  u <- c(0.07443716949081561, 0.2166976970646236, 0.3397047841495122, 
         0.43253168334449227, 0.48695326425858587)
  
  tp  <- 1/(2.*pi) # 0.159155
  tv1 <- .Machine$double.eps # 1e-19 in AS R65, 1E-35 in matlab code
                             # 2.220446e-16 on my machine
  tv2 <- 15.0                # is =13 in AS R65
  tv3 <- 15.0
  tv4 <- 1.0E-05
  #
  # Test for X (=h) near zero
  # this is superflous since it is handled in the higher level function
  #if (abs(x) < tv1) return(tp*atan(fx))
  #
  # Test for large values of abs(x) = abs(h)
  # May be this is also superflous?
  if (tv2 < abs(x)) return(0)
  #
  # Test for FX (= a) near zero.
  # this is superflous since it is handled in the higher level function
  #if (abs(fx) < tv1) return(0)
  #
  # Test whether abs(FX) is so large that it must be truncated
  # Is this really necessary? Since we call this function with fx = a <=1
  # or otherwise with 1/a.
  xs  <- -0.5*x*x
  x2  <- fx
  fxs <- fx * fx
  # Computation of truncation point by Newton iteration
  if (tv3 <= log(1.0 + fxs) - xs*fxs){
    x1  <- 0.5 * fx
    fxs <- 0.25 * fxs
    while (1){
      rt  <- fxs + 1.0
      x2  <- x1 + (xs*fxs + tv3 - log(rt))/(2.0*x1*(1.0/rt - xs))
      fxs <- x2 * x2;
      if (abs(x2 - x1) < tv4) break
      x1 <- x2
    }
  }
  #
  # 10 point Gaussian quadrature.
  # original via loop
  # rt <- 0.0
  # for (i in 1:ng) {
  #   r1 <- 1.0 + fxs*(0.5 + u[i])^2;
  #   r2 <- 1.0 + fxs*(0.5 - u[i])^2;
  # 
  #   rt <- rt + r[i]*(exp(xs*r1)/r1 + exp(xs*r2)/r2)
  # }
  # vectorized form, gives a run-time boost of ~50% compared to loop
  # example: system.time(for (i in 1:10000) OwensT(h=2, a=0.7))
  # loop: 0.45 sec, vectorized: 0.22, OwensT2(integrate() solution): 0.34
  r1 <- 1.0 + fxs*(0.5 + u)^2
  r2 <- 1.0 + fxs*(0.5 - u)^2
  rt <- sum(r*(exp(xs*r1)/r1 + exp(xs*r2)/r2))
  return(rt * x2 * tp)
}

# ----------------------------------------------------------------------------
# Owen's Q-function
# Calculates Owen's Q-function via repeated integration by parts
# formulas as given in 
# Owen, D B (1965)
# "A Special Case of a Bivariate Non-central t-Distribution"
# Biometrika Vol. 52, pp.437-446.

OwensQOwen <- function(nu, t, delta, a=0, b)
{
  if (nu<1) stop("nu must be >=1!")
  if (a != 0) stop("Only a=0 implemented!")
  
  # return nct if b is infinite
  if (!is.finite(b)) return(pt(t, df=nu, ncp=delta))
  # use a large delta if delta is infinite, using Inf results in an error
  if (!is.finite(delta)) delta <- sign(delta)*1e20
  
  A   <- t/sqrt(nu)
  B   <- nu/(nu + t*t)
  upr <- nu-2           # upper index of integration by parts
  # the coefficients a(k)
  av  <- vector(mode="numeric", length=nu)
  for (k in seq_along(av)){
    if (k==1 | k==2) av[k] <- 1 else av[k] <- 1/((k-2)*av[k-1])
  }
  ll <- ifelse((upr-1)>0, upr-1, 0)
  L  <- vector(mode="numeric", length=ll)
  # k-1 of the formulas transformes to k here
  if(is.finite(b)){
    # all L[k] are NaN if b==Inf, since dnorm(Inf)==0 and 0*Inf=NaN
    # we use L[k]=0 instead
    for (k in seq_along(L)){
                       # gives NaN 
      if (k==1) L[1] <- 0.5*A*B*b*dnorm(b)*dnorm(A*b-delta)
        else L[k] <- av[k+3]*b*L[k-1]
    }
  } 
  # the series 0 to k is stored as 1 to k+1
  ll <- ifelse((upr+1)>0, upr+1, 0)
  H  <- vector(mode="numeric", length=ll)
  # k+1 of the formulas transformes to k here
  if(is.finite(b)){
    # since H[1]==0 we also get NaN here if b==Inf  
    for (k in seq_along(H)){
      if (k==1) H[1] <- - dnorm(b)*pnorm(A*b-delta)
        else H[k] <- av[k+1]*b*H[k-1]
    }
  }  
  M    <- vector(mode="numeric", length=ll)
  sB   <- sqrt(B)
  # k+1 in the formulas transformes to k here
  for (k in seq_along(M)){
    if (k==1) M[1] <- A*sB*dnorm(delta*sB)*( pnorm(delta*A*sB) - 
                                             pnorm((delta*A*B-b)/sB) )
    if (k==2) M[2] <- B*( delta*A*M[1] + A*dnorm(delta*sB)*
                        ( dnorm(delta*A*sB)-dnorm((delta*A*B-b)/sB)) )
    if (k>2)  M[k] <- ((k-2)/(k-1))*B*(av[k-1]*delta*A*M[k-1] + M[k-2]) - L[k-2]
  }

  sumt <- 0 
  if (2*(nu%/%2)!= nu){
    # odd values of nu
    # sum over odd indices 1, 3, ..., nu-2
    # they are stored in k+1 (0 as index is not allowed), 
    # i.e. at the even indices
    if (upr>=1){
      # for (k in seq(1, upr, by=2)) sumt <- sumt + M[k+1] + H[k+1]
      # vectorized form of a for loop
      k    <- seq(1, upr, by=2)
      sumt <- sum(M[k+1]) + sum(H[k+1])
    }
      # if b==Inf what then?
      # rewrite second argument in first OwensT call (A*b-delta)/b which 
      # gives NaN to A-delta/b which gives A if b==Inf
      qv <- pnorm(b) - 2*OwensT(b, A-delta/b) - 
                       2*OwensT(delta*sB, (delta*A*B-b)/B/delta) +
                       2*OwensT(delta*sB, A) - (delta>=0) + 2*sumt
  } else {
    # even values of nu
    # sum over even indices 0, 2, ..., nu-2
    # they are stored in k+1
    if (upr>=0){
      # vectorized form of the for loop
      k    <- seq(0, upr, by=2)
      sumt <- sum(M[k+1]) + sum(H[k+1])
    }
    qv <- pnorm(-delta) + sqrt(2*pi)*sumt
  }
  return(qv)
}

# ----------------------------------------------------------------------------
# raw power function using OwensQOwen
.power.TOST.Q0 <- function(alpha=0.05, ltheta1, ltheta2, diffm, se, n, df, bk=2)
{
  tval   <- qt(1 - alpha, df, lower.tail = TRUE)
  # if alpha>0.5 (very unusual) then b=R is negative if not using abs()
  # in the application of OwensQ the upper integration limit 
  # is lower then the lower integration limit!
  # SAS OwenQ gives missings if b or a are negative!
  
  delta1 <- (diffm-ltheta1)/(se*sqrt(bk/n))
  delta2 <- (diffm-ltheta2)/(se*sqrt(bk/n))
  # 0/0 -> NaN in case diffm=ltheta1 or diffm=ltheta2 and se=0!
  delta1[is.nan(delta1)] <- 0
  delta2[is.nan(delta2)] <- 0
    
  # R is infinite in case of alpha>=0.5
  R <- (delta1-delta2)*sqrt(df)/(2.*tval)
  R[R<0] <- Inf
  # to avoid numerical errors in OwensQ implementation
  if (min(df)>10000) { 
    # Joulious formula (57) or (67), normal approximation
    p1 <- pnorm( (abs(delta1)-tval), lower.tail = TRUE)
    p2 <- pnorm( (abs(delta2)-tval), lower.tail = TRUE)
    
    return(p1 + p2 -1.)
  }
  # attempt to vectorize (it vectorizes properly if diffm is a vector
  # OR se OR n,df are vectors) 
  nel <- length(delta1)
  dl <- length(tval)
  p1 <- c(1:nel)	
  p2 <- p1
#  print(df)
#  print(tval)
#  print(delta1)
#  print(delta2)
#  print(R)
  for (i in seq_along(delta1)) {
    if (dl>1) {ddf <- df[i]; ttt <- tval[i]} 
        else {ddf <- df[1]; ttt <- tval[1]}
    p1[i] <- OwensQOwen(ddf,  ttt, delta1[i], 0, R[i])
    p2[i] <- OwensQOwen(ddf, -ttt, delta2[i], 0, R[i])
  }
  return( p2-p1 )
}

# check the function(s)
#require(PowerTOST)
#se <- CV2se(0.075)
#.power.TOST.Q0(alpha=0.05, ltheta1=log(0.8), ltheta2=log(1.25), 
#               diffm=0, se=se, n=4, df=2, bk=2)
#PowerTOST:::.power.TOST(alpha=0.05, ltheta1=log(0.8), ltheta2=log(1.25), 
#                        diffm=0, se=se, n=4, df=2, bk=2)
#.power.TOST.Q0(alpha=0.05, ltheta1=log(0.8), ltheta2=log(1.25), 
#               diffm=0, se=se, n=6, df=4, bk=2)
#PowerTOST:::.power.TOST(alpha=0.05, ltheta1=log(0.8), ltheta2=log(1.25), 
#                        diffm=0, se=se, n=6, df=4, bk=2)
#           
#
#OwensQ(2,-2.919986,-4.213542,0,2.040712)
#OwensQ(2,2.919986,4.213542,0,2.040712)
#
#OwensQOwen(2, 2.919986, 4.213542,0,2.040712)
