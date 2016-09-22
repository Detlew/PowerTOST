#------------------------------------------------------------------------------
# Owen's Q-function 
# author: dlabes
#-------------------------------------------------------------------------------
# a, b must be a scalar numeric
# nu, t and delta also, no vectors allowed
# This is a variant with variables transformation
# and using the fact that Integral(0,b)+Integral(b,Inf) gives p(nct)
OwensQ <- function (nu, t, delta, a=0, b)
{
  if (a!=0) stop("Only a==0 implemented.")
  
  if (length(nu)>1 | length(t)>1 | length(delta)>1 | length(a)>1 | length(b)>1) 
    stop("Input must be scalars!")
  if (nu<1) stop("nu must be >=1.")
  
  if(a==b) return(0)

  # in case of alpha>=0.5 b is infinite
  # also in case of se=0 and diffm != ltheta1 or !=ltheta2
  # for that case see:
  # A bivariate noncentral T-distibution with applications
  # Youn Min Chou
  # Communications in Statistics - Theory and Methods, 21:12, 3427-3462, 1992
  # DOI: 10.1080/03610929208830988
  # There are sometimes warnings regarding precision of nct, suppress them?
  if(a==0 && is.infinite(b)) return(suppressWarnings(pt(t, df=nu, ncp=delta)))
  # should also work for sufficient high b, but what is sufficient?
  if(a==0 && b>150) return(suppressWarnings(pt(t, df=nu, ncp=delta)))
  
  # we calculate Owen's Q via
  # pt(t, df=nu, ncp=delta) - Integral(b,Inf)
  # the Integral(b,Inf) is via transformation of the variables x=b+y/(1-y)
  # remapped to an integral(0,1)
  # maybe this gives us more precision and numerical stability
  i_fun <- function(y){
    .Q.integrand(b+y/(1-y), nu, t, delta)/(1-y)^2
  }
  Integral01 <- integrate(i_fun, lower=0, upper=1, subdivisions = 1000, 
                          rel.tol = 1.e-8, stop.on.error = TRUE)
  
  # suppress the warning w.r.t. precision of nct
  OQ <- suppressWarnings(pt(t, df=nu, ncp=delta)) - Integral01[[1]]
  
  OQ
                            
}
#-------------------------------------------------------------------------------
# Integrand of the definit integral in Owen's Q. Used in the call of integrate()
# Not useful alone, I think ? Leading . hides this function 
# function must give a vectorized answer in respect to x
.Q.integrand <- function(x, nu, t, delta)
{ #version without for - loop, it works without
  lnQconst <- -((nu/2.0)-1.0)*log(2.0) - lgamma(nu/2.0)
  
  # what if x<0? Should here not possible, but ...
  # simple x^(nu-1) doesnt work for high nu because  = Inf 
  # and then exp( -0.5*x^2 + lnQconst )*x^(nu-1) -> NaN
  # (nu-1)*log(abs(x)) is NaN if nu=1, x=0! 0*(-Inf) -> NaN
  
  dens <- x  # assures that dens=0 if x=0
  x1   <- x[x!=0]
  dens[x!=0] <- sign(x1)^(nu-1) *
    pnorm( t*x1/sqrt(nu) - delta, mean = 0, sd = 1, log.p = FALSE) * 
    exp( (nu-1)*log(abs(x1)) - 0.5*x1^2 + lnQconst )
  # return    
  dens
}

# Test cases:
# Craig Zupke's observations:
# power.TOST(0.410,FALSE,-5.97,5.97,8.5448,1,14,"parallel",TRUE) #!old call
# power.TOST(0.410,FALSE,-5.97,5.97,8.5448,1,14,"parallel","exact")
# gave an error; high b/delta
# should give: 2.335633e-07

# Jul 2012: Helmuts observation
# n=4, CV=1E-5(=se) gives power=1 (delta1=24303.3, delta2=-38811.23, R=b=15283.88
#      CV=1E-6 gives power=0      (      243033          -388112.3   R  152838.8
#      CV=0    gives power=1             Inf              -Inf       Inf
# tval=2.919986
# for CV=1e-6: erroneous in versions pre 0.9-9. The 2. call gave =0
# OwensQ(nu=2, t= 2.919986, delta= 243033,   0, 152838.8) ==0
# OwensQ(nu=2, t=-2.919986, delta=-388112.3, 0, 152838.8) ==1
# for CV=0 - no longer staistics
# OwensQ(nu=2, t=2.919986, delta=Inf, 0, Inf)  ==0
# OwensQ(nu=2, t=-2.919986, delta=-Inf, 0,Inf) ==1
# 
# Helmuts cases (ver 0.9-9) Jul 2012
# sampleN.TOST(theta0=1, CV=0.02, design="2x2", print=TRUE) # Ok
# next gave an error due to 0*-Inf in .Q.integrand()
# sampleN.TOST(theta0=1, CV=0.01, design="2x2", print=TRUE) 
#
# Jiri Hofmann's case(s) (V1.1-08, Jan 2014)
# power.TOST(CV=0.39, n=7000, theta0=1.24) gave power=0.3511889, correct
# power.TOST(CV=0.385, n=7000, theta0=1.24) gave power=0, correct is 0.3568784
# power.TOST(CV=0.38, n=7000, theta0=1.24) gave power=0, correct is 0.3627554
# power.TOST(CV=0.375, n=7000, theta0=1.24) gave power=0, correct is 0.3688277
# power.TOST(CV=0.37, n=7000, theta0=1.24) gave power=0.3751032, correct
#
# power.TOST(CV=0.37, n=7000, theta0=1.0) gave power=0, correct is 1
# power.TOST(CV=0.375, n=7000, theta0=1.0) gave power=0, correct is 1
# power.TOST(CV=0.38, n=7000, theta0=1.0) gave power=0, correct is 1
# power.TOST(CV=0.385, n=7000, theta0=1.0) gave power=0, correct is 1
# power.TOST(CV=0.39, n=7000, theta0=1.0) gave power=1, correct
