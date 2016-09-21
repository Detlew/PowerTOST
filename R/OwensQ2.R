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
  # Communications in Statistics - Theory and Methods, 21:12, 3427-3462, 
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
