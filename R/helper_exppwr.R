#------------------------------------------------------------------------
# helper functions for expected power
#------------------------------------------------------------------------
# Densitiy of inverse gamma distriubtion
# Adapted dinvgamma() from R package MCMCpack (1.3-3)
# Author B. Lang, slightly modified by D. Labes
#
# MCMCpack:dinvgamma() authored by Andrew D. Martin, Kevin M. Quinn, 
# and Jong Hee Park
#------------------------------------------------------------------------
my_dinvgamma <- function(x, shape, scale = 1) {
  stopifnot(is.numeric(shape), shape > 0, is.numeric(scale), scale > 0,
            is.numeric(x), x >= 0)  # include domain for x
  d <- exp(shape * log(scale) - lgamma(shape) - (shape + 1) * log(x) - (scale/x))
  # Support is x > 0, but we return the asymptotic d=0 if x==0
  d[x == 0] <- 0
  d
}

#------------------------------------------------------------------------
# Density of non-standardized Student's t-distribution
# See https://en.wikipedia.org/wiki/Location-scale_family
#------------------------------------------------------------------------
dt_ls <- function(x, df, mu, sigma) {
  stopifnot(is.numeric(sigma), sigma > 0)
  
  dt((x - mu)/sigma, df)/sigma
}

#------------------------------------------------------------------------
# Density of normal inverse gamma distribution
# https://en.wikipedia.org/wiki/Normal-inverse-gamma_distribution
# Author B. Lang
#------------------------------------------------------------------------
dninvgamma <- function(x, v, mu, lambda, alpha, beta) {
  # alpha = shape, beta = scale
  stopifnot(is.numeric(lambda), lambda > 0, is.numeric(alpha), alpha > 0,
            is.numeric(beta), beta > 0,
            is.numeric(x), is.numeric(v), v >= 0)
  d <- exp(log(sqrt(lambda)) - log(sqrt(v)) - log(sqrt(2*pi)) + 
           alpha * log(beta) - lgamma(alpha) - (alpha + 1) * log(v) -
           (2*beta + lambda*(x - mu)^2) / (2*v))
  # Support wrt v is v > 0, but we return the asymptotic d=0 if v==0
  d[v == 0] <- 0
  d
}

# DL: according to
# http://www.cs.utah.edu/~fletcher/cs6190/lectures/GaussianPriors.pdf
dnig <- function(x, v, mu, lambda, alpha, beta) 
{
  dnorm (x, mean = mu, sd =sqrt(v/lambda)) * my_dinvgamma (v, alpha, beta)
} 

#------------------------------------------------------------------------
# Function to calculate degrees of freedom and/or 
# the prefactor of the standard error of the mean difference
# Input: Either sample size and design or degrees of freedom and design
# Author B. Lang
#------------------------------------------------------------------------
get_df_sefac <- function(n = NULL, v = NULL, design, robust = FALSE) 
{
  # Check if design is implemented
  d.no <- .design.no(design)
  if (is.na(d.no))
    stop("Design ",design, " unknown!", call. = FALSE)
  ades <- .design.props(d.no)
  if (!is.null(n)) {
    # Input is sample size and design
    if (!is.null(v))
      warning("Both n and v are given. Input v was ignored.", call. = FALSE)
    if (length(n) == 1) {
      # total n given    
      # for unbalanced designs we divide the ns by ourself
      # to have only small imbalance (function nvec() from Helper_dp.R)
      n <- nvec(n = n, grps = ades$steps)
      if (n[1] != n[length(n)]) {
        message("Unbalanced design. n(i)=", paste(n, collapse = "/"), " assumed.")
      }
    } else {
      if (length(n) != ades$steps) {
        stop("Length of n vector must be ", ades$steps, "!")
      }
    }
    nc <- sum(1/n)
    dfe  <- .design.df(ades, robust = robust)
    n <- sum(n)
    return(list(df = eval(dfe), sefac = sqrt(ades$bkni * nc)))
  }
  if (!is.null(v)) {
    # Input is degrees of freedom and design
    # This feature is currently not used!
    # Note:
    #  - Only gives exact results in case of no missing data
    #  - For all other cases this is only an approximation (which gets worse
    #    the more data are missing)
    if (v <= 4)
      stop("v has to be >4", call. = FALSE)
    f <- function(n) {
      eval(dfe[[1]]) - v
    }
    ssize <- uniroot(f, lower = 4, upper = 100, extendInt = "upX")$root
    # Check for ssize being an integer, see ?is.integer
    tol <- sqrt(.Machine$double.eps)
    if (!(abs(ssize - round(ssize)) < tol)) {
      ssize <- max(ceiling(ssize), ades$steps)
      if (ssize %% ades$steps != 0) 
        ssize <- ades$steps * floor(ssize / ades$steps) + ades$steps
    }
    # Here, df = v is not needed as return value
    # However, for consistency we return the same list structure as above 
    return(list(df = v, sefac = sqrt(ades$bk / ssize)))
  }
}
