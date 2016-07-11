\name{sampleN.2TOST}
\alias{sampleN.2TOST}
\title{
Sample size based on power of 2 TOSTs
}
\description{
Calculates the necessary sample size to have at least a given power when two
parameters are being tested simultaneously.
}
\usage{
sampleN.2TOST(alpha = c(0.05, 0.05), targetpower = 0.8, logscale = TRUE, 
              theta0, theta1, theta2, CV, rho, design = "2x2", setseed = TRUE,
              robust = FALSE, print = TRUE, details = FALSE, imax = 100)
}
\arguments{
  \item{alpha}{
  Vector; contains one-sided significance level for each of the two TOSTs.\cr
  For one TOST, by convention mostly set to 0.05.
  }
  \item{targetpower}{
Power to achieve at least. Must be >0 and <1.\cr
Typical values are 0.8 or 0.9.
  }
  \item{logscale}{
Should the data used on log-transformed or on original scale? TRUE or FALSE.\cr
Defaults to TRUE.
  }
  \item{theta0}{
Vector; contains 'true' assumed bioequivalence ratio for each of the two TOSTs. \cr
In case of \code{logscale=TRUE} each element must be given as ratio, \cr
otherwise as difference to 1. See examples.\cr
Defaults to \code{c(0.95, 0.95)} if \code{logscale=TRUE} or to 
\code{c(0.05, 0.05)} if \code{logscale=FALSE}.
  }
  \item{theta1}{
Vector; contains lower bioequivalence limit for each of the two TOSTs.\cr
In case of \code{logscale=TRUE} it is given as ratio, otherwise as diff. to 1.\cr
Defaults to \code{c(0.8, 0.8)} if \code{logscale=TRUE} or to \code{c(-0.2, -0.2)}
if \code{logscale=FALSE}.
  }
  \item{theta2}{
Vector; contains upper bioequivalence limit for each of the two TOSTS.\cr
If not given theta2 will be calculated as \code{1/theta1} if \code{logscale=TRUE}\cr
or as \code{-theta1} if \code{logscale=FALSE}.
  }
  \item{CV}{
Vector of coefficient of variations (given as as ratio, e.g. 0.2 for 20\%).\cr
In case of cross-over studies this is the within-subject CV, \cr
in case of a parallel-group design the CV of the total variability.\cr
In case of \code{logscale=FALSE} CV is assumed to be the respective standard 
deviation.
  }
  \item{rho}{
Correlation between the two parameters under consideration. This is defined 
as correlation between the estimator of the treatment difference of
parameter one and the estimator of the treatment difference of parameter two.
}
  \item{design}{
Character string describing the study design.\cr
See \code{known.designs()} for designs covered in this package.
  }
  \item{setseed}{
Calculation depends on \code{pmvt()} which is based on randomized quasi Monte Carlo
methods. If \code{setseed=TRUE} a seed value is set, the default.
}
  \item{robust}{
Defaults to \code{FALSE}. With that value the usual degrees of freedom will be used.\cr
Set to \code{TRUE} will use the degrees of freedom according to the 'robust' evaluation
(aka Senn's basic estimator). These df are calculated as \code{n-seq}.\cr
See \code{known.designs()$df2} for designs covered in this package.\cr
Has only effect for higher-order crossover designs.
  }
  \item{print}{
If \code{TRUE} (default) the function prints its results.\cr 
If \code{FALSE} only the result list will be returned. 
  }
  \item{details}{
If \code{TRUE} the design characteristics and the steps during
sample size calculations will be shown. \cr 
Defaults to \code{FALSE}.
  }
  \item{imax}{
Maximum number of steps in sample size search. \cr 
Defaults to 100.
  }
}
\details{
The sample size is calculated via iterative evaluation of power of the 2 TOSTs.\cr
Start value for the sample size search is taken from a large sample approximation
(1 TOST) according to Zhang, modified.\cr
The sample size is bound to 4 as minimum.
}
\value{
A list with the input and results will be returned.\cr
The element name "Sample size" contains the total sample size.
}
\note{
If both theta0 are near the acceptance limits then the starting value may not
be a good approximation resulting in a lot of iteration steps; imax may need
to be increased to obtain the required sample size.\cr\cr
}
\section{Warning}{
The function does not vectorize properly.\cr
If you need sample sizes with varying CVs f.i. use for-loops or the apply-family.
}
\references{
Hua S. Y., Xu S., and D'Agostino Sr. R. B.\cr
"Multiplicity adjustments in testing for bioequivalence"\cr
Statistics in Medicine, Vol. 34, Issue 2, 215-231 (2015)

Lang B., Fleischer F.\cr
"Letter to the Editor: Comments on 'Multiplicity adjustments in testing for bioequivalence'".\cr
Statistics in Medicine, Vol. 35, Issue 14, 2479-2480 (2016)

Zhang P.\cr
"A Simple Formula for Sample Size Calculation in Equivalence Studies"\cr
J. Biopharm. Stat. 13(3), 529-538 (2003)
}
\author{
Benjamin Lang
}
\seealso{
  \code{\link{power.2TOST}, \link{known.designs}}
}
\examples{
# Sample size for 2x2x2 cross-over design, intra-subject CV = 30\% and assumed
# ratios of 0.95 for both parameters, and correlation 0.9 between parameters
# (using all the other default values)
# Should give n=44 with power=0.808840
sampleN.2TOST(theta0 = rep(0.95, 2), CV = rep(0.3, 2), rho = 0.9)

# Sample size for a parallel group design,
# evaluation on the original (untransformed) scale
# BE limits 80 ... 120\% = -20\% ... +20\% of reference,
# assumed true BE ratio 0.95\% = -5\% to reference mean for both parameters,
# total CV=20\% for both parameters, and correlation 0.9 between parameters
# should give n=52 with power=0.801250
sampleN.2TOST(logscale=FALSE, theta0 = rep(-0.05, 2), CV = c(0.2, 0.2), 
              rho = 0.9, design = "parallel")
}
