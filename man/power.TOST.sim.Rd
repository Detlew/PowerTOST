\name{power.TOST.sim}
\alias{power.TOST.sim}
\title{
Power of the TOST procedure obtained via simulations
}
\description{
Power is calculated by simulations of studies (PE via it's normal distribution, 
MSE via it's associated chi-squared distribution) and application of the two 
one-sided t-tests. Power is obtained via ratio of studies found BE to 
# of simulated studies. 
}
\usage{
power.TOST.sim(alpha = 0.05, logscale = TRUE, theta1, theta2, theta0, CV, n, 
               design = "2x2", robust = FALSE, setseed = TRUE, nsims = 1e+05)
}
\arguments{
  \item{alpha}{
Type I error probability, significance level. By convention mostly set to 0.05.
}
  \item{logscale}{
Should the data used on log-transformed or on original scale? TRUE or FALSE.\cr
Defaults to TRUE.
}
  \item{theta1}{
Lower bioequivalence limit.\cr
In case of \code{logscale=TRUE} it is given as ratio, otherwise as diff. to 1.\cr
Defaults to 0.8 if \code{logscale=TRUE} or to -0.2 if \code{logscale=FALSE}.
}
  \item{theta2}{
Upper bioequivalence limit.\cr
If not given theta2 will be calculated as \code{1/theta1} if \code{logscale=TRUE} \cr
or as \code{-theta1} if \code{logscale=FALSE}.
}
  \item{theta0}{
'True' or assumed bioequivalence ratio. \cr
In case of \code{logscale=TRUE} it must be given as ratio, \cr
otherwise as difference to 1. See examples.\cr
Defaults to 0.95 if \code{logscale=TRUE} or to 0.05 if \code{logscale=FALSE}
}
  \item{CV}{
Coefficient of variation as ratio.\cr
In case of cross-over studies this is the within-subject CV,
in case of a parallel-group design the CV of the total variability.
}
  \item{n}{
Number of subjects under study.\cr
Is total number if given as scalar, else number of subjects in the (sequence) 
groups. In the latter case the length of n vector has to be equal to the 
number of (sequence) groups.
}
  \item{design}{
Character string describing the study design.\cr
See \code{known.designs()} for designs covered in this package.
}
  \item{robust}{
Defaults to FALSE. With that value the usual degrees of freedom will be used.\cr
Set to \code{TRUE} will use the degrees of freedom according to the 'robust' evaluation
(aka Senn's basic estimator). These df are calculated as \code{n-seq}.
See \code{known.designs()$df2} for designs covered in this package.\cr
Has only effect for higher-order crossover designs.
}
  \item{setseed}{
Simulations are dependent on the starting point of the (pseudo) random number 
generator. To avoid differences in power for different runs a \code{set.seed(1234567)} 
is issued if \code{setseed=TRUE}, the default.\cr
Set this argument to \code{FALSE} to view the variation in power between 
different runs.
}
  \item{nsims}{
Number of studies to simulate. Defaults to 1E5 = 100 000.
}
}
\value{
Value of power according to the input arguments.
}
\author{
D. Labes
}
\note{
This function was intended for internal check of the analytical power
calculation methods. Use of the analytical power calculation methods 
(\code{power.TOST()}) for real problems is recommended.\cr
For sufficient precision nsims > 1E5 (default) may be necessary. 
Be patient if using nsims=1E6. May take some seconds.
}

\seealso{
\code{\link{power.TOST}}, 
}
\examples{
# using the default design 2x2, BE range 0.8 ... 1.25, logscale, theta0=0.95
power.TOST.sim(alpha=0.05, CV=0.3, n=12)
# should give 0.15054, with nsims=1E6 it will be 0.148533
# exact analytical is
power.TOST(alpha=0.05, CV=0.3, n=12)
# should give 0.1484695

# very unusual alpha setting
power.TOST.sim(alpha=0.9, CV=0.3, n=12)
# should give the same (within certain precision) as
power.TOST(alpha=0.95, CV=0.3, n=12)
# or also within certain precision equal to  
power.TOST(alpha=0.95, CV=0.3, n=12, method="mvt")
# SAS Proc Power gives here the incorrect value 0.60525
}
