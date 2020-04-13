\encoding{utf-8}
\name{power.TOST.sim}
\alias{power.TOST.sim}
\title{
Power of the TOST procedure obtained via simulations
}
\description{
Power is calculated by simulations of studies (PE via its normal distribution, 
MSE via its associated \ifelse{html}{\out{&chi;<sup>2</sup>}}{\eqn{\chi^2}{chi-squared}} distribution) and application of the two one-sided \emph{t}-tests. Power is obtained via ratio of studies found BE to 
the number of simulated studies. 
}
\usage{
power.TOST.sim(alpha = 0.05, logscale = TRUE, theta1, theta2, theta0, CV, n, 
               design = "2x2", robust = FALSE, setseed = TRUE, nsims = 1e+05)
}
\arguments{
  \item{alpha}{
Significance level (one-sided). Commonly set to 0.05.
}
  \item{logscale}{
Should the data used on log-transformed or on original scale? \code{TRUE} (default) or \code{FALSE}.
}
  \item{theta0}{
\sQuote{True} or assumed T/R ratio or difference.\cr
In case of \code{logscale=TRUE} it must be given as ratio T/R.\cr
If \code{logscale=FALSE}, the difference in means. In this case, the difference may be expressed in two ways: relative to the same (underlying) reference mean, i.e. as (T-R)/R = T/R - 1; or as difference in means T-R. Note that in the former case the units of \code{CV}, \code{theta1} and \code{theta2} need also be given relative to the reference mean (specified as ratio).\cr
Defaults to 0.95 if \code{logscale=TRUE} or to 0.05 if \code{logscale=FALSE}
}
  \item{theta1}{
Lower (bio-)equivalence limit.\cr
In case of \code{logscale=TRUE} it is given as ratio.\cr
If \code{logscale=FALSE}, the limit may be expressed in two ways:
difference of means relative to the same (underlying) reference mean or in units of the difference of means.
Note that in the former case the units of \code{CV}, \code{theta0} and \code{theta2} need also be given relative to the reference mean (specified as ratio).\cr
Defaults to 0.8 if \code{logscale=TRUE} or to -0.2 if \code{logscale=FALSE}.
}
  \item{theta2}{
Upper (bio-)equivalence limit.\cr
In case of \code{logscale=TRUE} it is given as ratio.
If \code{logscale=FALSE}, the limit may be expressed in two ways:
difference of means relative to the same (underlying) reference mean or in units of the difference of means.
Note that in the former case the units of \code{CV}, \code{theta0} and \code{theta1} need also be given relative to the reference mean (specified as ratio).\cr
If not given, \code{theta2} will be calculated as \code{1/theta1} if \code{logscale=TRUE} or as \code{-theta1} if \code{logscale=FALSE}.
}
  \item{CV}{
In case of \code{logscale=TRUE} the (geometric) coefficient of variation given as ratio.\cr
If \code{logscale=FALSE} the argument refers to (residual) standard deviation of the response. In this case, standard deviation may be expressed two ways: relative to a reference mean (specified as ratio sigma/muR), i.e. again as a coefficient of variation; or untransformed, i.e. as standard deviation of the response. Note that in the former case the units of \code{theta0}, \code{theta1} and \code{theta2} need also be given relative to the reference mean (specified as ratio).\cr

In case of cross-over studies this is the within-subject CV, in case of a parallel-group design the CV of the total variability.
}
  \item{n}{
Number of subjects under study.\cr
Is total number if given as scalar, else number of subjects in the (sequence) 
groups. In the latter case the length of \code{n} vector has to be equal to the 
number of (sequence) groups.
}
  \item{design}{
Character string describing the study design.\cr
See \code{known.designs()} for designs covered in this package.
}
  \item{robust}{
Defaults to \code{FALSE}. With that value the usual degrees of freedom will be used.\cr
Set to \code{TRUE} will use the degrees of freedom according to the \sQuote{robust} evaluation
(aka \enc{Senn’s}{Senn's} basic estimator). These degrees of freedom are calculated as \code{n-seq}.
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
Number of studies to simulate. Defaults to 100,000 = 1E5.
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
power.TOST.sim(alpha = 0.05, CV = 0.3, n = 12)
# should give 0.15054, with nsims=1E6 it will be 0.148533
# exact analytical is
power.TOST(alpha = 0.05, CV = 0.3, n = 12)
# should give 0.1484695

# very unusual alpha setting
power.TOST.sim(alpha = 0.9, CV = 0.3, n = 12)
# should give the same (within certain precision) as
power.TOST(alpha = 0.95, CV = 0.3, n = 12)
# or also within certain precision equal to  
power.TOST(alpha = 0.95, CV = 0.3, n = 12, method = "mvt")
# SAS Proc Power gives here the incorrect value 0.60525
}
