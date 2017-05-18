\name{sampleN.noninf}
\alias{sampleN.noninf}
\title{
Sample size for the non-inferiority t-test
}
\description{
Function for calculating the sample size needed to have a pre-specified power
for the one-sided non-inferiority t-test for normal or log-normal distributed data.
}
\usage{
sampleN.noninf(alpha = 0.025, targetpower = 0.8, logscale = TRUE, margin, 
               theta0, CV, design = "2x2", robust = FALSE, 
               details = FALSE, print = TRUE, imax=100)
}
\arguments{
  \item{alpha}{
Type I error probability, significance level. Defaults here to 0.025.
}
  \item{targetpower}{
Power to achieve at least. Must be >0 and <1.\cr
Typical values are 0.8 or 0.9.
}
  \item{logscale}{
Should the data used on log-transformed or on original scale? TRUE or FALSE.\cr
Defaults to TRUE.
}
  \item{margin}{
Non-inferiority margin.
In case of \code{logscale=TRUE} it must be given as ratio, otherwise as diff. to 1.\cr
Defaults to 0.8 if \code{logscale=TRUE} or to -0.2 if \code{logscale=FALSE}.
}
  \item{theta0}{
'True' or assumed bioequivalence ratio or difference. \cr
In case of \code{logscale=TRUE} it must be given as ratio, \cr
otherwise as difference to 1. See examples.\cr
Defaults to 0.95 if \code{logscale=TRUE} or to 0.05 if \code{logscale=FALSE}
}
  \item{CV}{
Coefficient of variation as ratio.\cr
In case of cross-over studies this is the within-subject CV, \cr
in case of a parallel-group design the CV of the total variability.
}
  \item{design}{
Character string describing the study design.\cr
See \code{\link{known.designs}} for designs covered in this package.
}
  \item{robust}{
Defaults to FALSE. With that value the usual degrees of freedom will be used.\cr
Set to \code{TRUE} will use the degrees of freedom according to the 'robust' evaluation
(aka Senn's basic estimator). These df are calculated as \code{n-seq}.\cr
See \code{known.designs()$df2} for designs covered in this package.\cr
Has only effect for higher-order crossover designs.
}
  \item{details}{
If TRUE the design characteristics and the steps during
sample size calculations will be shown. \cr 
Defaults to FALSE.
}
  \item{print}{
If TRUE (default) the function prints its results.\cr 
If FALSE only the data.frame with the results will be returned. 
  }
  \item{imax}{
Maximum number of steps in sample size search. \cr 
Defaults to 100. Adaption only in rare cases needed.
  }
  
}
\details{
The sample size is calculated via iterative evaluation of \code{power.noninf()}.\cr
Start value for the sample size search is taken from a large sample approximation.\cr
The sample size is bound to 4 as minimum.\cr\cr

\strong{Notes on the underlying hypotheses}\cr
If the supplied margin is < 0 (\code{logscale=FALSE}) or < 1 (\code{logscale=TRUE}), 
then it is assumed higher response values are better. The hypotheses are\cr
\code{H0: theta0 <= margin vs. H1: theta0 > margin}\cr
where \code{theta0 = mean(test)-mean(reference)} if \code{logscale=FALSE}\cr
or\cr
\code{H0: log(theta0) <= log(margin) vs. H1: log(theta0) > log(margin)}\cr
where \code{theta0 = mean(test)/mean(reference)} if \code{logscale=TRUE}.\cr\cr

If the supplied margin is > 0 (\code{logscale=FALSE}) or > 1 (\code{logscale=TRUE}), 
then it is assumed lower response values are better. The hypotheses are\cr
\code{H0: theta0 >= margin vs. H1: theta0 < margin}\cr
where \code{theta0 = mean(test)-mean(reference)} if \code{logscale=FALSE}\cr
or\cr
\code{H0: log(theta0) >= log(margin) vs. H1: log(theta0) < log(margin)}\cr
where \code{theta0 = mean(test)/mean(reference)} if \code{logscale=TRUE}.\cr
This latter case may also be considered as 'non-superiority'.
}
\value{
A data.frame with the input settings and results will be returned.\cr
Explore it with \code{str(sampleN.noninf(...)}
}
\section{Warning }{
The function does not vectorize properly.\cr
If you need sample sizes with varying CVs f.i. use for-loops or the apply-family.
}
\references{
Julious SA.\cr
\emph{Sample sizes for clinical trials with Normal data}\cr
Stat Med. 2004;23(12):1921--86. \doi{10.1002/sim.1783}
}
\author{
D. Labes
}
\seealso{
  \code{\link{known.designs}}, \code{\link{power.noninf}}
}
\examples{
# using all the defaults: margin=0.8, theta0=0.95, alpha=0.025
# log-transformed, design="2x2"
sampleN.noninf(CV=0.3)
# should give n=48
#
# 'non-superiority' case, log-transformed data
# with assumed 'true' ratio somewhat above 1
sampleN.noninf(CV=0.3, targetpower=0.9, margin=1.25, theta0=1.05)
# should give n=62
}
