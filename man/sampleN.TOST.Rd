\encoding{utf-8}
\name{sampleN.TOST}
\alias{sampleN.TOST}
\title{
Sample size based on power of TOST
}
\description{
Estimates the necessary sample size to have at least a given power.
}
\usage{
sampleN.TOST(alpha = 0.05, targetpower = 0.8, logscale = TRUE, 
             theta0, theta1, theta2, CV, design = "2x2", method="exact", 
             robust=FALSE, print = TRUE, details = FALSE, imax=100)
}
\arguments{
  \item{alpha}{
Significance level (one-sided). Commonly set to 0.05.
  }
  \item{targetpower}{
Power to achieve at least. Must be >0 and <1.\cr
Typical values are 0.8 or 0.9.
  }
  \item{logscale}{
Should the data used on log-transformed (\code{TRUE}) or on original
scale (\code{FALSE})? Defaults to \code{TRUE}.
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
  \item{design}{
Character string describing the study design.\cr
See \code{known.designs()} for designs covered in this package.
  }
  \item{method}{
Method for calculation of the power.\cr  
Defaults to \code{"exact"} in which case the calculation is done based on formulas 
with \enc{Owen’s}{Owen's} Q. The calculation via \enc{Owen’s}{Owen's} Q can also be choosen with 
\code{method="owenq"}.\cr
Another exact method via direct use of the bivariate non-central \emph{t}-distribution
may be chosen with \verb{method="mvt"}. This may have somewhat lower precision
compared to \enc{Owen’s}{Owen's} Q and has a much longer run-time.\cr  
Approximate calculations can be choosen via \verb{method="noncentral"} or 
\verb{method="nct"} for the approximation using the non-central \emph{t}-distribution.
With \verb{method="central"} or \verb{method="shifted"} the relatively crude
approximation via the \sQuote{shifted} central \emph{t}-distribution is chosen.\cr
The strings for \verb{method} may be abbreviated.
  }
  \item{robust}{
Defaults to \code{FALSE}. With that value the usual degrees of freedom will be used.\cr
Set to \code{TRUE} will use the degrees of freedom according to the \sQuote{robust} evaluation
(aka \enc{Senn’s}{Senn's} basic estimator). These df are calculated as \code{n-seq}.\cr
See \code{known.designs()$df2} for designs covered in this package.\cr
Has only effect for higher-order crossover designs.
  }
  \item{print}{
If \code{TRUE} (default) the function prints its results. If \code{FALSE} only the data.frame with the results will be returned. 
  }
  \item{details}{
If \code{TRUE} the design characteristics and the steps during
sample size calculations will be shown. Defaults to \code{FALSE}.
  }
  \item{imax}{
Maximum number of steps in sample size search.\cr 
Defaults to 100. Adaption only in rare cases needed.\cr
  }
}
\details{
The sample size is estimated via iterative evaluation of power of the \acronym{TOST} procedure.\cr
Start value for the sample size search is taken from a large sample approximation
according to Zhang, modified.\cr
The sample size is bound to 4 as minimum.\cr\cr
The estimated sample size gives always the \emph{total} number of subjects (not subject/sequence in crossovers or subjects/group in parallel designs -- like in some other software packages).
}
\value{
A data.frame with the input and results will be returned.\cr
The \code{Sample size} column contains the total sample size.
}
\note{
Of course it is highly recommended to use the default \verb{method="exact"}. :-)\cr
There is no reason besides testing and for comparative purposes to use an 
approximation if the exact method is available at no extra costs.
}
\section{Warning}{
The function does not vectorize properly.\cr
If you need sample sizes with varying CVs, use f.i. for-loops or the apply-family.
}
\references{
Phillips KF. \emph{Power of the Two One-Sided Tests Procedure in Bioequivalence.} J Pharmacokin Biopharm. 1990;18:137--44. \doi{10.1007/BF01063556}

Diletti D, Hauschke D, Steinijans VW. \emph{Sample Size Determination for Bioequivalence Assessment by Means of Confidence Intervals.} Int J Clin Pharmacol Ther Toxicol. 1991;29(1):1--8.\cr

Diletti D, Hauschke D, Steinijans VW. \emph{Sample size determination: Extended tables for the multiplicative model and bioequivalence ranges of 0.9 to 1.11 and 0.7 to 1.43.} Int J Clin Pharmacol Ther Toxicol. 1992;30(Suppl 1):S59--62.

Zhang P. \emph{A Simple Formula for Sample Size Calculation in Equivalence Studies.} J Biopharm Stat. 2003;13(3):529--538. \doi{10.1081/BIP-120022772}
}
\author{
D. Labes
}
\seealso{
  \code{\link{power.TOST}, \link{known.designs}}
}
\examples{
# Exact calculation for a classical 2x2 cross-over (TR/RT),
# BE limits 80 ... 125\%, assumed true BE ratio 0.95, intra-subject CV=30\%,
# using all the default values
# should give n=40 power=0.815845
sampleN.TOST(CV = 0.3)

# Exact calculation for a parallel group design
# evaluation on the original (untransformed) scale
# BE limits 80 ... 120\% = -20\% ... +20\% of reference,
# assumed true BE ratio 0.95\% = -5\% to reference mean,
# total CV=20\%
# should give n=48 (total) power=0.815435
sampleN.TOST(logscale = FALSE, theta1 = -0.2, theta0 = -0.05,
             CV = 0.2, design = "parallel")

# A rather strange setting of theta0! Have a look at n.
# It would be better this is not the sample size but the running total
# of my bank account. But the first million is the hardest. ;-)
sampleN.TOST(CV = 0.2, theta0 = 0.8005, theta1 = 0.8)
}
