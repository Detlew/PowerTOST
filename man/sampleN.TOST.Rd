\name{sampleN.TOST}
\alias{sampleN.TOST}
\title{
Sample size based on power of TOST
}
\description{
Calculates the necessary sample size to have at least a given power.
}
\usage{
sampleN.TOST(alpha = 0.05, targetpower = 0.8, logscale = TRUE, 
             theta0, theta1, theta2, CV, design = "2x2", method="exact", 
             robust=FALSE, print = TRUE, details = FALSE, imax=100)
}
\arguments{
  \item{alpha}{
Type I error probability. Per convention mostly set to 0.05.
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
'True' or assumed bioequivalence ratio. \cr
In case of \code{logscale=TRUE} it must be given as ratio, otherwise as 
difference to 1. See examples.\cr
Defaults to 0.95 if \code{logscale=TRUE} or to 0.05 if \code{logscale=FALSE}.
  }
  \item{theta1}{
Lower bioequivalence limit.\cr
In case of \code{logscale=TRUE} it is given as ratio, otherwise as diff. to 1.\cr
Defaults to 0.8 if \code{logscale=TRUE} or to -0.2 if \code{logscale=FALSE}.
  }
  \item{theta2}{
Upper bioequivalence limit.\cr
If not given theta2 will be calculated as \code{1/theta1} if \code{logscale=TRUE}\cr
or as \code{-theta1} if \code{logscale=FALSE}.
  }
  \item{CV}{
Coefficient of variation as ratio.
  }
  \item{design}{
Character string describing the study design.\cr
See \code{known.designs()} for designs covered in this package.
  }
  \item{method}{
Method for calculation of the power.\cr  
Defaults to "exact" in which case the calculation is done based on formulas 
with Owen's Q. The calculation via Owen's Q can also be choosen with 
\code{method="owenq"}.\cr
Another exact method via direct use of the bivariate non-central t-distribution
may be chosen with \code{method="mvt"}. This may have somewhat lower precision
compared to Owen's Q and has much longer run-time.\cr  
Approximate calculations can be choosen via \code{method="noncentral"} or 
\code{method="nct"} for the approximation using the non-central t-distribution.
With \code{method="central"} or \code{method="shifted"} the relative crude
approximation via 'shifted' central t-distribution is chosen.\cr
The strings for \code{method} may be abbreviated.
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
If \code{FALSE} only the data.frame with the results will be returned. 
  }
  \item{details}{
If \code{TRUE} the design characteristics and the steps during
sample size calculations will be shown. \cr 
Defaults to \code{FALSE}.
  }
  \item{imax}{
Maximum number of steps in sample size search. \cr 
Defaults to 100. Adaption only in rare cases needed.
  }
}
\details{
The sample size is calculated via iterative evaluation of power of the TOST procedure.\cr
Start value for the sample size search is taken from a large sample approximation
according to Zhang, modified.\cr
The sample size is bound to 4 as minimum.
}
\value{
A data.frame with the input and results will be returned.\cr
The "Sample size" column contains the total sample size.
}
\note{
Of course it is highly recommended to use the default \code{method="exact"} :-)). \cr
There is no reason beside testing and comparative purposes to use an 
approximation if the exact method is available at no extra costs.
}
\section{Warning }{
The function does not vectorize properly.\cr
If you need sample sizes with varying CVs f.i. use for-loops or the apply-family.
}
\references{
Phillips, K. F. (1990)\cr 
"Power of the Two One-Sided Tests Procedure in Bioequivalence"\cr
J. of Pharmacokinetics and Biopharmaceutics, 18, 137-144.

Diletti, D., Hauschke, D., and Steinijans, V. W. (1991)\cr
"Sample Size Determination for Bioequivalence Assessment by Means of Confidence 
Intervals"\cr
Int. J. of Clinical Pharmacology, Therapy and Toxicology, 29(1), 1-8 (1991)\cr 
30 Suppl.No.1, S51-58 (1992)

Diletti, D., Hauschke, D., and Steinijans, V. W. (1992)\cr
"Sample size determination : Extended tables for the multiplicative model\cr
and bioequivalence ranges of 0.9 to 1.11 and 0.7 to 1.43"\cr
Int. J. of Clinical Pharmacology, Therapy and Toxicology, 30 Suppl.No.1, S59-62

Zhang P.\cr
"A Simple Formula for Sample Size Calculation in Equivalence Studies"\cr
J. Biopharm. Stat. 13(3), 529-538 (2003)
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
sampleN.TOST(CV=0.3)

# Exact calculation for a parallel group design
# evaluation on the original (untransformed) scale
# BE limits 80 ... 120\% = -20\% ... +20\% of reference,
# assumed true BE ratio 0.95\% = -5\% to reference mean,
# total CV=20\%
# should give n=48 (total) power=0.815435
sampleN.TOST(logscale=FALSE, theta1=-0.2, theta0=-0.05, CV=0.2, design="parallel")

# A rather strange setting of theta0! Have a look at n.
# It would be better this is not the sample size but the running total
# of my bank account. But the first million is the hardest ;-).
sampleN.TOST(CV=0.2, theta0=0.8005, theta1=0.8)
}
