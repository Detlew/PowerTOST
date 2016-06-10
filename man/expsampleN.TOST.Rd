\name{expsampleN.TOST}
\alias{expsampleN.TOST}
\title{
Sample size based on 'expected' power
}
\description{
Calculates the sample size based on 'expected' power for a variety of study 
designs used in bioequivalence studies. See known.designs() for the study 
designs covered.
}
\usage{
expsampleN.TOST(alpha = 0.05, targetpower = 0.8, logscale=TRUE, theta0, 
                theta1, theta2, CV, dfCV, design = "2x2", robust=FALSE,
                method=c("exact", "approx"), print = TRUE, details = FALSE, 
                imax=100)
}
\arguments{
  \item{alpha}{
Error probability. Typically set to 0.05.
  }
  \item{targetpower}{
Power to achieve at least. Must be >0 and <1.
Typical values are 0.8 or 0.9.
  }
  \item{logscale}{
Should the data used on log-transformed or on original scale? TRUE or FALSE.\cr
Defaults to TRUE.
  }  
  \item{theta0}{
'True' or assumed bioequivalence ratio or difference. \cr
Maybe missing. Defaults the to 0.95 if \code{logscale=TRUE} or to 0.05 if 
\code{logscale=FALSE}.
  }
  \item{theta1}{
Lower bioequivalence limit as ratio if \code{logscale=TRUE} or as difference.\cr
Can be missing. Defaults then to 0.8 if \code{logscale=TRUE} or to -0.2 if 
\code{logscale=FALSE}.
  }
  \item{theta2}{
Upper bioequivalence limit as ratio if \code{logscale=TRUE} or as difference.\cr
If not given \code{theta2} will be calculated as \code{1/theta1} if \code{logscale=TRUE},
else as \code{-theta1}.
  }
  \item{CV}{
Coefficient of variation as ratio. Only values > 0 ate allowed.\cr
May be given as vector. Then the CVs were pooled as weighted mean with
their df=degrees of freedom as weights.
  }
  \item{dfCV}{
Degrees of freedom for the CVs. Must be a vector of same length as CV.\cr
\code{dfCV=Inf} is allowed but makes only sense if \code{CV/dfCV} are scalars.\cr
The (pooled) \code{dfCV} has to be >4. 
  }
  \item{design}{
Character string describing the study design.\cr
See \code{known.designs()} for designs covered in this package.
  }
  \item{robust}{
Defaults to FALSE. With that value the usual degrees of freedom will be used.\cr
Set to \code{TRUE} will use the degrees of freedom according to the 'robust' evaluation
(aka Senn's basic estimator). These df are calculated as \code{n-seq}.\cr
See \code{known.designs()$df2} for designs covered in this package.
  }
  \item{method}{
Defaults to \code{method="exact"}.\cr
In that case the expected power will be calculated as expected value of the
power with respect to the (prior) distribution of sigma^2 (inverse gamma
distribution).\cr
Set to \code{method="approx"} will calculate the expected power according to
the formula given in the book from Julious or according to Julious/Owen paper.
}
  \item{print}{
If TRUE (default) the function prints its results.\cr
If FALSE only a data.frame with the results will be returned.
  }
  \item{details}{
If TRUE the design characteristics and the steps during
sample size calculations will be shown.\cr Defaults to FALSE.
  }
  \item{imax}{
Maximum number of steps in sample size search.\cr 
Defaults to 100. Adaption only in very rare cases needed.\cr
Never seen a need for adaption up to now.
  }
}
\details{
The sample size is calculated based on iterative evaluation of 'expected' power 
exact as expected value of the power with respect to the (prior) distribution of 
sigma^2 (inverse gamma distribution) or via Julious formulas based on 
non-central t-distribution.\cr
The start value of the sample size search is taken from a large sample approximation.
}
\value{
A data.frame with the input values and the result of the sample size estimation.\cr
The "Sample size" column contains the \bold{total} sample size in case of all design
implemented.
}
\references{
A.P. Grieve\cr 
"Confidence Intervals and Sample Sizes"\cr
Biometrics 47, 1597-1603 December 1991

O'Hagan et al.\cr
"Assurance in Clinical Trial Design"\cr
Pharmaceut. Statist. 2005; 4: 187-201

S.A. Julious, R.J. Owen\cr
"Sample size calculations for clinical studies allowing for 
uncertainty in variance"\cr
Pharmaceutical Statistics (2006), 5, 29-37

S.A. Julious\cr
"Sample sizes for Clinical Trials"\cr
CRC Press, Chapman & Hall 2010

Bertsche et al.\cr
"The predictive distribution of the residual variability in the 
linear-fixed effects model for clinical cross-over trials"\cr
Biometrical Journal 00 (2016) 0, 1-13\cr
Article first published online: 22 MAR 2016

S. Senn\cr
"Cross-over Trials in Clinical Research"\cr
Second edition\cr
Wiley, Chichester 2002
}
\author{
D. Labes
}
\seealso{
 \code{\link{exppower.TOST}, \link{known.designs}, \link{sampleN.TOST}}
}
\examples{
# Classical 2x2 cross-over, target power = 80\%,
# BE limits 80 ... 125\%, assumed true BE ratio = 95%,
# intra-subject CV=30\% estimated with 10 df
# using all the defaults
expsampleN.TOST(CV=0.3, dfCV=10)
# -> gives n=48 with achieved expected power 0.804617
# Compare this to the usual sample size with CV assumed known ('carved in stone')
sampleN.TOST(CV=0.3)

# More then one CV with corresponding degrees of freedom 
# other settings as above
CVs <- c(0.25, 0.3)
dfs <- c(  22, 10)
expsampleN.TOST(CV=CVs, dfCV=dfs)
# -> gives a pooled CV=0.2664927 with df=32
# and a sample size n=34 with achieved expected power 0.812653 exact
# achieved expected power 0.815019 approximate acc. Julious
}