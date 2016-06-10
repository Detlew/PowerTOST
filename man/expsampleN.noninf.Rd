\name{expsampleN.noninf}
\alias{expsampleN.noninf}
\title{
Sample size based on 'expected' power for the non-inferiority test
}
\description{
Calculates the sample size based on 'expected' power for a
variety of study designs used in bioequivalence studies. \cr
See \code{known.designs()} for the study designs covered.
}
\usage{
expsampleN.noninf(alpha = 0.025, targetpower = 0.8, logscale=TRUE,
                  theta0, margin, CV, dfCV, design = "2x2", robust=FALSE, 
                  method=c("exact", "approx"), print = TRUE, details = FALSE, 
                  imax=100)
}
\arguments{
  \item{alpha}{
Error probability. Typically set to 0.025 for one-sided test.
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
Maybe missing. Defaults then to 0.95 if \code{logscale=TRUE} or to -0.05 if 
\code{logscale=FALSE}.
  }
  \item{margin}{
Non-inferiority margin.\cr
In case of \code{logscale=TRUE} it must be given as ratio, otherwise as diff.\cr
Defaults to 0.8 if \code{logscale=TRUE} or to -0.2 if \code{logscale=FALSE}.
  }
  \item{CV}{
Coefficient of variation as ratio. Only values > 0 ate allowed.\cr
May be given as vector. Then the CVs were pooled as weighted mean (of s2) with
their df (degrees of freedom) as weights.
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
exact as expected value of the power with respect to the (prior) distribution 
of sigma^2 (inverse gamma distribution) or via Julious formulas based on 
non-central t-distribution.\cr
The start value of the sample size search is taken from a large sample approximation.\cr
The sample size is bound to 4 as minimum.
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

S.A. Julious\cr
"Sample sizes for Clinical Trials"\cr
CRC Press, Chapman & Hall, Boca Raton 2010
}
\author{
D. Labes
}
\seealso{
 \code{\link{exppower.noninf}, \link{sampleN.noninf}} 
}
\examples{
# Classical 2x2 cross-over, target power = 80\%, alpha=0.025
# logscale=TRUE, 'non-superiority' margin 125\%, assumed true BE ratio = 105\%,
# intra-subject CV=30\% estimated with 10 df
# using all the other defaults
expsampleN.noninf(theta0=1.05, margin=1.25, CV=0.3, dfCV=10)
# -> gives n=56 with achieved expected power 0.0.806862
# Compare this to the usual sample size with CV known as 'carved in stone'
sampleN.noninf(theta0=1.05, margin=1.25, CV=0.3)

# More then one CV with corresponding degrees of freedom 
# other parameters as above
CVs <- c(0.25, 0.3)
dfs <- c(  22, 10)
expsampleN.noninf(theta0=1.05, margin=1.25, CV=CVs, dfCV=dfs)
# -> gives a pooled CV=0.2664927 with df=32
# and a sample size n=40 with achieved expected power 0.806073 exact
# sample size n=40 with achieved expected power 0.808157 acc. to Julious/Owen
}