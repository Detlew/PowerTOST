\name{exppower.noninf}
\alias{exppower.noninf}
\title{
'Expected' power of non-inferiority test
}
\description{
Calculates the 'expected' power for a variety of study
designs used in bioequivalence studies.
}
\usage{
exppower.noninf(alpha = 0.025, logscale=TRUE, theta0, margin, CV, dfCV, n, 
                design = "2x2", robust=FALSE, method=c("exact", "approx"))
}
\arguments{
  \item{alpha}{
Type I error probability, significance level. Defaults here to 0.025.
  }
  \item{logscale}{
Should the data used on log-transformed or on original scale? TRUE or FALSE.\cr
Defaults to TRUE.
  }  
  \item{theta0}{
'True' or assumed bioequivalence ratio or difference. 
Typically set to 0.95 (default if missing) if \code{logscale=TRUE}.\cr
Defaults to -0.05 if \code{logscale=FALSE}.
  }
  \item{margin}{
Non-inferiority margin.\cr
In case of \code{logscale=TRUE} it must be given as ratio, otherwise as diff.\cr
Defaults to 0.8 if \code{logscale=TRUE} or to -0.2 if \code{logscale=FALSE}.
  }
  \item{CV}{
Coefficient of variation as ratio.\cr
Only values > 0 ate allowed.
  }
  \item{dfCV}{
Degrees of freedom for the CV (error/residual degree of freedom). 
\code{dfCV} has to be >4.\cr
\code{dfCV=Inf} is allowed and will give the same result as
\code{power.noninf(...)}.
}
  \item{n}{
Number of subjects under study.\cr
Is total number if given as scalar, else number of subjects in the (sequence) 
groups. In the latter case the length of n vector has to be equal to the 
number of (sequence) groups.
}
  \item{design}{
Character string describing the study design.
See \code{known.designs()} for designs covered in this package.
}
  \item{robust}{
Defaults to FALSE.\cr
Set to \code{TRUE} will use the degrees of freedom according to the 'robust' evaluation
(aka Senn's basic estimator). These df are calculated as \code{n-seq}.\cr
See \code{known.designs()$df2} for designs covered in this package.
}
  \item{method}{
Defaults to \code{method="exact"}.\cr
In that case the expected power will be calculated as expected value of the
power with respect to the (prior) distribution of sigma^2 (inverse gamma
distribution).\cr
Set to \code{method="approx"} the expected power according to the 
approximate formulas (usin non-central t) given in the book from Julious or in 
the Julious/Owen paper will be calculated.
}
}
\details{
This function calculates the so-called 'expected' power taking into account that 
usually the CV is not known but estimated from a previous study / studies with 
an uncertainty. See references.
}
\value{
Value of expected power according to the input.
}
\section{Warning}{
Extreme settings of \code{CV} and/or \code{dfCV} may give incorrect results if
\code{method="exact"} due to the use of \code{integrate()} for calculating the 
expected value of the power with respect to the (prior) distribution of sigma^2. 
The integrand may in such cases become approximately zero over nearly all 
its range. Then \code{integrate()} may fail seriously. See \code{?integrate}.
}
\author{
D. Labes & B. Lang
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
}
\seealso{
 \code{\link{expsampleN.noninf}, \link{power.noninf}, \link{exppower.TOST}} 
}
\examples{
# expected power for non-inferiority test of a 2x2 crossover
# CV 30\% known from a pilot study with 12 subjects (-> dfCV=10)
# using all the defaults for other parameters
# should give: [1] 0.6761068
exppower.noninf(CV=0.3, dfCV=10, n=40)
# Julious approximation
exppower.noninf(CV=0.3, dfCV=10, n=40, method="approx")
# should give: [1] 0.6751358

# Compare this to the usual power (CV assumed known, "carved in stone")
power.noninf(CV=0.3, n=40)
# should give: [1] 0.7228685
# same as 'expected' power with dfCV=Inf
exppower.noninf(CV=0.3, dfCV=Inf, n=40)

# may also calculated via exppower.TOST() after setting upper acceptance limit 
# to Inf and alpha=0.025: settings of first example above
exppower.TOST(CV=0.3, dfCV=10, n=40, theta2=Inf, alpha=0.025)

}
