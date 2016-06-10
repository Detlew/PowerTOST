\name{CVCL}
\alias{CVCL}
\title{
Confidence limits of a CV for log-normal data
}
\description{
The function calculates the 1-alpha confidence limits (either 1-sided
or 2-sided) via the chi-squared distribution of the error variance the CV 
is based on.
}
\usage{
CVCL(CV, df, side = c("upper", "lower", "2-sided"), alpha = 0.05)
}
\arguments{
  \item{CV}{Coefficient of variation}
  \item{df}{degrees of freedom of the CV (error variance)}
  \item{side}{Side(s) to calculate the confidence limits for}
  \item{alpha}{Type I error probability, aka significance level}
}
\value{
Numeric vector of the confidence limits named as 'lower CL' and 'upper CL'.\cr
In case of the one-sided upper confidence limit the 'lower CL' is = 0.\cr 
In case of the one-sided lower confidence limit the 'upper CL' is = Inf. 
}
\author{D. Labes}
\examples{
# upper one-sided 95\% CL of a CV=0.3 
# from a study with df=22 (f.i. a 2x2 crossover with n=24)
# side="upper" is standard if not explicitly given
CVCL(0.3, df=22)
# should give:
# lower CL  upper CL 
#0.0000000 0.4075525 
}
