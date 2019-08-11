\encoding{utf-8}
\name{CVCL}
\alias{CVCL}
\title{
Confidence limits of a CV for log-normal data
}
\description{
The function calculates the \ifelse{html}{\out{1&ndash;&alpha;}}{\eqn{1-\alpha}{1-alpha}} confidence limits
(either 1-sided or 2-sided) via the \ifelse{html}{\out{&chi;<sup>2</sup>}}{\eqn{\chi^2}{chi-squared}}
distribution of the error variance the CV is based on.
}
\usage{
CVCL(CV, df, side = c("upper", "lower", "2-sided"), alpha = 0.05)
}
\arguments{
  \item{CV}{Coefficient of variation}
  \item{df}{degrees of freedom of the CV (error variance)}
  \item{side}{Side(s) to calculate the confidence limits for, defaults to \code{upper}}
  \item{alpha}{Type I error probability, aka significance level}
}
\value{
Numeric vector of the confidence limits named as \code{lower CL} and \code{upper CL}.\cr
In case of the one-sided upper confidence limit the \code{lower CL} is = 0.\cr 
In case of the one-sided lower confidence limit the \code{upper CL} is = Inf. 
}
\author{D. Labes}
\examples{
# upper one-sided 95\% CL of a CV=0.3 
# from a study with df=22 (f.i. a 2x2 crossover with n=24)
# default side="upper" since not explicitly given
CVCL(0.3, df = 22)
# should give:
#  lower CL  upper CL 
# 0.0000000 0.4075525 
}
