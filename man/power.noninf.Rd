\name{power.noninf}
\alias{power.noninf}
\title{
Power of the one-sided non-inferiority t-test
}
\description{
Function calculates of the power of the one-sided non-inferiority t-test
for normal or log-normal distributed data.
}
\usage{
power.noninf(alpha = 0.025, logscale = TRUE, margin, theta0, CV, n, 
             design = "2x2", robust = FALSE)
}
\arguments{
  \item{alpha}{
Type I error probability, significance level. Defaults here to 0.025.
}
  \item{logscale}{
Should the data used on log-transformed or on original scale? TRUE or FALSE.\cr
Defaults to TRUE.
}
  \item{margin}{
Non-inferiority margin.\cr
In case of \code{logscale=TRUE} it must be given as ratio, otherwise as diff. to 1.\cr
Defaults to 0.8 if \code{logscale=TRUE} or to -0.2 if \code{logscale=FALSE}.
}
  \item{theta0}{
'True' or assumed bioequivalence ratio or difference. \cr
In case of \code{logscale=TRUE} it must be given as ratio, \cr
otherwise as difference to 1. See examples.\cr
Defaults to 0.95 if \code{logscale=TRUE} or to -0.05 if \code{logscale=FALSE}.
}
  \item{CV}{
Coefficient of variation as ratio.\cr
In case of cross-over studies this is the within-subject CV, \cr
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
See \code{\link{known.designs}} for designs covered in this package.
}
  \item{robust}{
Defaults to FALSE. With that value the usual degrees of freedom will be used.\cr
Set to \code{TRUE} will use the degrees of freedom according to the 'robust' evaluation
(aka Senn's basic estimator). These df are calculated as \code{n-seq}.\cr
See \code{known.designs()$df2} for designs covered in this package.\cr
Has only effect for higher-order crossover designs.
}
}
\details{
The power is calculated via non-central t-distribution.
}
\value{
Value of power according to the input arguments.
}
\section{Warning }{
The function does not vectorize if design is a vector.\cr
The function vectorize properly if CV or theta0 are vectors.\cr
Other vector input is not tested yet.
}
\references{
S.A. Julious\cr
"TUTORIAL IN BIOSTATISTICS\cr
Sample sizes for clinical trials with Normal data"\cr
Statist. Med. 2004; 23: 1921-1986}
\author{
D. Labes
}
\note{
This function does not rely on TOST but may be useful in planning BE studies 
if the question is not equivalence but 'non-superiority'. \cr
Hint: Evaluation of Fluctuation in the EMA MR NfG (1999) between modified release 
formulation and immediate release product.
}

\seealso{
\code{\link{known.designs}, \link{sampleN.noninf}}
}
\examples{
# using all the defaults: margin=0.8, theta0=0.95, alpha=0.025
# log-transformed, design="2x2"
# should give: 0.4916748
power.noninf(CV=0.3, n=24)
}
