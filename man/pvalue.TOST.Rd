\name{pvalue.TOST}
\alias{pvalue.TOST}
\alias{pvalues.TOST}
\title{
p-value(s) of the TOST procedure
}
\description{
Calculates the p-value(s) of the TOST procedure via students t-distribution
given pe, CV and n.
}
\usage{
pvalue.TOST(pe, CV, n, logscale = TRUE, theta1, theta2, design = "2x2", 
            robust = FALSE, both = FALSE)
pvalues.TOST(pe, CV, n, logscale = TRUE, theta1, theta2, design = "2x2", 
             robust = FALSE, both = TRUE)
}
\arguments{
  \item{pe}{
Observed point estimate of the ratio Test vs. Reference (if \code{logscale=TRUE}) 
or of the difference (if \code{logscale=FALSE}).
}
  \item{CV}{
Obseved coefficient of variation as ratio or error standard deviation.
}
  \item{n}{
Total number of subjects if given as scalar.\cr
Number of subjects in (sequence) groups if given as vector.
}
  \item{logscale}{
Should the data used after log-transformation or on original scale? \cr
\code{TRUE} or \code{FALSE}. Defaults to \code{TRUE}.\cr
\code{logscale=FALSE} is useful if you have the data of analysis after 
log-transformation only.
}
  \item{theta1}{
Lower bioequivalence limit.\cr
In case of \code{logscale=TRUE} it has to be given as ratio, otherwise as value < 0.\cr
Defaults to 0.8 if \code{logscale=TRUE} or to log(0.8)= -0.2231 if \code{logscale=FALSE}.
}
  \item{theta2}{
Upper bioequivalence limit.\cr
If not given theta2 will be calculated as \code{1/theta1} if \code{logscale=TRUE}
or as \code{-theta1} if \code{logscale=FALSE}.
}
  \item{design}{
Character string describing the study design.\cr
See \code{known.designs()} for designs covered in this package.
}
  \item{robust}{
If set to \code{TRUE} triggers the use of degrees of freedom according to the 'robust' 
evaluation (aka Senn's basic estimator). These df are calculated as \code{n-seq}.\cr
See \code{known.designs()$df2}.Has only effect for higher-order crossover designs.\cr
Defaults to \code{FALSE}. With that value the usual degrees of freedom will be used.
}
  \item{both}{
Indicates if both p-values (t-tests of pe>=theta1 and pe<=theta2) shall be given 
back or only the maximum.\cr
Defaults to \code{FALSE} for the function \code{pvalue.TOST()} and to 
\code{TRUE} for the function \code{pvalue}\bold{s}\code{.TOST()}.
}
}
\value{
Returns the p-value(s).\cr
Returns a vector with named elements "p.left", "p.right" if arguments pe and CV 
are scalars, else a matrix with columns "p.left", "p.right". \cr
p.left names the p-value of testing HA1: theta>=theta1, p.right the p-value of 
testing HA2: theta<=theta2 against their respective Nulls.
}
\references{
Schuirmann, DJ\cr
"A comparison of the two one-sided tests procedure and the power approach for
assessing the equivalence of average bioavailability"\cr
J. Pharmacokin. Biopharm. 1987; 15: 657-680.

Hauschke D, Steinijans V, Pigeot I\cr
"Bioequivalence Studies in Drug Development"\cr
John Wiley & Sons Ltd, Chichester 2007
}
\author{
Benjamin Lang\cr
Man page by D. Labes
}
\note{
The formulas implemented cover balanced and unbalanced designs.\cr\cr
In case of argument \code{n} given as n(total) and is not divisible by the number
of (sequence) groups the total sample size is partitioned to the (sequence) 
groups to have small imbalance only. A message is given in such cases.\cr\cr
SAS procedure TTEST with the TOST option names p.left = Upper, p.right= Lower
according to the tail of the t-distribution to be used for obtaining the 
p-values.
}
\seealso{
\code{\link{CI.BE}}
}
\examples{
# Defaults: 2x2 crossover, log-transformation
# BE acceptance limits 0.8 ... 1.25, usual df's
# interested in both p-values
pvalues.TOST(pe=0.95, CV=0.3, n=12)
# gives the vector (named elements)
#     p.left    p.right
# 0.09105601 0.02250985
# i.e. 'left' hypothesis H01: theta<=theta1 can't be rejected
# 'right' hypothesis H02: theta>=theta2 can be rejected

# max. p-value only as 'overall' pvalue, preferred by Benjamin
pvalue.TOST(pe=0.912, CV=0.333, n=24)
# should give 0.08777621, i.e inequivalence can't be rejected
# this is operationally identical to 
CI.BE(pe=0.912, CV=0.333, n=24)
# lower limit = 0.7766 outside 0.8 ... 1.25, i.e inequivalence can't be rejected
}

