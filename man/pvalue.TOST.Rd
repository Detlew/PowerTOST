\encoding{utf-8}
\name{pvalue.TOST}
\alias{pvalue.TOST}
\alias{pvalues.TOST}
\title{
p-value(s) of the TOST procedure
}
\description{
Calculates the \emph{p}-value(s) of the TOST procedure via students \emph{t}-distribution
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
Observed point estimate of the T/R ratio (if \code{logscale=TRUE}) 
or of the difference T--R (if \code{logscale=FALSE}).
}
  \item{CV}{
Obseved coefficient of variation as ratio (if \code{logscale=TRUE}) 
or residual error standard deviation (if \code{logscale=FALSE}).
}
  \item{n}{
Total number of subjects if given as scalar.\cr
Number of subjects in (sequence) groups if given as vector.
}
  \item{logscale}{
Should the data be used after log-transformation or on original scale? \cr
\code{TRUE} or \code{FALSE}. Defaults to \code{TRUE}.
}
  \item{theta1}{
Lower bioequivalence limit.\cr
In case of \code{logscale=TRUE} it has to be given as ratio, otherwise as value < 0.\cr
Defaults to 0.8 if \code{logscale=TRUE} or to \verb{log(0.8) ~ -0.2231} if \code{logscale=FALSE}.
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
If set to \code{TRUE} triggers the use of degrees of freedom according to the \sQuote{robust} 
evaluation (aka Senn’s basic estimator). These degrees of freedom are calculated as \code{n-seq}.\cr
See \code{known.designs()$df2}. Has only effect for higher-order crossover designs.\cr
Defaults to \code{FALSE}. With that value the usual degrees of freedom will be used.
}
  \item{both}{
Indicates if both \emph{p}-values (\emph{t}-tests of pe >= theta1 and pe <= theta2) shall be given 
back or only the maximum.\cr
Defaults to \code{FALSE} for the function \code{pvalue.TOST()} and to 
\code{TRUE} for the function \code{pvalue}\bold{s}\code{.TOST()}.
}
}
\value{
Returns the \emph{p}-value(s).\cr
Returns a vector with named elements \code{p.left}, \code{p.right} if arguments \code{pe} and \code{CV} 
are scalars, else a matrix with columns \code{p.left}, \code{p.right}. \cr
p.left gives the \emph{p}-value of testing\cr
\verb{  HA1: theta >= theta1}\cr
p.right the \emph{p}-value of testing\cr
\verb{  HA2: theta <= theta2}\cr
against their respective Nulls.
}
\references{
Schuirmann DJ. \emph{A comparison of the two one-sided tests procedure and the power approach for
assessing the equivalence of average bioavailability.} J Pharmacokin Biopharm. 1987;15:657--80. \doi{10.1007/BF01068419}

Hauschke D, Steinijans V, Pigeot I. \emph{Bioequivalence Studies in Drug Development.} Chichester: Wiley; 2007.
}
\author{
B. Lang, man page by D. Labes
}
\note{
The formulas implemented cover balanced and unbalanced designs.\cr\cr
In case of argument \code{n} given as n(total) and is not divisible by the number
of (sequence) groups the total sample size is partitioned to the (sequence) 
groups to have small imbalance only. A message is given in such cases.\cr\cr
SAS procedure TTEST with the TOST option names p.left = Upper, p.right= Lower
according to the tail of the \emph{t}-distribution to be used for obtaining the 
\emph{p}-values.
}
\seealso{
\code{\link{CI.BE}}
}
\examples{
# Defaults: 2x2 crossover, log-transformation
# BE acceptance limits 0.8 ... 1.25, usual dfs
# interested in both p-values
pvalues.TOST(pe = 0.95, CV = 0.3, n = 12)
# gives the vector (named elements)
#     p.left    p.right
# 0.09105601 0.02250985
# i.e. 'left' hypothesis H01: theta<=theta1 can't be rejected
# 'right' hypothesis H02: theta>=theta2 can be rejected

# max. p-value only as 'overall' pvalue, preferred by Benjamin
pvalue.TOST(pe = 0.912, CV = 0.333, n = 24)
# should give 0.08777621, i.e., inequivalence can't be rejected
# this is operationally identical to 
CI.BE(pe = 0.912, CV = .333, n = 24)
# lower limit = 0.7766 outside 0.8 ... 1.25, i.e., inequivalence can't be rejected
}

