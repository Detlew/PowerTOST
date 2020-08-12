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
Observed point estimate of the T/R ratio or difference.\cr
In case of \code{logscale=TRUE} it must be given as ratio T/R.\cr
If \code{logscale=FALSE}, the observed difference in means. In this case, the difference may be expressed in two ways: relative to the same (underlying) reference mean, i.e. as (T-R)/R = T/R - 1; or as difference in means T-R. Note that in the former case the units of \code{CV}, \code{theta1} and \code{theta2} need also be given relative to the reference mean (specified as ratio).\cr
}
  \item{CV}{
In case of \code{logscale=TRUE} the observed (geometric) coefficient of variation given as ratio.\cr
If \code{logscale=FALSE} the argument refers to the observed (residual) standard deviation of the response. In this case, standard deviation may be expressed two ways: relative to a reference mean (specified as ratio sigma/muR), i.e. again as a coefficient of variation; or untransformed, i.e. as standard deviation of the response. Note that in the former case the units of \code{pe}, \code{theta1} and \code{theta2} need also be given relative to the reference mean (specified as ratio).\cr

In case of cross-over studies this is the within-subject CV, in case of a parallel-group design the CV of the total variability.
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
Lower (bio-)equivalence limit.\cr
In case of \code{logscale=TRUE} it is given as ratio.\cr
If \code{logscale=FALSE}, the limit may be expressed in two ways:
difference of means relative to the same (underlying) reference mean or in units of the difference of means.
Note that in the former case the units of \code{CV}, \code{pe} and \code{theta2} need also be given relative to the reference mean (specified as ratio).\cr
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
  \item{design}{
Character string describing the study design.\cr
See \code{known.designs()} for designs covered in this package.
}
  \item{robust}{
If set to \code{TRUE} triggers the use of degrees of freedom according to the \sQuote{robust} 
evaluation (aka \enc{Senn’s}{Senn's} basic estimator). These degrees of freedom are calculated as \code{n-seq}.\cr
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
\code{p.left} gives the \emph{p}-value of testing\cr
\verb{  HA1: theta >= theta1}\cr
and \code{p.right} the \emph{p}-value of testing\cr
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

