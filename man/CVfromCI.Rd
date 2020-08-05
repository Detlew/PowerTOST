\encoding{utf-8}
\name{CVfromCI}
\alias{CVfromCI}
\alias{CI2CV}
\title{CV from a given Confidence interval}
\description{
Calculates the CV (coefficient of variation) from a known confidence interval
of a BE study.\cr 
Useful if no CV but the 90\% CI was given in literature.
}
\usage{
CVfromCI(pe, lower, upper, n, design = "2x2", alpha = 0.05, robust = FALSE)
CI2CV(pe, lower, upper, n, design = "2x2", alpha = 0.05, robust = FALSE)
}
\arguments{
  \item{pe}{
Point estimate of the T/R ratio.\cr
The \code{pe} may be missing. In that case it will be calculated as geometric 
mean\cr
of \code{lower} and \code{upper}.
  }
  \item{lower}{
Lower confidence limit of the BE ratio.
  }
  \item{upper}{
Upper confidence limit of the BE ratio.
  }
  \item{n}{
Total number of subjects under study if given as scalar.\cr
Number of subjects in (sequence) groups if given as vector.
  }
  \item{design}{
Character string describing the study design.\cr
See \code{known.designs()} for designs covered in this package.
  }
  \item{alpha}{
Error probability. Set it to \code{(1-confidence)/2} (\emph{i.e.} to 0.05 for the usual 90\% confidence intervals).
  }
  \item{robust}{
With \code{robust=FALSE} (the default) usual degrees of freedom of the designs are used.\cr
With \code{robust=TRUE} the degrees of freedom for the so-called robust
evaluation (df2 in known.designs()) will be used. This may be helpful
if the CI was evaluated via a mixed model or via intra-subject contrasts 
(aka \enc{Senn’s}{Senn's} basic estimator).
  }
}
\details{
See Helmut \enc{Schütz’}{Schuetz'} \href{https://bebac.at/lectures/Bucharest2013WS1.pdf#page=26}{presentation} for the algebra underlying this function.
}
\value{
Numeric value of the CV as ratio.
}

\references{
Yuan J, Tong T, Tang M-L. \emph{Sample Size Calculation for Bioequivalence Studies Assessing Drug Effect and Food Effect at the Same Time With a 3-Treatment Williams Design.} Regul Sci. 2013;47(2):242--7. \doi{10.1177/2168479012474273}
}

\author{Original by D. Labes with suggestions by H. \enc{Schütz}{Schuetz}.\cr
Reworked and adapted to unbalanced studies by B. Lang.}
\note{
The calculations are based on the assumption of evaluation via log-transformed values.\cr
The calculations are further based on a common variance of Test and Reference 
treatments in replicate crossover studies or parallel group study, respectively.\cr\cr
In case of argument \code{n} given as n(total) and is not divisible by the number
of (sequence) groups the total sample size is partitioned to the (sequence) groups
to have small imbalance only. A message is given in such cases.\cr
The estimated CV is conservative (\emph{i.e.}, higher than actually observed) in case of 
unbalancedness.\cr\cr
\code{CI2CV()} is simply an alias to \code{CVfromCI()}.
}

\examples{
# Given a 90\% confidence interval (without point estimate) 
# from a classical 2x2 crossover with 22 subjects
CVfromCI(lower=0.91, upper=1.15, n=22, design="2x2")
# will give [1] 0.2279405, i.e a CV ~ 23\%
#
# unbalanced 2x2 crossover study, but not reported as  such
CI2CV(lower=0.89, upper=1.15, n=24)
# will give a CV ~ 26.3\%
# unbalancedness accounted for
CI2CV(lower=0.89, upper=1.15, n=c(16,8))
# should give CV ~ 24.7\%
}