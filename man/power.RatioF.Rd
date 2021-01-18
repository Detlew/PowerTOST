\encoding{utf-8}
\name{power.RatioF}
\alias{power.RatioF}
\title{
Power for equivalence of the ratio of two means with normality on original scale
}
\description{
Calculates the power of the test of equivalence of the ratio of two means 
with normality on original scale.\cr 
This test is based on Fieller’s confidence (\sQuote{fiducial}) interval and Sasabuchi’s 
test (a TOST procedure as well).
}
\usage{
power.RatioF(alpha = 0.025, theta1 = 0.8, theta2, theta0 = 0.95,
             CV, CVb, n, design = "2x2", setseed=TRUE)
}
\arguments{
  \item{alpha}{
Type I error probability, aka significance level.\cr
Defaults here to 0.025 because this function is intended for studies
with clinical endpoints. 
}
  \item{theta1}{
Lower bioequivalence limit. Typically 0.8 (default).
}
  \item{theta2}{
Upper bioequivalence limit. Typically 1.25.\cr
Is set to \code{1/theta1} if missing.
}
  \item{theta0}{
\sQuote{True} or assumed T/R ratio. Typically set to 0.95 for planning.
}
  \item{CV}{
Coefficient of variation as ratio. In case of \code{design="parallel"} this is
the CV of the total variability, in case of \code{design="2x2"} the 
intra-subject CV (CVw in the reference).
}
  \item{CVb}{
CV of the between-subject variability. Only necessary for \code{design="2x2"}.
}
  \item{n}{
Number of subjects to be planned.\cr
\code{n} is for both designs implemented the \bold{total} number of subjects.\cr
}
  \item{design}{
A character string describing the study design.\cr 
\code{design="parallel"} or \code{design="2x2"} allowed for a two-parallel 
group design or a classical TR|RT crossover design.
}
  \item{setseed}{
If set to \code{TRUE} the dependence of the power from the state of the random number
generator is avoided. With \code{setseed = FALSE} you may see the dependence
from the state of the random number generator.
}
}
\details{
The power is calculated exact using the bivariate non-central \emph{t}-distribution
via function \code{\link[mvtnorm]{pmvt}} of the package \code{mvtnorm}.\cr
Due to the calculation method of the used package mvtnorm -- randomized 
Quasi-Monte-Carlo -- these probabilities are dependent from the state of the 
random number generator within the precision of the power. 
See argument \code{setseed}.
}
\value{
Value of power according to the input.
}
\references{
Fieller EC. \emph{Some Problems in Interval Estimation.} J Royal Stat Soc B. 1954;16(2):175--85. \doi{10.1111/j.2517-6161.1954.tb00159.x}

Sasabuchi S. \emph{A test of a multivariate normal mean with composite hypotheses determined by linear inequalities.} Biometrika. 1980;67(2):429--39. \doi{10.1093/biomet/67.2.429}

Hauschke D, Kieser M, Diletti E, Burke M. \emph{Sample size determination for proving equivalence based on the ratio of two means for normally distributed data.} Stat Med. 1999;18(1):93--105.

Hauschke D, Steinijans V, Pigeot I. \emph{Bioequivalence Studies in Drug Development.} Chichester: Wiley; 2007. Chapter 10.

European Agency for the Evaluation of Medicinal Products, CPMP. \emph{Points to Consider on Switching between Superiority and Non-Inferiority.} London, 27 July 2000. \href{https://www.ema.europa.eu/en/documents/scientific-guideline/points-consider-switching-between-superiority-non-inferiority_en.pdf}{CPMP/EWP/482/99}
}
\author{
D. Labes
}
\note{
This function is intended for studies with clinical endpoints where the 95\% confidence intervals are usually used for equivalence testing.\cr
Therefore, alpha defaults here to 0.025 (see EMEA 2000).\cr\cr
The formulas given in the references rely on the assumption of equal variances
in the two treatment groups for the parallel group design or on assuming equal 
within-subject and between-subject variabilities for the 2×2 crossover design.
}

\seealso{
\code{\link{sampleN.RatioF}}
}
\examples{
# power for alpha=0.025, ratio0=0.95, theta1=0.8, theta2=1/theta1=1.25
# within-subject CV=0.2, between-subject CV=0.4 
# 2x2 crossover study, n=24
# using all the defaults:
power.RatioF(CV = 0.2, CVb = 0.4, n = 24)
# gives [1] 0.7315357
}
