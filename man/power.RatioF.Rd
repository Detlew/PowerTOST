\name{power.RatioF}
\alias{power.RatioF}
\title{
Power for equivalence of the ratio of two means with normality on original scale
}
\description{
Calculates the power of the test of equivalence of the ratio of two means 
with normality on original scale.\cr 
This test is based on Fieller's confidence ('fiducial') interval and Sasabuchi's 
test (again a TOST procedure).
}
\usage{
power.RatioF(alpha = 0.025, theta1 = 0.8, theta2, theta0 = 0.95, CV, CVb, n,   
             design = "2x2", setseed=TRUE)
}
\arguments{
  \item{alpha}{
Type I error probability, aka significance level. \cr
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
True ('null') assumed bioequivalence ratio. Typically set to 0.95 for planning.
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
n is for both designs implemented the \bold{total} number of subjects.\cr
}
  \item{design}{
A character string describing the study design.\cr 
\code{design="parallel"} or \code{design="2x2"} allowed for a two-parallel 
group design or a classical TR/RT crossover design.
}
  \item{setseed}{
If set to TRUE the dependence of the power from the state of the random number
generator is avoided. With \code{setseed = FALSE} you may see the dependence
from the state of the random number generator.
}
}
\details{
The power is calculated exact using the bivariate non-central t-distribution
via function \code{pmvt()} from the package mvtnorm.\cr
Due to the calculation method of the used package mvtnorm - randomized 
Quasi-Monte-Carlo - these probabilities are dependent from the state of the 
random number generator within the precision of the power. 
See argument \code{setseed}.
}
\value{
Value of power according to the input.
}
\references{
Hauschke D, Kieser M, Diletti E, Burke M.\cr
\emph{Sample size determination for proving equivalence based on the ratio of two means for normally distributed data}\cr
Stat Med. 1999;18(1):93--105. \doi{10.1002/(SICI)1097-0258(19990115)18:1<93::AID-SIM992>3.0.CO;2-8}

Hauschke D, Steinijans V, Pigeot I.\cr
\emph{Bioequivalence Studies in Drug Development}\cr
Chichester: Wiley; 2007. Chapter 10.
}
\author{
D. Labes
}
\note{
This function is intended for studies with clinical endpoints. \cr 
In such studies the 95\% confidence intervals are usually used for equivalence testing.\cr
Therefore alpha defaults here to 0.025.\cr
See CPMP/EWP/482/99 "Points to consider on switching between superiority and 
non-inferiority" EMEA, London (2000).\cr\cr
The formulas given in the references rely on the assumption of equal variances
in the two treatment groups for the parallel group design or on assuming equal 
within-subject and between-subject variabilities for the 2x2 crossover design.
}

\seealso{
\code{\link{sampleN.RatioF}}
}
\examples{
# power for alpha=0.025, ratio0=0.95, theta1=0.8, theta2=1/theta1=1.25
# within-subject CV=0.2, between-subject CV=0.4 
# 2x2 crossover study, n=24
# using all the defaults:
power.RatioF(CV=0.2, CVb=0.4, n=24)
# gives [1] 0.7315357
}
