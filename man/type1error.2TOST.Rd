\name{type1error.2TOST}
\alias{type1error.2TOST}
\title{
Type I error rate for two simultaneous TOST procedures
}
\description{
Calculates the exact type I error rate of two simultaneous TOST procedures 
(where the two parameters of the two TOSTs are correlated with some correlation)
for various study designs used in BE studies 
}
\usage{
type1error.2TOST(alpha = c(0.05, 0.05), logscale = TRUE, theta1, theta2, 
                 CV, n, rho, design = "2x2", robust = FALSE, nsims, 
                 setseed = TRUE, details = FALSE)
}
\arguments{
  \item{alpha}{
Vector; contains one-sided significance level for each of the two TOSTs.\cr
For one TOST, by convention mostly set to 0.05.
}
  \item{logscale}{
Should the data used on log-transformed or on original scale? TRUE or FALSE.\cr
Defaults to TRUE.
}
  \item{theta1}{
Vector; contains lower bioequivalence limit for each of the two TOSTs.\cr
In case of \code{logscale=TRUE} it is given as ratio, otherwise as diff. to 1.\cr
Defaults to \code{c(0.8, 0.8)} if \code{logscale=TRUE} or to \code{c(-0.2, -0.2)}
if \code{logscale=FALSE}.
}
  \item{theta2}{
Vector; contains upper bioequivalence limit for each of the two TOSTS.\cr
If not given theta2 will be calculated as \code{1/theta1} if \code{logscale=TRUE} \cr
or as \code{-theta1} if \code{logscale=FALSE}.
}
  \item{CV}{
Vector of coefficient of variations (given as as ratio, e.g. 0.2 for 20\%).\cr
In case of cross-over studies this is the within-subject CV, \cr
in case of a parallel-group design the CV of the total variability.\cr
In case of \code{logscale=FALSE} CV is assumed to be the respective standard 
deviation.
}
  \item{n}{
Number of subjects under study.\cr
Is total number if given as scalar, else number of subjects in the (sequence) 
groups. In the latter case the length of n vector has to be equal to the 
number of (sequence) groups.
}
  \item{rho}{
Correlation between the two parameters under consideration. This is defined 
as correlation between the estimator of the treatment difference of
parameter one and the estimator of the treatment difference of parameter two.
}
  \item{design}{
Character string describing the study design.\cr
See \code{known.designs()} for designs covered in this package.
}
  \item{robust}{
Defaults to FALSE. With that value the usual degrees of freedom will be used.\cr
Set to \code{TRUE} will use the degrees of freedom according to the 'robust' 
evaluation (aka Senn's basic estimator). These df are calculated as \code{n-seq}.\cr
See \code{known.designs()$df2} for designs covered in this package.\cr
Has only effect for higher-order crossover designs.
  }
  \item{setseed}{
Calculation depends on \code{pmvt()} which is based on randomized quasi Monte Carlo
methods. If \code{setseed=TRUE} a seed value is set, the default.
}
  \item{details}{
logical; if TRUE, P(Type I error) will be returned for all eight intersection 
null sets, otherwise the maximum (default).
}
}
\details{
The exact calculations of the type 1 error rate are performed via integration of
the 4-dimensional non-central t-distribution via function \code{pmvt()} of package 
\code{mvtnorm}. An absolute error tolerance of 1e-05 is set within \code{pmvt()}.\cr\cr
The formulas cover balanced and unbalanced studies w.r.t (sequence) groups.\cr\cr
In case of parallel group design and higher order crossover designs 
(replicate crossover or crossover with more than two treatments) the calculations 
are based on the assumption of equal variances for Test and Reference products
under consideration.\cr\cr
The formulas for the paired means 'design' do not take an additional correlation
parameter into account. They are solely based on the paired t-test 
(TOST of differences = zero).
}
\value{
Value of Type I Error rate if argument \code{details = FALSE}.\cr
A data.frame with the TIE for each of the eight nullhypothesis regions if
\code{details = TRUE}.
}
\note{
If n is given as scalar (total sample size) and this number is not 
divisible by the number of (sequence) groups of the design an unbalanced design 
with small imbalance is assumed. A corresponding message is thrown showing the 
assumed numbers of subjects in (sequence) groups.\cr
The function does not vectorize properly if design is a vector. Moreover,
CV must be of length two, thus further vectorizing is not possible.\cr
Other vector input is not tested yet.\cr\cr
The calculation of the TIE may last some seconds. Be patient.
}
\references{
Hua SY, Xu S, D'Agostino RB Sr.\cr
\emph{Multiplicity adjustments in testing for bioequivalence}\cr
Stat Med. 2015;34(2):215--31. \doi{10.1002/sim.6247}

Lang B, Fleischer F.\cr
\emph{Letter to the Editor: Comments on 'Multiplicity adjustments in testing for bioequivalence'}\cr
Stat Med. 2016;35(14):2479--80. \doi{10.1002/sim.6488}
}
\author{
Benjamin Lang\cr
}
\seealso{
 \code{\link{sampleN.2TOST}, \link{known.designs}}
}
\examples{
\dontrun{
# Replicate type 1 error rate for scenario S2 from Hua et al.
# runs 6-7 seconds, more than allowed for examples on CRAN
n <- 24
cv <- se2CV(c(0.0490 / sqrt(2/n), 0.0657 / sqrt(2/n)))
type1error.2TOST(CV = cv, n = n, rho = 0.6794, details = FALSE)}
}
