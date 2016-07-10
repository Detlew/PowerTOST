\encoding{utf-8}
\name{scABEL.ad}
\alias{scABEL.ad}
\title{
Iteratively adjusted alpha for ABEL
}
\description{
This function iteratively adjusts alpha for the \acronym{BE} decision
via Average Bioequivalence with Expanding Limits (\acronym{ABEL}) based
on simulations in order to maintain the consumer risk at the nominal level.
}
\usage{
scABEL.ad(alpha = 0.05, theta0, theta1, theta2, CV = 0.3,
          design = c("2x3x3", "2x2x4", "2x2x3"), regulator, n, alpha.pre = 0.05,
          imax = 100, tol, print = TRUE, details = FALSE, 
          setseed = TRUE, nsims = 1e6)
}
\arguments{
  \item{alpha}{
Type I error (\acronym{TIE}) probability (nominal level of the test). Per
convention commonly set to 0.05.\cr
A lower value might be specified (e.g., in order to adjust for multiplicity
in dose proportionality studies).
}
  \item{theta0}{
\sQuote{True} or assumed bioavailability ratio.\cr
Defaults to 0.90 if not given explicitly.
}
  \item{theta1}{
Conventional lower \acronym{ABE} limit to be applied in the mixed procedure
if \code{CVwR} = \code{CVswitch}. Also lower limit for the point estimate
constraint.\cr
Defaults to 0.80 if not given explicitly.
}
  \item{theta2}{
Conventional upper \acronym{ABE} limit to be applied in the mixed procedure
if \code{CVwR} = \code{CVswitch}. Also upper limit for the point estimate
constraint.\cr
Defaults to 1.25 if not given explicitly.
}
  \item{CV}{
Coefficient(s) of variation as ratio (not percent). Defaults to
0.30 (i.e., the location of maximum inflation of the \acronym{TIE} for the
\acronym{EMA}'s \acronym{ABEL}.\cr
If \code{length(CV) = 1} the same CV of Test and Reference is assumed
(\code{CVwT} = \code{CVwR}).\cr
If \code{length(CV) = 2} the CV of Test \emph{must} be given in \code{CV[1]}
and the one of Reference in \code{CV[2]}.
}
  \item{design}{
Design of the study to be planned.\cr
\code{"2x3x3"} is the partial replicate design (RRT|RTR|TRR).\cr
\code{"2x2x3"} is the 2-sequence 3-period full replicate design (RTR|TRT).\cr
\code{"2x2x4"} is the 2-sequence 4-period full replicate design (RTRT|TRTR).\cr
Defaults to \code{design = "2x3x3"}.
}
  \item{regulator}{
Regulatory settings for the widening of the \acronym{BE} acceptance limits.\cr
Choose from \code{c("EMA", "HC")}. Defaults to \code{regulator = "EMA"}.\cr
This argument may be given also in lower case.\cr
The former \code{regulator="ANVISA"} is also allowed but is internally changed 
to "EMA" since ANVISA now recommends the use of EMA regulatory settings. 
Don't use "ANVISA" any longer since it will be removed in next versions.
A corresponding message informs about this.
}
  \item{n}{
Total sample size of the study or a vector of sample size / sequences.
If \code{n} leads to an unbalanced design (i.e., is not a multiple of two
in the full replicate designs or not a multiple of three in the
partial replicate), the code tries to keep subjects / sequence as
balanced as possible.\cr
In evaluating a particular \emph{unbalanced} study \strong{always} give \code{n}
as a vector.\cr
If n is missing a sample size is estimated with targetpower = 0.8 and pre-specified
alpha if defined. Else alpha is used.
}
  \item{nsims}{
Number of simulations to be performed to estimate the (empirical)
\acronym{TIE} error and in each iteration of adjusting alpha.\cr
The default value 1,000,000 = 1e+6 should not be lowered.
}
  \item{imax}{
Maximum number of steps in sample size search. Defaults to 100.
}
  \item{tol}{
Desired accuracy (convergence tolerance).\cr
Defaults to 1e-6.
}
  \item{print}{
If \code{TRUE} (default), the function prints its results.
}
  \item{details}{
If \code{TRUE}, the relative change of the consumer risk is shown.
Additionally information about the impact on power (for specified \code{theta0}
and target power 0.80), runtime, and number of simulations (iterations)
are given.\cr
Defaults to \code{FALSE}.
}
  \item{alpha.pre}{
Pre-specified alpha (optional). Must be <=\code{alpha}. \acronym{ABEL} will be
performed at level \code{alpha.pre} and the \acronym{TIE} assessed at level \code{alpha}.\cr
Less powerful than adjusting alpha but an alternative in the critical
region of maximum inflation of the \acronym{TIE}. Not recommended for \code{CVwR}
>=0.45 due to poor power characteristics.
}
\item{setseed}{
Simulations are dependent on the starting point of the (pseudo)
random number generator. To avoid differences in power for different
runs a \code{set.seed(123456)} is issued if \code{setseed=TRUE} (default).
}
}
\details{
The simulations are done via the distributional properties of the
statistical quantities necessary for assessing \acronym{BE} based on \acronym{ABEL}.\cr
Simulations for the \acronym{TIE} are performed at the upper (expanded) limit \emph{U}
of the acceptance range.  Examples of \emph{U} at the \acronym{EMA}'s 
\code{CVswitch} and \code{CVcap}:\cr
\code{scABEL(CV=0.3, regulator="EMA")[["upper"]]}\cr
\code{[1] 1.25}\cr
\code{scABEL(CV=0.5, regulator="EMA")[["upper"]]}\cr
\code{[1] 1.43191}\cr
Due to the symmetry around 1 results are valid for the lower (expanded)
limit \emph{L} as well.\cr
If an inflation of the \acronym{TIE} is expected (i.e., > \code{alpha}), alpha 
is iteratively adjusted until the consumer risk is maintained (<= \code{alpha}).\cr
}
\value{
Returns a list with the input, adjusted alpha, and type I error
(for nominal and adjusted alpha).\cr
If no adjustment is necessary, \code{NA} will be returned for
\code{adj. alpha} and the \acronym{TIE} for \code{alpha0} (or \code{alpha1}, 
if applicable) in \code{TIE.unadj}.
}
\references{
Wonnemann, M., \enc{Frömke}{Froemke}, C., Koch, A.\cr
Inflation of the Type I Error: Investigations on Regulatory Recommendations for 
Bioequivalence of Highly Variable Drugs\cr
Pharm Res. 2015;32(1):135-43. \href{http://dx.doi.org/10.1007/s11095-014-1450-z}{doi 10.1007/s11095-014-1450-z}

\enc{Muñoz}{Munoz} J., Alcaide D., \enc{Ocaña}{Ocana} J.\cr 
Consumer’s risk in the EMA and FDA regulatory approaches for bioequivalence
in highly variable drugs\cr
Stat Med. 2015;35(12):1933-43. \href{http://dx.doi.org/10.1002/sim.6834}{doi 10.1002/sim.6834}

Labes, D., \enc{Schütz}{Schuetz}, H.\cr
Inflation of Type I Error in the Evaluation of Scaled Average Bioequivalence, 
and a Method for its Control\cr
submitted to Pharm. Res. 2016
}
\author{
H. \enc{Schütz}{Schuetz}
}
\note{
The type I error in \acronym{ABEL} depends only on \code{CVwR} and - to
a minor degree - on the sample size. The \acronym{TIE} is assessed twice:\cr
In the first step based on \code{alpha} (or \code{alpha.pre}) and compared to \code{alpha}.
If no inflation of the TIE is found, the algo stops.\cr
Otherwise, alpha is iteratively adjusted (alway towards lower values)
until no more inflation is detected. In the final assessment \code{alpha.adj} is
compared to the nominal level of the test \code{alpha}.\cr

Specifying \code{theta0} is not necessary.\cr
If \code{theta0} \emph{is not} given, achievable power for the
common target of 0.80 (both for \code{alpha} and adjusted alpha) will be
estimated. If \code{theta0} \emph{is} specified, its value will be used; again for
target power 0.80.\cr
If you are interested in other levels of power, use \code{\link{sampleN.scABEL.ad}}.
}

\section{Warning }{
See the Warning section of the function \code{\link{power.scABEL}} concerning
the power value agreement to those obtained from simulations via subject data.
}

\seealso{
\code{\link{sampleN.scABEL.ad}}, \code{\link{power.scABEL}}, \code{\link{scABEL}}
}
\examples{
# partial replicate design, target power=80\%,
# true assumed ratio = 0.90, 1E+6 simulated studies
# EMA regulatory settings (ABE limits, PE constraint 0.8 - 1.25)
# Not run: due to timing policy of CRAN for examples
\dontrun{
scABEL.ad(CV = 0.3)}
# should result in adjusted alpha 0.03389 (TIE 0.5000, TIE for nominal alpha 0.07189).
#
# EMA, full replicate design, CV 0.35, sample size 33 (unbalanced)
# Not run: due to timing policy of CRAN for examples
\dontrun{
scABEL.ad(CV = 0.35, design = "2x2x4", n = c(16, 17))}
# should result in adjusted alpha 0.03632 (TIE 0.5000, TIE for nominal alpha 0.06544).
}
