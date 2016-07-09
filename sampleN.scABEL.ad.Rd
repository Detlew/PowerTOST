\encoding{utf-8}
\name{sampleN.scABEL.ad}
\alias{sampleN.scABEL.ad}
\title{
Sample size estimation for ABEL and iteratively adjusted alpha
}
\description{
This function performs a sample size estimation for the BE decision
via Average Bioequivalenc with Expanding Limits (\acronym{ABEL}) based
on simulations. Simultaneously alpha is iteratively adjusted in order
to maintain the consumer risk at the nominal level.
}
\usage{
sampleN.scABEL.ad(alpha = 0.05, targetpower = 0.8, theta0, theta1, theta2,
                  CV = 0.3, design = c("2x3x3", "2x2x4", "2x2x3"), regulator,
                  nstart = NA, nsims = 1e6, imax = 100, tol, print = TRUE,
                  details = FALSE, alpha.pre = 0.05, setseed = TRUE)
}
\arguments{
  \item{alpha}{
Type I error (\acronym{TIE}) probability (nominal level of the test). Per
convention commonly set to 0.05.\cr
A lower value might be specified (e.g., in order to adjust for multiplicity
in dose proportionality studies).
}
  \item{targetpower}{
Power to achieve at least. Must be \{=>0, <=1\}. Typical
values are 0.80 to 0.90 (i.e., 80\% to 90\%).\cr
Defaults to 0.80 if not given explicitly.
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
0.30 (i.e., the location of maximum inflation of the TIE for the
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
  \item{nstart}{
Best \dQuote{guess} sample size. If not given (default), simulations
start with the sample size estimated for \code{alpha} (or \code{alpha.pre},
if given), \code{theta0}, and \code{targetpower}.\cr
Can also be set to start the sample size search if a previous run
failed.\cr
According to regulatory requirements must be >=12 for the \acronym{EMA} and
>=24 for \acronym{ANVISA}.
}
  \item{nsims}{
Number of simulations to be performed to estimate the (empirical)
\acronym{TIE} and in each iteration of adjusting alpha.\cr
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
If \code{TRUE} (default), the steps during sample size search are
shown. Additionally information about the impact on power by
adjusting alpha and change of study costs due to the increased sample
size is given.
}
  \item{alpha.pre}{
Pre-specified alpha (optional). Must be <=\code{alpha}. \acronym{ABEL} will be
performed at level \code{alpha.pre} and the \acronym{TIE} assessed at level \code{alpha}.\cr
Less powerful than adjusting alpha but an alternative in the critical
region of maximum inflation of the \acronym{TIE}. Not recommended for \code{CVwR}
>= 0.45 due to poor power characteristics.
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
of the acceptance range. Examples of \emph{U} at \acronym{EMA}'s \code{CVswitch} and \code{CVcap}:\cr
\code{scABEL(CV=0.3, regulator="EMA")[["upper"]]}\cr
\code{[1] 1.25}\cr
\code{scABEL(CV=0.5, regulator="EMA")[["upper"]]}\cr
\code{[1] 1.43191}\cr
Due to the symmetry around 1 results are valid for the lower (expanded)
limit \emph{L} as well.\cr

If an inflation of the \acronym{TIE} is expected (i.e., > \code{alpha}), alpha is iteratively adjusted
until at least the target power is reached and the consumer risk is maintained (<= \code{alpha})

}
\value{
Returns a data.frame with the input and results for adjusted alpha,
type I error, sample size, and achieved power.\cr
The \dQuote{Sample size} column contains the total sample size.
If no adjustment is necessary, \code{NA} will be returned in the
\dQuote{adj. alpha} column and other results are identical to the ones
obtained by \code{\link{sampleN.scABEL}}.
}
\references{
\enc{Tóthfalusi}{Tothfalusi}, L., \enc{Endrényi}{Endrenyi}, L.\cr
Sample Sizes for Designing Bioequivalence Studies for Highly Variable Drugs\cr
J Pharm Pharmaceut Sci. 2011;15(1):73-84. \href{http://ejournals.library.ualberta.ca/index.php/JPPS/article/download/11612/9489}{free download}

Wonnemann, M., \enc{Frömke}{Froemke}, C., Koch, A.\cr
Inflation of the Type I Error: Investigations on Regulatory Recommendations for 
Bioequivalence of Highly Variable Drugs\cr
Pharm Res. 2015;32(1):135-43. \href{http://dx.doi.org/10.1007/s11095-014-1450-z}{doi 10.1007/s11095-014-1450-z}

\enc{Muñoz}{Munoz} J., Alcaide D., \enc{Ocaña}{Ocana} J.\cr 
Consumer’s risk in the EMA and FDA regulatory approaches for bioequivalence
in highly variable drugs\cr
Stat Med. 2016;35(12):1933-43. \href{http://dx.doi.org/10.1002/sim.6834}{doi 10.1002/sim.6834}

Labes, D., \enc{Schütz}{Schuetz}, H.\cr
Inflation of Type I Error in the Evaluation of Scaled Average Bioequivalence, 
and a Method for its Control\cr
submitted to Pharm. Res. 2016
}
\author{
H. \enc{Schütz}{Schuetz}
}
\note{
We are doing the sample size estimation only for balanced designs
since the break down of the total subject number in case of unbalanced
sequences is not unique. Moreover the formulas used are only for
balanced designs.\cr
}

\section{Warning }{
The sample size estimation for extreme \code{theta0} (<0.83 or >1.21) may be time
consuming and will eventually also fail since the start values chosen are
not really reasonable in that ranges.\cr
If you really need sample sizes in that range be prepared to restart the sample
size estimation with \code{nstart} above the last one before failure.\cr
Since the dependence of power from n is very flat in the mentioned region you may
also consider to adapt the number of simulations not to tap in the simulation
error trap.\cr\cr
See also the Warning section of the function \code{\link{power.scABEL}} concerning
the power value agreement to those obtained from simulations via subject data.
}

\seealso{
\code{\link{scABEL.ad}}, \code{\link{sampleN.scABEL}}, \code{\link{power.scABEL}},
\code{\link{scABEL}}
}
\examples{
# --- Not run due to timing policy of CRAN for examples
# each may run some ten seconds or more
# using all the defaults:
# partial replicate design, target power=80\%,
# true assumed ratio = 0.90, 1E+6 simulated studies
# EMA regulatory settings (ABE limits, PE constraint 0.8 - 1.25)
\dontrun{
sampleN.scABEL.ad(CV = 0.3)}
# should result in n 60, power 0.8022.
# Note: Without adjustment by sampleN.scABEL(): n 54, power 0.8159
# Easier to show the details:
\dontrun{
sampleN.scABEL.ad(CV = 0.3, details = TRUE)}
#
# full replicate design, target power = 90\%, pre-specified alpha 0.025
\dontrun{
sampleN.scABEL.ad(CV = 0.3, targetpower = 0.8, design = "2x2x4", alpha.pre = 0.025)}
# should result in n 44, power 0.8040; pre-specified alpha justified.
}
