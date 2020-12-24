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
                  CV, design = c("2x3x3", "2x2x4", "2x2x3"), regulator,
                  nstart = NA, nsims = 1e+06, imax = 100, tol, print = TRUE,
                  details = FALSE, alpha.pre = 0.05, setseed = TRUE,
                  sdsims = FALSE, progress)
}
\arguments{
  \item{alpha}{
Type I error (\acronym{TIE}) probability (nominal level of the test). Per
convention commonly set to 0.05.
}
  \item{targetpower}{
Power to achieve at least. Must be >0 and <1. Typical values are 0.80
to 0.90 (\emph{i.e.}, 80\% to 90\%). Defaults to 0.80 if not given explicitly.
}
  \item{theta0}{
\sQuote{True} or assumed T/R ratio. Defaults to 0.90 according to the two \enc{Lászlós}{Laszlos} if not given explicitly.
}
  \item{theta1}{
Conventional lower \acronym{ABE} limit to be applied in the mixed procedure
if \code{CVwR==CVswitch}. Also lower limit for the point estimate
constraint. Defaults to 0.80 if not given explicitly.
}
  \item{theta2}{
Conventional upper \acronym{ABE} limit to be applied in the mixed procedure
if \code{CVwR==CVswitch}. Also upper limit for the point estimate
constraint. Defaults to 1.25 if not given explicitly.
}
  \item{CV}{
Intra-subject coefficient(s) of variation as ratio (not percent).
  \itemize{
    \item If given as a scalar (\code{length(CV)==1}) the \emph{same} \ifelse{latex}{\out{$CV$}}{\ifelse{html}{\out{<em>CV</em>}}{CV}} of Test
          and Reference is assumed (homoscedasticity: \ifelse{latex}{\out{$CV_\textrm{wT}=CV_\textrm{wR}$}}{\ifelse{html}{\out{<em>CV</em><sub>wT</sub>&nbsp=&nbsp;<em>CV</em><sub>wR</sub>}}{CVwT = CVwR}}).
    \item If given as a vector (\code{length(CV)==2}) -- assuming
          heteroscedasticity --\cr
          the \ifelse{latex}{\out{$CV$}}{\ifelse{html}{\out{<em>CV</em>}}{CV}} of Test \strong{must} be given in the \emph{first} element
          and the one of Reference in the \emph{second}.
  }
}
  \item{design}{
Design of the study to be planned.\cr
\code{"2x3x3"} is the partial replicate design.\cr
\code{"2x2x4"} is a full replicate design with 2 sequences and 4 periods.\cr
\code{"2x2x3"} is a full replicate design with 2 sequences and 3 periods.\cr
Defaults to \code{"2x3x3"}. Details are given the section about Designs.
}
  \item{regulator}{
Regulatory settings for the widening of the \acronym{BE} acceptance limits.
Choose from \code{"EMA"} (default), \code{"HC"}, or \code{"GCC"}. This argument may be given
also in lower case.
}
  \item{nstart}{
Best \dQuote{guess} sample size. If not given (default), simulations
start with the sample size estimated for \code{alpha} (or \code{alpha.pre},
if given), \code{theta0}, and \code{targetpower}.\cr
Can also be set to start the sample size search if a previous run
failed.\cr
According to regulatory requirements must be \enc{≥}{>=}12 for the \acronym{EMA} and
\enc{≥}{>=}24 for \acronym{ANVISA}.
}
  \item{nsims}{
Number of simulations to be performed to estimate the (empirical)
\acronym{TIE} and in each iteration of adjusting alpha. The default value
1,000,000 = 1E+6 should not be lowered.
}
  \item{imax}{
Maximum number of steps in sample size search. Defaults to 100.
}
  \item{tol}{
Desired accuracy (convergence tolerance). Defaults to 1E-6.
}
  \item{print}{
If \code{TRUE} (default), the function sends its results to the console.
}
  \item{details}{
If \code{TRUE} (default), the steps during sample size search are
shown. Additionally information about the impact on power by
adjusting alpha and change of study costs due to the increased sample
size is given.
}
  \item{alpha.pre}{
Pre-specified alpha (optional). Must be \code{<=alpha}. \acronym{ABEL} will be
performed at level \verb{alpha.pre} and the \acronym{TIE} assessed at level \code{alpha}.\cr
Less powerful than adjusting alpha but an alternative in the critical
region of maximum inflation of the \acronym{TIE}. In certain scenarios
Bonferroni\enc{’}{'}s 0.025 is not sufficient to preserve the Type I Error.\cr
Not recommended if \ifelse{latex}{\out{$CV_\textrm{wR}\geq0.45$}}{\ifelse{html}{\out{<em>CV</em><sub>wR</sub>&nbsp;&ge;&nbsp;0.45}}{CVwR >= 0.45}} due to poor power characteristics.
}
  \item{setseed}{
Simulations are dependent on the starting point of the (pseudo)
random number generator. To avoid differences in power for different
runs \code{set.seed(123456)} is issued if \code{setseed=TRUE} (default).
}
  \item{sdsims}{
If \code{FALSE} (default) power is estimated by the respective \sQuote{key} statistics.
Recommended for speed reasons.\cr
Set to \code{TRUE} if results of \code{\link{power.scABEL}} are expected to
be inaccurate (partial replicate design with unbalanced sequences and/or heteroscedasticity
where \ifelse{latex}{\out{$CV_\textrm{wT}>CV_\textrm{wR}$}}{\ifelse{html}{\out{<em>CV</em><sub>wT</sub>&nbsp;&gt;&nbsp;<em>CV</em><sub>wR</sub>}}{CVwT > CVwR}}) and subject data via \code{\link{power.scABEL.sdsims}} should
be simulated instead. Very time consuming (easily 100times slower)! Subject data
simulations are only supported for \code{regulator="EMA"} and \code{regulator="GCC"}.
}
  \item{progress}{
Set to \code{TRUE} if a progress bar should be displayed. Defaults to \code{FALSE}.\cr
Ignored if \code{sdsims=FALSE}.
}
}
\details{
The simulations are done via the distributional properties of the statistical
quantities necessary for assessing \acronym{BE} based on \acronym{ABEL}.
Simulations of the \acronym{TIE} are performed at the upper (expanded) limit \emph{U}
of the acceptance range. Due to the symmetry around 1 results are valid for the lower
(expanded) limit \emph{L} as well.\cr
\emph{U} at the \acronym{EMA}\enc{’}{'}s and Health Canada\enc{’}{'}s \code{CVcap}, the \acronym{GCC}\enc{’}{'}s for \emph{any} \ifelse{latex}{\out{$CV_\textrm{wR}>0.30$}}{\ifelse{html}{\out{<em>CV</em><sub>wR</sub>&nbsp;&gt;&nbsp;0.30}}{CVwR > 0.30}}:\preformatted{
scABEL(CV=0.5, reg="EMA")[["upper"]]
[1] 1.43191
scABEL(CV=0.57382, reg="HC")[["upper"]]
[1] 1.5
scABEL(CV=0.31, reg="GCC")[["upper"]]
[1] 1.333333}
Simulated studies are evaluated by ANOVA (Method A) as recommended in the
\acronym{EMA}\enc{’}{'}s Q&A-document and by intra-subject contrasts if \code{regulator="HC"}.
Health Canada requires a mixed effects model which cannot be implemented in R. However,
intra-subjects contrasts are a sufficiently close approximation.\cr

If an inflation of the \acronym{TIE} is expected (\emph{i.e.}, \code{>alpha}), alpha is
iteratively adjusted until at least the target power is reached and the consumer
risk is maintained (\code{<=alpha}). For details about the algorithm see the
respective section of \code{\link{scABEL.ad}}.\cr\cr
The estimated sample size gives always the \emph{total} number of subjects (not subject/sequence -- like in some other software packages).
}
\value{
Returns a data.frame with the input and results for adjusted alpha,
type I error, sample size, and achieved power.\cr
The \code{Sample size} column contains the total sample size.
If no adjustment is necessary, \code{NA} will be returned in the
\code{alpha.adj} column and other results are identical to the ones
obtained by \code{\link{sampleN.scABEL}}.
}
\references{
\enc{Tóthfalusi}{Tothfalusi} L, \enc{Endrényi}{Endrenyi} L. \emph{Sample Sizes for Designing Bioequivalence Studies for Highly Variable Drugs.} J Pharm Pharmaceut Sci. 2011;15(1):73--84. \href{http://ejournals.library.ualberta.ca/index.php/JPPS/article/download/11612/9489}{open access}

Wonnemann M, \enc{Frömke}{Froemke} C, Koch A. \emph{Inflation of the Type I Error: Investigations on Regulatory Recommendations for Bioequivalence of Highly Variable Drugs.} Pharm Res. 2015;32(1):135--43. \doi{10.1007/s11095-014-1450-z}

\enc{Muñoz}{Munoz} J, Alcaide D, \enc{Ocaña}{Ocana} J. \emph{Consumer’s risk in the EMA and FDA regulatory approaches for bioequivalence in highly variable drugs.} Stat Med. 2015;35(12):1933--43. \doi{10.1002/sim.6834}

Labes D, \enc{Schütz}{Schuetz} H. \emph{Inflation of Type I Error in the Evaluation of Scaled Average Bioequivalence, 
and a Method for its Control.} Pharm Res. 2016;33(11):2805--14. \doi{10.1007/s11095-016-2006-1}
}
\author{
H. \enc{Schütz}{Schuetz}
}
\note{
We are doing the sample size estimation only for balanced designs
since the break down of the total subject number in case of unbalanced
sequences is not unique. Moreover the formulas used are only for
balanced designs.
}
\section{Designs}{
Although some designs are more \sQuote{popular} than others, sample size estimations are valid for \emph{all} of the following designs:
  \tabular{ll}{
    \code{"2x2x4"} \tab TRTR | RTRT\cr
    \tab TRRT | RTTR\cr
    \tab TTRR | RRTT\cr
    \code{"2x2x3"} \tab TRT | RTR\cr
    \tab TRR | RTT\cr
    \code{"2x3x3"} \tab TRR | RTR | RRT
  }
}
\section{Warning}{
The sample size estimation for extreme \code{theta0} (<0.83 or >1.21) may be time
consuming and will eventually also fail since the start values chosen are
not really reasonable in that ranges.\cr
If you really need sample sizes in that range be prepared to restart the sample
size estimation with \code{nstart} above the last one before failure.\cr
Since the dependence of power from \code{n} is very flat in the mentioned region you may
also consider to adapt the number of simulations not to tap in the simulation
error trap.\cr\cr
See also the Warning section of the function \code{\link{power.scABEL}} concerning
the power value agreement to those obtained from simulations via subject data.\cr\cr
For the \acronym{GCC} and \ifelse{latex}{\out{$CV_\textrm{wR}\leq0.30$}}{\ifelse{html}{\out{<em>CV</em><sub>wR</sub>&nbsp;&le;&nbsp;0.30}}{CVwR <= 0.30} simulations will be time consuming and may result in large sample sizes.
}
}

\seealso{
\code{\link{scABEL.ad}}, \code{\link{sampleN.scABEL}}, \code{\link{power.scABEL}},
\code{\link{scABEL}}
}
\examples{
# --- Not run due to timing policy of CRAN for examples
# each may run some ten seconds or more
# using all the defaults:
# TRR|RTR|RRT, target power 80\%, assumed ratio 0.90, 1E+6 simulated studies,
# EMA regulatory settings (ABE limits, PE constraint 0.8 - 1.25)
\donttest{
sampleN.scABEL.ad(CV = 0.3)}
# should result in n 60, power 0.8022.
# Note: Without adjustment by sampleN.scABEL(): n 54, power 0.8159
# Easier to show the details:
\donttest{
sampleN.scABEL.ad(CV = 0.3, details = TRUE)}
#
# TRTR|RTRT, target power 90\%, pre-specified alpha 0.025
\donttest{
sampleN.scABEL.ad(CV = 0.3, targetpower = 0.9, design = "2x2x4", alpha.pre = 0.025)}
# should result in n 60, power 0.9021; pre-specified alpha justified.
}
