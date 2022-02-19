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
scABEL.ad(alpha = 0.05, theta0, theta1, theta2, CV,
          design = c("2x3x3", "2x2x4", "2x2x3"), regulator,
          n, alpha.pre = 0.05, imax = 100, tol, print = TRUE,
          details = FALSE, setseed = TRUE, nsims = 1e+06,
          sdsims = FALSE, progress)
}
\arguments{
  \item{alpha}{
Type I Error (\acronym{TIE}) probability (nominal level of the test). Per
convention commonly set to 0.05.
}
  \item{theta0}{
\sQuote{True} or assumed T/R ratio. Defaults to 0.90 according to the two \enc{Lászlós}{Laszlos} if not given explicitly.
}
  \item{theta1}{
Conventional lower \acronym{ABE} limit to be applied in the mixed procedure
if \code{CVwR == CVswitch}. Also lower limit for the point estimate
constraint. Defaults to 0.80 if not given explicitly.
}
  \item{theta2}{
Conventional upper \acronym{ABE} limit to be applied in the mixed procedure
if \code{CVwR == CVswitch}. Also upper limit for the point estimate
constraint. Defaults to 1.25 if not given explicitly.
}
  \item{CV}{
Intra-subject coefficient(s) of variation as ratio (not percent).
  \itemize{
    \item If given as a scalar (\code{length(CV) == 1}) the \emph{same} CV of Test
          and Reference is assumed (homoscedasticity, \code{CVwT==CVwR}).
    \item If given as a vector (\code{length(CV) == 2}), \emph{i.e.}, assuming
          heteroscedasticity, the CV of the Test \strong{must} be given in \code{CV[1]} and the one of the Reference in the \code{CV[2]}.
  }
}
  \item{design}{
Design of the study.\cr
\code{"2x3x3"} is the partial replicate design.\cr
\code{"2x2x4"} is a full replicate design with 2 sequences and 4 periods.\cr
\code{"2x2x3"} is a full replicate design with 2 sequences and 3 periods.\cr
Defaults to \code{"2x3x3"}. Details are given the section about Designs.
}
  \item{regulator}{
Regulatory settings for the expanding of the \acronym{BE} acceptance limits.
Choose from \code{"EMA"} (default), \code{"HC"}, \code{"GCC"}, or \code{"FDA"}. This argument
may also be given in lower case.
}
  \item{n}{
Total sample size of the study or a vector of sample size / sequences.
If \code{n} leads to an unbalanced design (\emph{i.e.}, is not a multiple of two
in the full replicate designs or not a multiple of three in the partial
replicate), the code tries to keep subjects / sequence as balanced as possible.\cr
In evaluating a particular \emph{unbalanced} study \strong{always} give \code{n}
as a vector.\cr
Only if \code{design = "2x2x3"} (TRT|RTR) the \emph{order} of sample sizes
is important. \code{n[1]} is for sequence TRT and \code{n[2]} for sequence RTR.\cr
If \code{n} is missing, a sample size is estimated with target power 0.80 and pre-specified
alpha if defined. Otherwise, alpha is used.
}
  \item{alpha.pre}{
Pre-specified alpha (optional). Must be \code{<=alpha}. \acronym{ABEL} will be
performed at level \verb{alpha.pre} and the \acronym{TIE} assessed at level \verb{alpha}.\cr
Less powerful than adjusting alpha but an alternative in the critical region
of maximum inflation of the \acronym{TIE}. In certain scenarios \enc{Bonferroni’s}{Bonferroni's}
0.025 is not sufficient to preserve the Type I Error (\emph{e.g.}, the third example).\cr
Not recommended if \code{CVwR >= 0.45} due to poor power characteristics.
}
  \item{imax}{
Maximum number of steps in sample size search. Defaults to 100.
}
  \item{tol}{
Desired accuracy (convergence tolerance). Defaults to 1E-6.
}
  \item{print}{
If \verb{TRUE} (default), the function sends its results to the console.
}
  \item{details}{
If \code{TRUE}, the \emph{relative} change of the consumer risk in percent is shown.
Additionally information about the impact on power (for specified \code{theta0}
and target power 0.80), runtime, and number of simulations (iterations)
are given. Defaults to \code{FALSE}.
}
  \item{setseed}{
Simulations are dependent on the starting point of the (pseudo)
random number generator. To avoid differences in power for different
runs \code{set.seed(123456)} is issued if \code{setseed=TRUE} (default).
}
  \item{nsims}{
Number of simulations to be performed to estimate the (empirical)
\acronym{TIE} error and in each iteration of adjusting alpha. The default value
1,000,000 = 1E+6 should not be lowered.
}
  \item{sdsims}{
If \code{FALSE} (default) power is estimated by the respective \sQuote{key} statistics.
Recommended for speed reasons.\cr
Set to \code{TRUE} if results of \code{\link{power.scABEL}} are expected to
be inaccurate (partial replicate design with unbalanced sequences and/or heteroscedasticity
where CVwT > CVwR) and subject data via \code{\link{power.scABEL.sdsims}} should
be simulated instead. Very time consuming (easily 100times slower)! Subject data
simulations are only supported for \code{regulator = "EMA"} and \code{regulator = "GCC"}.
}
  \item{progress}{
Set to \code{TRUE} if a progress bar should be displayed. Ignored if \code{sdsims = FALSE}.
}
}
\details{
The simulations are done via the distributional properties of the statistical
quantities necessary for assessing \acronym{BE} based on \acronym{ABEL}.
Simulations for the \acronym{TIE} are performed at the upper (expanded) limit \emph{U}
of the acceptance range. Due to the symmetry around 1 results are valid for the lower
(expanded) limit \emph{L} as well.\cr
\emph{U} at the \acronym{EMA}\enc{’}{'}s and Health Canada\enc{’}{'}s \code{CVcap}, the \acronym{GCC}\enc{’}{'}s for \emph{any} CVwR > 0.3:\preformatted{scABEL(CV = 0.5, reg = "EMA")[["upper"]]
[1] 1.43191
scABEL(CV = 0.57382, reg = "HC")[["upper"]]
[1] 1.5
scABEL(CV = 0.5, reg = "GCC")[["upper"]]
[1] 1.333333}
Simulated studies are evaluated by ANOVA (Method A) as recommended in the
\acronym{EMA}\enc{’}{'} Q&A-document and by intra-subject contrasts if \code{regulator = "HC"}.
Health Canada requires a mixed-effects model which cannot be implemented in R.
However, intra-subjects contrasts are a sufficiently close approximation.\cr
The Type I Error in \acronym{ABEL} depends only on \code{CVwR} and -- to a
minor degree -- the sample size. Algorithm:
\enumerate{
  \item The \acronym{TIE} is assessed based on \code{alpha} (or \code{alpha.pre})
  and compared to the nominal level of the test \code{alpha}.
  \item If no inflation of the \acronym{TIE} is found, the algorithm stops.
  \item Otherwise, alpha is iteratively adjusted (\emph{i.e.}, \code{alpha.adj <alpha})
  until no more relevant inflation of the \acronym{TIE} is detected (\emph{i.e.},
  \code{abs(TIE - alpha) <= tol}).
}
}
\value{
Sends results to the console if argument \code{print=TRUE} (default).\cr
Returns a list with the input, adjusted alpha, and Type I Error (for nominal
and adjusted alpha) if argument \code{print=FALSE}.\cr
If no adjustment is necessary, \code{NAs} will be returned for the respective
variables (\code{alpha.adj}, \code{TIE.adj}, \code{rel.change}, \code{pwr.adj}, \code{rel.loss}).
}
\references{
Wonnemann M, \enc{Frömke}{Froemke} C, Koch A. \emph{Inflation of the Type I Error: Investigations on Regulatory Recommendations for Bioequivalence of Highly Variable Drugs.} Pharm Res. 2015;32(1):135--43. \doi{10.1007/s11095-014-1450-z}

\enc{Muñoz}{Munoz} J, Alcaide D, \enc{Ocaña}{Ocana} J. \emph{Consumer’s risk in the EMA and FDA regulatory approaches for bioequivalence in highly variable drugs.} Stat Med. 2015;35(12):1933--43. \doi{10.1002/sim.6834}

Labes D, \enc{Schütz}{Schuetz} H. \emph{Inflation of Type I Error in the Evaluation of Scaled Average Bioequivalence, and a Method for its Control.} Pharm Res. 2016;33(11):2805--14. \doi{10.1007/s11095-016-2006-1}

\enc{Tóthfalusi}{Tothfalusi} L, \enc{Endrényi}{Endrenyi} L. \emph{Algorithms for Evaluating Reference Scaled Average Bioequivalence: Power, Bias, and Consumer Risk.} Stat Med. 2017;36(27):4378--90. \doi{10.1002/sim.7440}

Molins E, Cobo E, \enc{Ocaña}{Ocana} J. \emph{Two-Stage Designs Versus European Scaled Average Designs in Bioequivalence Studies for Highly Variable Drugs: Which to Choose?} Stat Med. 2017;36(30):4777--88. \doi{10.1002/sim.7452}

European Medicines Agency, Committee for Medicinal Products for Human Use. \emph{Guideline on the Investigation of Bioequivalence.} London, 20 January 2010. \href{https://www.ema.europa.eu/en/documents/scientific-guideline/guideline-investigation-bioequivalence-rev1_en.pdf}{CPMP/EWP/QWP/1401/98 Rev. 1/ Corr **}

European Medicines Agency, Committee for Medicinal Products for Human Use. \emph{Questions & Answers: positions on specific questions addressed to the Pharmacokinetics Working Party (PKWP).} London, 19 November 2015. \href{https://www.ema.europa.eu/en/documents/scientific-guideline/questions-answers-positions-specific-questions-addressed-pharmacokinetics-working-party_en.pdf}{EMA/618604/2008 Rev. 13}

Health Canada, Therapeutic Products Directorate. \emph{Comparative Bioavailability Standards: Formulations Used for Systemic Effects, 2.1.1.8 Highly variable drug products} Ottawa, 08 June 2018. \href{https://www.canada.ca/en/health-canada/services/drugs-health-products/drug-products/applications-submissions/guidance-documents/bioavailability-bioequivalence/comparative-bioavailability-standards-formulations-used-systemic-effects.html#a2.1.1.8}{online}

Executive Board of the Health Ministers’ Council for GCC States. \emph{The GCC Guidelines for Bioequivalence.} May 2021. \href{https://www.sfda.gov.sa/sites/default/files/2021-10/GCC_Guidelines_Bioequivalence.pdf}{Version 3.0}
}

\author{
H. \enc{Schütz}{Schuetz}
}
\note{
Specifying \code{theta0} is not necessary.\cr
If \code{theta0} \emph{is not} given, achievable power for the common target
of 0.80 (both for \code{alpha} and adjusted alpha) will be estimated. If
\code{theta0} \emph{is} specified, its value will be used; again for target power 0.80.\cr
If you are interested in other levels of power, use \code{\link{sampleN.scABEL.ad}}.

The EMA\enc{’s}{'s} method is currently recommended in other jurisdictions as well (\emph{e.g.}, by the WHO;
in ASEAN States, Australia, Brazil, Egypt, the Eurasian Economic Union, New Zealand, and the East African Community).\cr
If CVwR > 30\%, \emph{fixed} wider limits of 0.7500--1.3333 are recommended by the Gulf Cooperation Council (Bahrain, Kuwait, Oman, Qatar, Saudi Arabia, United Arab Emirates).
}
\section{Designs}{
Although some designs are more \sQuote{popular} than others, power calculations are valid for \emph{all} of the following designs:
  \tabular{ll}{
    \code{"2x2x4"} \tab TRTR | RTRT\cr
    \tab TRRT | RTTR\cr
    \tab TTRR | RRTT\cr
    \code{"2x2x3"} \tab TRT | RTR\cr
    \tab TRR | RTT\cr
    \code{"2x3x3"} \tab TRR | RTR | RRT
  }
}
\section{Warning }{
See the Warning section of the function \code{\link{power.scABEL}} concerning
the power value agreement to the one obtained by simulations via subject data.
}

\seealso{
\code{\link{sampleN.scABEL.ad}}, \code{\link{power.scABEL}}, \code{\link{power.scABEL.sdsims}}, \code{\link{scABEL}}
}
\examples{
# Using all defaults:
# TRR|RTR|RRT, target power 80\% for assumed ratio 0.90 (estimated sample size 54),
# EMA regulatory settings (ABE limits and PE constraint 0.80 - 1.25),
# 1E+6 simulated studies.
# Not run: due to timing policy of CRAN for examples
\donttest{
scABEL.ad(CV = 0.3)}
# Should result in adjusted alpha 0.03389 (TIE 0.5000, TIE for nominal alpha 0.07189).
#
# As above but subject data simulations.
\donttest{
scABEL.ad(CV = 0.3, sdsims = TRUE)}
# Should result in adjusted alpha 0.03336 (TIE 0.5000, TIE for nominal alpha 0.07237).
#
# TRT|RTR, heteroscedasticity, sample size 48 (unbalanced), subject data simulations.
\donttest{
scABEL.ad(CV = c(0.25, 0.3), design = "2x2x3", n = c(23, 25), sdsims = TRUE)}
# Should result in adjusted alpha 0.02465 (TIE 0.5000, TIE for nominal alpha 0.09050).
#
# TRTR|RTRT, CV 0.35, sample size 33 (unbalanced).
\donttest{
scABEL.ad(CV = 0.35, design = "2x2x4", n = c(16, 17))}
# Should result in adjusted alpha 0.03632 (TIE 0.5000, TIE for nominal alpha 0.06544).
}
