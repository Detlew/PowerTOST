\encoding{utf-8}
\name{sampleN.NTIDFDA}
\alias{sampleN.NTIDFDA}
\title{
Sample size estimation for BE decision via the FDA's method for narrow therapeutic index drugs (NTIDs)
}
\description{
This function performs the sample size estimation for the BE decision for 
the FDA\enc{’}{'}s method for NTIDs based on simulations. The study design is the full 
replicate design 2x2x4 (TRTR|RTRT) or the 3-period replicate design with sequences TRT|RTR.
}
\usage{
sampleN.NTIDFDA(alpha = 0.05, targetpower = 0.8, theta0, theta1, theta2, CV, 
                design = c("2x2x4", "2x2x3"), nsims = 1e+05, nstart, imax = 100,
                print = TRUE, details = TRUE, setseed = TRUE)
}
\arguments{
  \item{alpha}{
Type I error probability. Per convention mostly set to 0.05.
}
  \item{targetpower}{
Power to achieve at least. Must be >0 and <1.\cr
Typical values are 0.8 or 0.9.
}
  \item{theta0}{
\sQuote{True} or assumed T/R ratio.\cr
Attention! Defaults here to 0.975 if not given explicitly. The value was chosen 
closer to 1 because the potency (contents) settings for NTIDs are tightened 
by the FDA.
}
  \item{theta1}{
Conventional lower ABE limit to be applied in the FDA procedure.\cr
Defaults to 0.8 if not given explicitly.
}
  \item{theta2}{
Conventional upper ABE limit to be applied in the FDA procedure.\cr
Defaults to 1.25 if not given explicitly.
}
  \item{CV}{
Intra-subject coefficient(s) of variation as ratio (not percent).
  \itemize{
    \item If given as a scalar (\code{length(CV)==1}) the \emph{same} CV of Test
          and Reference is assumed (homoscedasticity, \code{CVwT==CVwR}).
    \item If given as a vector (\code{length(CV)==2}), \emph{i.e.}, assuming
          heteroscedasticity, the CV of the Test \strong{must} be given in \code{CV[1]} and the one of the Reference in the \code{CV[2]}.
  }
}
  \item{design}{
Design of the study to be planned.\cr
\code{"2x2x4"} is the full replicate with 2 sequences and 4 periods (TRTR|RTRT).\cr
\code{"2x2x3"} is the full replicate with 2 sequences and 3 periods (TRT|RTR).\cr
Defaults to \code{design="2x2x4"}.
}
  \item{nsims}{
Number of simulations to be performed to obtain the empirical power.
Defaults to 100,000 = 1e+5.
}
  \item{nstart}{
Set this to a start value for the sample size if a previous run failed.\cr
May be missing.
}
  \item{imax}{
Maximum number of steps in sample size search. Defaults to 100.
}
  \item{print}{
If \code{TRUE} (default) the function prints its results. If \code{FALSE} only the resulting dataframe will be returned.
}
  \item{details}{
If set to \code{TRUE}, the default, the steps during sample size search are shown.
Moreover the details of the method settings are printed.
}
  \item{setseed}{
Simulations are dependent on the starting point of the (pseudo) random number 
generator. To avoid differences in power values for different runs a 
\code{set.seed(123456)} is issued if \code{setseed=TRUE}, the default.
}
}
\details{
The linearized scaled ABE criterion is calculated according to the SAS code
given in the FDA Warfarine guidance. For deciding BE the study must pass that criterion, 
the conventional ABE test and additional the test that the ratio of 
\code{sWT/sWR <= 2.5.}\cr\cr
The simulations are done via the distributional properties of the statistical
quantities necessary for deciding BE based on these method.\cr
Details can be found in a document \verb{Implementation_scaledABE_sims} located in 
the \code{/doc} sub-directory of the package.\cr\cr
The estimated sample size gives always the \emph{total} number of subjects (not subject/sequence -- like in some other software packages).
}
\value{
Returns a data.frame with the input settings and sample size results.\cr
The \code{Sample size} column contains the total sample size.\cr
The \code{nlast} column contains the last \code{n} value. May be useful for re-starting.
}
\references{
Food and Drug Administration, Office of Generic Drugs (OGD). \emph{Draft Guidance on Warfarin Sodium.} Recommended Dec 2012. \href{https://www.accessdata.fda.gov/drugsatfda_docs/psg/Warfarin_Sodium_tab_09218_RC12-12.pdf}{download}

Yu LX, Jiang W, Zhang X, Lionberger R, Makhlouf F, Schuirmann DJ, Muldowney L, Chen ML, Davit B, Conner D, Woodcock J. \emph{Novel bioequivalence approach for narrow therapeutic index drugs.} Clin Pharmacol Ther. 2015;97(3):286--91. \doi{10.1002/cpt.28}

Jiang W, Makhlouf F, Schuirmann DJ, Zhang X, Zheng N, Conner D, Yu LX, Lionberger R. \emph{A Bioequivalence Approach for Generic Narrow Therapeutic Index Drugs: Evaluation of the Reference-Scaled Approach and Variability Comparison Criterion.} AAPS J. 2015;17(4):891--901. \doi{10.1208/s12248-015-9753-5}

\enc{Endrényi}{Endrenyi} L, \enc{Tóthfalusi}{Tothfalusi} L. \emph{Determination of Bioequivalence for Drugs with Narrow Therapeutic Index: Reduction of the Regulatory Burden.} J Pharm Pharm Sci. 2013;16(5):676--82. \href{https://journals.library.ualberta.ca/jpps/index.php/JPPS/article/download/20900/15927/0}{open access}
}
\author{
D. Labes
}
\note{
The design recommended by the FDA is the full replicate design \code{"2x2x4"}.\cr 
The sample size estimation is done only for balanced studies since the 
break down of the total subject number in case of unbalanced sequence groups
is not unique. Moreover the formulas used are only valid for balanced designs.\cr
The FDA method is described for the ABE limits 0.8 ... 1.25 only. Setting theta1, 
theta2 to other values may not be reasonable and is not tested.\cr\cr
The results for the design \code{"2x2x3"} are to be considered as experimental since
at present not thorougly tested.\cr
The minimum sample size is 6, even if the power is higher than the intended 
targetpower.
}
\section{Warning}{
For some input constellations the sample size search may be very time 
consuming and will eventually also fail since the start values chosen may
not really reasonable for them. This applies especially for theta0 values
near to the implied scaled (tightened/widened) ABE limits according to
\code{exp(±1.053605*swR)}.\cr
In case of a failed sample size search you may restart with setting the argument 
\verb{nstart}.\cr
In case of theta0 values outside the implied scaled (tightened/widened) ABE limits
no sample size estimation is possible and the function throws an error 
(f.i. CV=0.04, theta0=0.95).
}

\seealso{ 
\code{\link{power.NTIDFDA}}
and \code{\link{power.HVNTID}}, \code{\link{sampleN.HVNTID}} for NTIDs with
high variability
}
\examples{
sampleN.NTIDFDA(CV = 0.04,theta0 = 0.975)
# should give
# n=54 with an (empirical) power of 0.809590
#
# Test formulation with lower variability
sampleN.NTIDFDA(CV = c(0.04,0.06),theta0 = 0.975)
# should give
# n=20 with an (empirical) power of 0.0.814610
#
# alternative 3-period design
sampleN.NTIDFDA(CV = 0.04,theta0 = 0.975, design="2x2x3")
# should give
# n=86 with power = 0.80364
}

