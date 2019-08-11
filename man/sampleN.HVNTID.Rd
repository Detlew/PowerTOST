\encoding{utf-8}
\name{sampleN.HVNTID}
\alias{sampleN.HVNTID}
\title{
Sample size estimation for BE decision via FDA method for highly variable (HV) 
narrow therapeutic index drugs (NTIDs)
}
\description{
This function performs the sample size estimation for the BE decision via 
the FDA’s method for highly variable NTIDs as described in respective guidances based on simulations. \cr
The study designs may be the full replicate design 2x2x4 with 4 periods (TRTR|RTRT) and 
the 3-period replicate design 2x2x3 with sequences RTR|TRT.
}
\usage{
sampleN.HVNTID(alpha = 0.05, targetpower = 0.8, theta0, theta1, theta2, CV, 
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
Defaults to 0.95 if not given explicitly.
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
Coefficient(s) of variation as ratio.\cr
If \code{length(CV) = 1} the same CV is assumed for Test and Reference.\cr
If \code{length(CV) = 2} the CV for Test must be given in \code{CV[1]} and for 
Reference in \code{CV[2]}.
}
  \item{design}{
Design of the study to be planned.\cr
\code{"2x2x4"} is the full replicate design with 2 sequences and 4 periods (TRTR|RTRT).\cr
\code{"2x2x3"} is the 3-period ful replicate design with sequences TRT|RTR.\cr
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
If \verb{TRUE} (default) the function prints its results. If \verb{FALSE} only the resulting dataframe will be returned. 
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
For deciding BE the study must pass the conventional ABE test and additionally
the test that the ratio of sWT/sWR is <= 2.5.\cr\cr
The simulations are done via the distributional properties of the statistical
quantities necessary for deciding BE based on these method.\cr
Details can be found in a document "Implementation_scaledABE_sims" located in 
the \code{\\doc} sub-directory of the package.
}
\value{
Returns a data.frame with the input and sample size results.\cr
The \code{Sample size} column contains the total sample size.\cr
The \code{nlast} column contains the last \code{n} value. May be useful for re-starting.
}
\references{
Food and Drug Administration, Office of Generic Drugs (OGD). \emph{Draft Guidance on Dabigatran Etexilate Mesylate.} Recommended Jun 2012; Revised Sep 2015, Jul 2017. \href{https://www.accessdata.fda.gov/drugsatfda_docs/psg/Dabigatran\%20etexilate\%20mesylate_oral\%20capsule_NDA\%20022512_RV05-17.pdf}{download}

Food and Drug Administration, Office of Generic Drugs (OGD). \emph{Draft Guidance on Rivaroxaban.} Recommended Sep 2015. \href{https://www.accessdata.fda.gov/drugsatfda_docs/psg/Rivaroxaban_oral\%20tablet_22406_RC09-15.pdf}{download}
}
\author{
D. Labes
}
\note{
The design recommended by the FDA is the full replicate design 2x2x4.\cr 
The sample size estimation is done only for balanced studies since the 
break down of the total subject number in case of unbalanced sequence groups
is not unique. Moreover the formulas used are only valid for balanced designs.\cr
The FDA method is described for the ABE limits 0.8 ... 1.25 only. Setting theta1, theta2
to other values may not be reasonable and is not tested.\cr
The minimum sample size is n=6, even if the power is higher than the intended 
targetpower.
}
\section{Warning }{
For some input constellations the sample size search may be very time 
consuming and will eventually also fail since the start values chosen may
not really reasonable for them.\cr
In case of a failed sample size search you may restart with setting the argument 
\verb{nstart}.\cr
}

\seealso{ 
\code{\link{power.HVNTID}} \cr 
and \code{\link{power.NTIDFDA}}, \code{\link{sampleN.NTIDFDA}} for NTIDs with
low variability
}
\examples{
# using all defaults but CV
sampleN.HVNTID(CV = 0.3)
# should give
# n=22 with an (empirical) power of 0.829700

# Test formulation with lower variability but same pooled CV
CVs <- CVp2CV(0.3, ratio = 0.25)
sampleN.HVNTID(CV = CVs)
# should give (no distinct difference to example above)
# n=22 with an (empirical) power of 0.837520
}

