\encoding{utf-8}
\name{power.HVNTID}
\alias{power.HVNTID}
\title{
(Empirical) Power for BE decision via FDA method for highly variable NTIDs
}
\description{
This function performs the power calculation of the BE decision via 
the FDA\enc{’}{'}s method for highly variable narrow therapeutic index drugs (NTIDs) 
as described in respective guidances based on simulations.
The study design could be the full replicate design 2x2x4 with 4-periods (TRTR|RTRT) or the
2x2x3 replicate design with 3-periods and sequences TRT|RTR.
}
\usage{
power.HVNTID(alpha = 0.05, theta1, theta2, theta0, CV, n, design = c("2x2x4", "2x2x3"), 
             nsims = 1e+05, details = FALSE, setseed = TRUE)
}
\arguments{
  \item{alpha}{
Type I error probability, significance level. Commonly set to 0.05.
}
  \item{theta1}{
Conventional lower ABE limit to be applied in the FDA procedure.\cr
Defaults to 0.8 if not given explicitly.
}
  \item{theta2}{
Conventional upper ABE limit to be applied in the FDA procedure.\cr
Defaults to 1.25 if not given explicitly.
}
  \item{theta0}{
\sQuote{True} or assumed T/R ratio. \cr
Defaults to 0.95 if not given explicitly.
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
  \item{n}{
Number of subjects under study.\cr
May be given as vector. In that case it is assumed that \code{n} contains the number
of subjects per sequence groups.\cr
Attention! In case of the \code{"2x2x3"} (TRT|RTR) design the order of sample sizes is important 
if given as vector. \code{n[1]} is for sequence group 'TRT' and \code{n[2]} is for 
sequence group 'RTR'.\cr 
If \code{n} is given as single number (total sample size) and this number is not 
divisible by the number of sequences of the design an unbalanced design is assumed.
A corresponding message is thrown showing the numbers of subjects in the sequence groups.
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
  \item{details}{
If set to \code{TRUE} the computational time is shown as well as the components
for the BE decision.\cr
p(BE-ABE) is the simulated probability for the conventional ABE test.
p(BE-sratio) is the probability that the upper 90\% confidence limit of the 
ratio of sWT/sWR is < 2.5.
}
  \item{setseed}{
Simulations are dependent on the starting point of the (pseudo) random number 
generator. To avoid differences in power for different runs a \code{set.seed(123456)} 
is issued if \code{setseed=TRUE}, the default. 
}
}
\details{
For deciding BE the study must pass the conventional ABE test (90\% CI within the
acceptance range) and additional the test that the ratio of sWT/sWR is < 2.5.\cr\cr
The simulations are done via the distributional properties of the statistical
quantities necessary for deciding BE based on this method.\cr
Details can be found in a document \verb{Implementation_scaledABE_sims} located in 
the \code{/doc} sub-directory of the package.
}
\value{
Returns the value of the (empirical) power if argument \code{details=FALSE}.\cr
Returns a named vector if argument \code{details=TRUE}.\cr
p(BE) is the power, p(BE-ABE) is the power of the ABE test alone and p(BE-sratio)
is the power of the criterion 'ratio of sWT/sWR is <= 2.5' alone.
}
\references{
Food and Drug Administration, Office of Generic Drugs (OGD). \emph{Draft Guidance on Dabigatran Etexilate Mesylate.} Recommended Jun 2012; Revised Sep 2015, Jul 2017. \href{https://www.accessdata.fda.gov/drugsatfda_docs/psg/Dabigatran\%20etexilate\%20mesylate_oral\%20capsule_NDA\%20022512_RV05-17.pdf}{download}

Food and Drug Administration, Office of Generic Drugs (OGD). \emph{Draft Guidance on Rivaroxaban.} Recommended Sep 2015. \href{https://www.accessdata.fda.gov/drugsatfda_docs/psg/Rivaroxaban_oral\%20tablet_22406_RC09-15.pdf}{download}
}
\author{
D. Labes
}
\note{
The FD\enc{’}{'}s guidances recommend only the full replicate design "2x2x4" (TRTR|RTRT).
The results for the design "2x2x3" (TRT|RTR) are to be considered as experimental since
at present not thorougly tested.
}
\seealso{
\code{\link{sampleN.HVNTID}}
and \code{\link{power.NTIDFDA}}, \code{\link{sampleN.NTIDFDA}} for NTIDs with
low variability
}
\examples{
# using the defaults:
# GMR=0.95, theta1=0.8, theta2=1.25, full replicate design 2x2x4, 100,000 simulations
# and a CV of 0.3 (=30\%) for both Reference and Test, with 24 subjects, balanced
power.HVNTID(CV = 0.3, n = 24)
# should give a power of 0.86354
}

