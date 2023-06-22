\encoding{utf-8}
\name{power.NTID}
\alias{power.NTID}
\alias{power.NTIDFDA}
\title{
(Empirical) Power for BE decision via FDA method for NTIDs
}
\description{
This function performs the power calculation of the BE decision via
the FDA\enc{’}{'}s method for narrow therapeutic index drugs (NTIDs) by simulations.
The study design could be the full replicate design 2x2x4 with 4-periods (TRTR|RTRT) or the
2x2x3 replicate design with sequences TRT|RTR.
}
\usage{
power.NTID(alpha = 0.05, theta1, theta2, theta0, CV, n, design=c("2x2x4", "2x2x3"),
           nsims = 1e+05, details = FALSE, setseed = TRUE)
}
\arguments{
  \item{alpha}{
Type I error probability, significance level. Conventionally mostly set to 0.05.
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
Attention! Defaults here to 0.975 if not given explicitly. The value was chosen
closer to 1 because the potency (contents) settings for NTIDs are tightened
by the FDA.
}
  \item{CV}{
Intra-subject coefficient(s) of variation as ratio (not percent).
  \itemize{
    \item If given as a scalar (\code{length(CV) == 1}) the \emph{same} CV of Test
          and Reference is assumed (homoscedasticity, \code{CVwT == CVwR}).
    \item If given as a vector (\code{length(CV) == 2}), \emph{i.e.}, assuming
          heteroscedasticity, the CV of the Test \strong{must} be given in \code{CV[1]} and the one of the Reference in the \code{CV[2]}.
  }
}
  \item{n}{
Number of subjects under study.\cr
May be given as vector. In that case it is assumed that n contains the number
of subjects per sequence groups.\cr
Attention! In case of the \code{"2x2x3"} (TRT|RTR) design the order of sample sizes important
if given as vector. \code{n[1]} is for sequence group \sQuote{TRT} and \code{n[2]} is for
sequence group \sQuote{RTR}.\cr\cr
If \code{n} is given as single number (total sample size) and this number is not
divisible by the number of sequences of the design an unbalanced design is assumed.
A corresponding message is thrown showing the numbers of subjects in the sequence groups.
}
  \item{design}{
Design of the study to be planned.\cr
\code{"2x2x4"} is the full replicate design with 2 sequences and 4 periods (TRTR|RTRT).\cr
\code{"2x2x3"} is the full replicate design with 2 sequences and 3 periods (TRT|RTR).\cr
Defaults to \code{design="2x2x4"}.
}
  \item{nsims}{
Number of simulations to be performed to obtain the empirical power.
Defaults to 100,000 = 1e+5.
}
  \item{details}{
If set to \code{TRUE} the computational time is shown as well as the components
for the BE decision.\cr
\verb{p(BE-ABE)} is the simulated probability for the conventional ABE test.\cr
\verb{p(BE-sABEc)} is the probability that the 95\% CI of the scaled ABE criterion is \ifelse{latex}{\out{$\leq0$}}{\ifelse{html}{\out{&leq;&nbsp;0}}{<= 0}}.\cr
\verb{p(BE-sratio)} is the probability that the upper confidence limit of
\ifelse{latex}{\out{$\sigma_{\rm{wT}}/\sigma_{\rm{wR}}\leq 2.5$}}{\ifelse{html}{\out{<em>&sigma;</em><sub>wT</sub>/<em>&sigma;</em><sub>wR</sub>&nbsp;&leq;&nbsp;2.5}}{sWT/sWR <= 2.5}}.
}
  \item{setseed}{
Simulations are dependent on the starting point of the (pseudo) random number
generator. To avoid differences in power for different runs a \code{set.seed(123456)}
is issued if \code{setseed = TRUE}, the default.
}
}
\details{
The linearized scaled ABE criterion is calculated according to the SAS code
given in the FDA\enc{’}{'}s guidances. For deciding BE the study must pass that criterion,
the conventional ABE test, and that the upper confidence limit of
\ifelse{latex}{\out{$\sigma_{\rm{wT}}/\sigma_{\rm{wR}}\leq 2.5$}}{\ifelse{html}{\out{<em>&sigma;</em><sub>wT</sub>/<em>&sigma;</em><sub>wR</sub>&nbsp;&leq;&nbsp;2.5}}{sWT/sWR <= 2.5}}.\cr

The simulations are done via the distributional properties of the statistical
quantities necessary for deciding BE based on these method.\cr
Details can be found in a document \verb{Implementation_scaledABE_sims} located in
the \code{/doc} sub-directory of the package.
}
\value{
Returns the value of the (empirical) power if argument \code{details = FALSE}.\cr
Returns a named vector if argument \code{details = TRUE}, where
\verb{p(BE)} is the (overall) power, \verb{p(BE-sABEc)} is the power of the BE test via scaled ABE criterion alone,
\verb{p(BE-ABE)} is the power of the conventional ABE test alone, and \verb{p(BE-sratio)}
is the power of the criterion \sQuote{upper confidence limit of
\ifelse{latex}{\out{$\sigma_{\rm{wT}}/\sigma_{\rm{wR}}\leq 2.5$}}{\ifelse{html}{\out{&sigma;<sub>wT</sub>/&sigma;<sub>wR</sub>&nbsp;&leq;&nbsp;2.5}}{sWT/sWR <= 2.5}}} alone.
}
\references{
Food and Drug Administration, Office of Generic Drugs (OGD). \emph{Draft Guidance on Warfarin Sodium.} Recommended Dec 2012. \href{https://www.accessdata.fda.gov/drugsatfda_docs/psg/Warfarin_Sodium_tab_09218_RC12-12.pdf}{download}

Food and Drug Administration, Center for Drug Evaluation and Research (CDER). \emph{Draft Guidance for Industry. Bioequivalence Studies with Pharmacokinetic Endpoints for Drugs Submitted Under an ANDA.} August 2021. \href{https://www.fda.gov/media/87219/download}{download}

Yu LX, Jiang W, Zhang X, Lionberger R, Makhlouf F, Schuirmann DJ, Muldowney L, Chen ML, Davit B, Conner D, Woodcock J. \emph{Novel bioequivalence approach for narrow therapeutic index drugs.} Clin Pharmacol Ther. 2015;97(3):286--91. \doi{10.1002/cpt.28}

Jiang W, Makhlouf F, Schuirmann DJ, Zhang X, Zheng N, Conner D, Yu LX, Lionberger R. \emph{A Bioequivalence Approach for Generic Narrow Therapeutic Index Drugs: Evaluation of the Reference-Scaled Approach and Variability Comparison Criterion.} AAPS J. 2015;17(4):891--901. \doi{10.1208/s12248-015-9753-5}

\enc{Endrényi}{Endrenyi} L, \enc{Tóthfalusi}{Tothfalusi} L. \emph{Determination of Bioequivalence for Drugs with Narrow Therapeutic Index: Reduction of the Regulatory Burden.} J Pharm Pharm Sci. 2013;16(5):676--82. \href{https://journals.library.ualberta.ca/jpps/index.php/JPPS/article/download/20900/15927/0}{open access}
}
\author{
D. Labes
}
\note{
The FDA\enc{’}{'}s method is described for the \acronym{ABE} limits 0.8 ... 1.25 only. Setting \code{theta1},
\code{theta2} to other values may not be reasonable and is not tested.\cr
The results for the design \code{"2x2x3"} are to be considered as experimental since at present not thorougly tested.\cr\cr
The method is also required by China\enc{’}{'}s Center of Drug Evaluation.
}
\seealso{
\code{\link{sampleN.NTID}}\cr
and \code{\link{power.HVNTID}}, \code{\link{sampleN.HVNTID}} for NTIDs with
high variability
}
\examples{
# using the all defaults:
# GMR=0.975, theta1=0.8, theta2=1.25, 100,000 simulations
# and a CV of 0.1 (= 10\%) with 12 subjects, balanced
power.NTID(CV = 0.1, n = 12)
# should give a power of 0.62553
}
