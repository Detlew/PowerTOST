\encoding{utf-8}
\name{pa.NTIDFDA}
\alias{pa.NTIDFDA}
\title{
Power analysis for scaled ABE for NTIDs according to FDA
}
\description{
An analysis tool for exploration/visualization of the impact of expected values
(CV, theta0, reduced sample size due to drop-outs) on power of BE decision via
scABE for narrow therapeutic drugs (NTIDs) if these values deviate from the ones
assumed in planning the sample size of the study.\cr
The only implemented design is the full replicate design "2x2x4" according to the
FDA Warfarin guidance.
}
\usage{
pa.NTIDFDA(CV, theta0 = 0.975, targetpower = 0.8, minpower = 0.7, ...)
}
\arguments{
  \item{CV}{
Coefficient of variation of the intra-subject variabilities of Test and Reference
as ratio.\cr
Here only the case CVwT=CVwR is implemented, i.e. CV has to be a scalar.
}
  \item{theta0}{
\sQuote{True} or assumed T/R ratio. Often named GMR.\cr
Must be given as ratio. Defaults here to 0.975.
}
  \item{targetpower}{
Power to achieve at least in sample size estimation. Must be >0 and <1.\cr
Typical values are 0.8 or 0.9. Defaults to 0.8.\cr
Note that targetpower < 0.5 doesn’t make much sense.
}
  \item{minpower}{
Minimum acceptable power to have if deviating from assumptions for sample size plan.\cr
Has to be lower than \code{targetpower}. Defaults to 0.7.\cr
\code{minpower} < 0.5 doesn’t make much sense.
}
  \item{\dots}{
More arguments to pass to \code{power.NTIDFDA()}.\cr
F.i. \code{alpha}, \code{theta1}, \code{theta2} or \code{nsims} if other values
then the defaults for these arguments are needed. \cr
See man page of \code{power.NTIDFDA()}.
}
}
\details{
Power calculations are done via \code{power.NTIDFDA()} and
calculations of CV and theta0 which result in \code{minpower} are derived via \code{uniroot()}.\cr
While one of the parameters (CV, theta0, n) is varied, the respective two others are
kept constant. The tool shows the relative impact of single parameters on power.\cr
The tool takes a minimum of 12 subjects into account as demanded in most BE guidances.
However, it should be kept in mind that the FDA requires at least 24 subjects to be enrolled
in studies intended for reference-scaling.\cr\cr
It should be kept in mind that this is \bold{not} a substitute for the \dQuote{Sensitivity Analysis}
recommended in ICH-E9. In a real study a combination of all effects occurs simultaneously.
It is up to \emph{you} to decide on reasonable combinations and analyze their respective power.
}
\value{
Returns a list with class \code{'pwrA'} with the components
  \item{plan}{A data.frame with the result of the sample size estimation.
              See output of \code{sampleN.NTIDFDA()}.}
  \item{paCV}{A data.frame with value pairs CV, pwr for impact of
              deviations from CV.}
  \item{paGMR}{A data.frame with value pairs theta0, pwr for impact of
               deviations from theta0 (GMR).}
  \item{paN}{A data.frame with value pairs N, pwr for impact of
             deviations from planned N (dropouts).}
  \item{method}{Method of BE decision. Here fix = "NTID FDA".}
  \item{regulator}{Here fix = "FDA".}
  \item{minpower}{Minimum acceptable power from the call of the function.}

The class \code{'pwrA'} has the S3 methods \code{print()} and \code{plot()}.
See \code{\link{pa.ABE}} for usage.
}
\references{
FDA. \emph{Draft Guidance on Warfarin Sodium.}\cr
Recommended Dec 2012\cr
\url{http://www.fda.gov/downloads/Drugs/GuidanceComplianceRegulatoryInformation/Guidances/UCM201283.pdf}

Yu LX et al.\cr
\emph{Novel bioequivalence approach for narrow therapeutic index drugs}\cr
Clin Pharmacol Ther. 2015;97(3):286--91. \doi{10.1002/cpt.28}

Jiang W et al.\cr
\emph{A Bioequivalence Approach for Generic Narrow Therapeutic Index Drugs: Evaluation of the Reference-Scaled Approach and Variability Comparison Criterion}\cr
AAPS J. 2015;17(4):891--901. \doi{10.1208/s12248-015-9753-5}

\enc{Endrényi}{Endrenyi} L, \enc{Tóthfalusi}{Tothfalusi} L.\cr
\emph{Determination of Bioequivalence for Drugs with Narrow Therapeutic Index: 
Reduction of the Regulatory Burden}\cr
J Pharm Pharm Sci. 2013;16(5):676--82. \href{https://ejournals.library.ualberta.ca/index.php/JPPS/article/viewFile/20900/15927}{free download}
}
\author{
D. Labes\cr
according to code by H. \enc{Schütz}{Schuetz} for pa.ABE() and pa.scABE()
}
\note{
The code for impact of deviations from planned sample size tries
to keep the degree of imbalance as low as possible between (sequence) groups.
This results in a lesser decrease of power than more extreme dropout-patterns.
}
\section{Warning }{
Be extremly carefull if your sample size plan has extremly small CV near or
below 0.05 (5\%). Adapt in that case your expected true ratio (\code{theta0})
to values nearer to 1 to not run into errors and/or long execution times.
}
\seealso{
  \code{\link{power.NTIDFDA}, \link{print.pwrA}, \link{plot.pwrA}, \link{pa.ABE}, \link{pa.scABE}}
}
\examples{
# using the defaults:
# targetpower=0.8, minpower=0.7, theta0/GMR=0.975
# BE margins from defaults of sampleN.NTIDFDA() 0.9002 ... 1.1108
# 1E5 sims in power.NTIDFDA()
# not run due to timing policy of CRAN for examples
# may run some ten seconds or more
\dontrun{
plot(pa.NTIDFDA(CV=0.1))}
}