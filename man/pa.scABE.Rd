\encoding{utf-8}
\name{pa.scABE}
\alias{pa.scABE}
\title{
Power analysis for scaled average bioequivalence (scABE)
}
\description{
An analysis tool for exploration/visualization of the impact of expected values
(CV, theta0, reduced sample size due to drop-outs) on power of BE decision via
scABE (for highly variable drugs) if these values deviate from the ones assumed
in planning the sample size of the study.
}
\usage{
pa.scABE(CV, theta0 = 0.9, targetpower = 0.8, minpower = 0.7,
         design = c("2x3x3", "2x2x4", "2x2x3"),
         regulator = c("EMA", "HC", "ANVISA", "FDA"), ...)
}
\arguments{
  \item{CV}{
Coefficient of variation of the intra-subject variability as ratio.\cr
Here only the case CVwT=CVwR is implemented, i.e. CV has to be a scalar.
}
  \item{theta0}{
'True' or assumed T/R ratio. Often named GMR.\cr
Must be given as ratio. Defaults to 0.9 here since HVD have a greater scatter
in point estimates of T/R.
}
  \item{targetpower}{
Power to achieve at least in sample size estimation. Must be >0 and <1.\cr
Typical values are 0.8 or 0.9. Defaults to 0.8.\cr
Note that targetpower < 0.5 doesn't make much sense.\cr
If package is run under R version < 3.1.0 \code{targetpower} has to be >= 0.5.
}
  \item{minpower}{
Minimum acceptable power to have if deviating from assumptions for sample size plan.\cr
Has to lower than \code{targetpower}. Defaults to 0.7.\cr
\code{minpower} < 0.5 doesn't make much sense.\cr
If package is run under R version < 3.1.0 \code{minpower} has to be >= 0.5.
}
  \item{design}{
Character string describing the study design.\cr
Defaults to \code{2x3x3}, the partial replicate design (TRR|RTR|RRT).
}
  \item{regulator}{
Character string describing the scaled ABE method recommended by the regulatory
bodies EMA, HC, ANVISA or FDA.\cr
Defaults to \code{EMA}, method of scaled (widened) bioequivalence limits.
}
  \item{\dots}{
More arguments to pass to \code{power.scABEL()} or \code{power.RSABE()}.\cr
F.i. \code{alpha}, \code{theta1}, \code{theta2} or \code{nsims} if other values
then the defaults for these arguments are needed. \cr
See man pages of \code{power.scABEL()} or \code{power.RSABE()}.
}
}
\details{
Power calculations are done via \code{power.scABEL()} or \code{power.RSABE()} and
calculations of CV and theta0 which result in \code{minpower} are derived via \code{uniroot()}.\cr
While one of the parameters (CV, GMR, N) is varied, the respective two others are
kept constant. The tool shows the relative impact of single parameters on power.\cr
The tool takes a minimum of 12 subjects as demanded in most BE guidances into account.\cr\cr
It should be kept in mind that this is \bold{not} a substitute for the "Sensitivity Analysis"
recommended in ICH-E9. In a real study a combination of all effects occurs simultaneously.
It is up to \bold{you} to decide on reasonable combinations and analyze their respective power.
}
\value{
Returns a list with class \code{'pwrA'} with the components
  \item{plan}{A data.frame with the result of the sample size estimation.
              See output of \code{sampleN.scABEL()} or \code{sampleN.RSABE()}.}
  \item{paCV}{A data.frame with value pairs CV, pwr for impact of
              deviations from CV.}
  \item{paGMR}{A data.frame with value pairs theta0, pwr for impact of
               deviations from theta0 (GMR).}
  \item{paN}{A data.frame with value pairs N, pwr for impact of
             deviations from planned N (drop-outs).}
  \item{method}{Method of BE decision. Here fix = "scABE".}
  \item{regulator}{"EMA" or "FDA".}
  \item{minpower}{Minimum acceptable power from the call of the function.}

The class \code{'pwrA'} has the S3 methods \code{print()} and \code{plot()}.
See \code{\link{pa.ABE}} for usage.
}
\references{
See \url{http://forum.bebac.at/mix_entry.php?id=13376}.
}
\author{
Idea and original code by H. \enc{Schütz}{Schuetz}\cr
with modifications by D. Labes to use PowerTOST infrastructure.
}
\note{
The code for impact of deviations from planned sample size tries
to keep the degree of imbalance as low as possible between (sequence) groups.
This results in a lesser drop of power than more extreme drop-out patterns.
}
\seealso{
  \code{\link{power.scABEL}, \link{power.RSABE}, \link{known.designs},
  \link{print.pwrA}, \link{plot.pwrA}, \link{pa.ABE}, \link{pa.NTIDFDA}}
}
\examples{
# using the defaults:
# design="2x3x3", targetpower=0.8, minpower=0.7, theta0/GMR=0.90
# widened BE margins from defaults of sampleN.scABEL() 0.7462 ... 1.3402
# 1E5 sims in power.scABEL()
# not run due to timing policy of CRAN, may run some ten seconds
\dontrun{
# implicit print & plot
pa.scABE(CV=0.4)}
}
