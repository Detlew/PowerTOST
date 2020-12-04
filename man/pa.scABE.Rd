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
         regulator = c("EMA", "HC", "FDA"), ...)
}
\arguments{
  \item{CV}{
Coefficient of variation of the intra-subject variability as ratio (not percent).\cr
Here only the case CVwT=CVwR is implemented, \emph{i.e.}, CV has to be a scalar.
}
  \item{theta0}{
\sQuote{True} or assumed T/R ratio. Often named GMR.\cr
Must be given as ratio. Defaults to 0.9 here since HVD have a greater scatter
in point estimates of T/R.
}
  \item{targetpower}{
Power to achieve at least in sample size estimation. Must be >0 and <1.\cr
Typical values are 0.8 or 0.9. Defaults to 0.8.\cr
Note that targetpower < 0.5 doesn’t make much sense.
}
  \item{minpower}{
Minimum acceptable power to have if deviating from assumptions for sample size plan.\cr
Has to lower than \code{targetpower}. Defaults to 0.7.\cr
\code{minpower} < 0.5 doesn’t make much sense.
}
  \item{design}{
Character string describing the study design.\cr
Defaults to \code{2x3x3}, the partial replicate design (TRR|RTR|RRT).
}
  \item{regulator}{
Character string describing the scaled ABE method recommended by the regulatory
bodies \code{"EMA"}, \code{"HC"}, \code{"FDA"} or \code{"GCC"}.\cr
Defaults to \code{"EMA"}, method of scaled (expanded) bioequivalence limits.
}
  \item{\dots}{
More arguments to pass to \code{power.scABEL()} or \code{power.RSABE()}.\cr
F.i., \code{alpha}, \code{theta1}, \code{theta2} or \code{nsims} if other values
than the defaults for these arguments are needed.\cr
See man pages of \code{power.scABEL()} or \code{power.RSABE()}.
}
}
\details{
Power calculations are done via \code{power.scABEL()} or \code{power.RSABE()} and
calculations of CV and theta0 which result in \code{minpower} derived via R base \code{\link{uniroot}}.\cr
While one of the parameters (CV, GMR, N) is varied, the respective two others are
kept constant. The tool shows the relative impact of single parameters on power.\cr\cr
The tool takes a minimum of 12 subjects as required in most BE guidances into account. 
However, the sample size will be increased from the estimated one if one of the 
following conditions is applicable:
  \itemize{
    \item The FDA requires at least 24 subjects \emph{enrolled} in studies intended for reference-scaling.
    \item The EMA requires at least 12 \emph{eligible} subjects in the sequence RTR of the TRT|RTR-design (hence the minimum sample size is 24).
  }
You should be aware that this is \bold{not} a substitute for the \dQuote{Sensitivity Analysis}
recommended in ICH-E9. In a real study a combination of all effects occurs simultaneously.
It is up to \emph{you} to decide on reasonable combinations and analyze their respective power.
}
\value{
Returns a list with class \code{'pwrA'} with the components
  \item{plan}{A data.frame with the result of the sample size estimation.\cr
              See output of \code{sampleN.scABEL()} or \code{sampleN.RSABE()}.}
  \item{paCV}{A data.frame with value pairs CV, pwr for impact of
              deviations from CV.}
  \item{paGMR}{A data.frame with value pairs theta0, pwr for impact of
               deviations from theta0 (GMR).}
  \item{paN}{A data.frame with value pairs N, pwr for impact of
             deviations from planned N (dropouts).}
  \item{method}{Method of BE decision. Here fix = "scABE".}
  \item{regulator}{"EMA", "HC", or "FDA".}
  \item{minpower}{Minimum acceptable power from the call of the function.}

The class \code{'pwrA'} has the S3 methods \code{print()} and \code{plot()}.
See \code{\link{pa.ABE}} for usage.
}
\references{
\enc{Schütz}{Schuetz} H. \emph{Deviating from assumptions.} August 08, 2014. 
\href{https://forum.bebac.at/mix_entry.php?id=13353}{BEBA Forum}
}
\author{
Idea and original code by H. \enc{Schütz}{Schuetz} with modifications by D. Labes 
to use PowerTOST infrastructure.
}
\note{
The code for impact of deviations from planned sample size tries to
keep the degree of imbalance as low as possible between (sequence) groups.
This results in a lesser decrease of power than more extreme dropout-patterns.
}
\seealso{
  \code{\link{power.scABEL}, \link{power.RSABE}, \link{known.designs},
  \link{print.pwrA}, \link{plot.pwrA}, \link{pa.ABE}, \link{pa.NTIDFDA}}
}
\examples{
# Implicitely using the defaults:
# design = "2x3x3", targetpower = 0.8, minpower = 0.7,
# theta0 = 0.9, GMR = 0.90, regulator = "EMA"
# widened BE margins from defaults of sampleN.scABEL() 0.7462 ... 1.3402
# 1E5 sims in power.scABEL()
# not run due to timing policy of CRAN, may run some ten seconds
\donttest{
# Implicit print & plot
pa.scABE(CV = 0.4)}
}