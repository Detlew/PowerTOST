\encoding{utf-8}
\name{pa.ABE}
\alias{pa.ABE}
\alias{print.pwrA}
\alias{plot.pwrA}
\title{
Power analysis for average bioequivalence (ABE)
}
\description{
An analysis tool for exploration/visualization of the impact of expected values
(CV, GMR, reduced sample size due to drop-outs) on power of BE decision via ABE
if these values deviate from the ones assumed in planning the sample size of
the study.
}
\usage{
pa.ABE(CV, theta0 = 0.95, targetpower = 0.8, minpower = 0.7, design = "2x2", \dots)
\method{print}{pwrA}(x, digits = 4, plotit = TRUE, \dots)
\method{plot}{pwrA}(x, pct = TRUE, ratiolabel = "theta0", cols = c("blue", "red"), \dots)
}
\arguments{
  \item{CV}{
Coefficient of variation as ratio.\cr
In case of cross-over studies this is the within-subject CV;\cr
in case of a parallel-group design the CV of the total (pooled) variability.
}
  \item{theta0}{
'True' or assumed T/R ratio. Often named GMR.\cr
Must be given as ratio.
}
  \item{targetpower}{
Power to achieve at least in sample size estimation. Must be >0 and <1.\cr
Typical values are 0.8 or 0.9. Defaults to 0.8.\cr
Note that targetpower < 0.5 doesn't make much sense.\cr
If package is run under R version < 3.1.0 \code{targetpower} has to be >= 0.5.
}
  \item{minpower}{
Minimum acceptable power to have if deviating from assumptions for sample size plan.\cr
Has to be lower than \code{targetpower}. Defaults to 0.7.\cr
\code{minpower} < 0.5 doesn't make much sense.\cr
If package is run under R version < 3.1.0 \code{minpower} has to be >= 0.5.
}
  \item{design}{
Character string describing the study design.\cr
See \code{known.designs()} for designs covered in this package.
}
  \item{\dots}{
More arguments to pass to \code{power.TOST()}.\cr
F.i. \code{alpha}, \code{theta1}, \code{theta2} or \code{robust} if other values
then the defaults for these arguments are needed. \cr
See man page of \code{power.TOST()}.\cr\cr
More arguments passed to the S3 methods. Here currently ignored.
}
Additional arguments of the S3 methods:
  \item{x}{
Object of class \code{'pwrA'}.
}
  \item{digits}{
Digits for rounding power in printing. The '\dots' argument is currently ignored
in \code{print()}.
}
  \item{plotit}{
If set to \code{TRUE}, the default, the print method calls \code{plot(x)} if R
is running interactively.
}
  \item{pct}{
If set to \code{TRUE} (the default) scales CV, theta0, and power in percent in
\code{plot()}. Else they will be given as ratios, the usual standard in PowerTOST.
}
  \item{ratiolabel}{
Label of the T/R-ratio. Can be set to any string, e.g. to \code{"GMR"}.
Defaults to \code{"theta0"}, the usual standard in PowerTOST.
}
  \item{cols}{
Colors for the plots. \code{cols[1]} gives the color for plotting points
with \code{power>targetpower}. From \code{targetpower} toward \code{minpower}
the color changes gradually to \code{cols[2]}.}
}
\details{
Power calculations are done via \code{power.TOST()} and calculations of CV and theta0
which gave a power=\code{minpower} are derived via R base \code{uniroot()}.
While one of the parameters (CV, theta0, N) is varied, the respective two others are
kept constant. The tool shows the relative impact of single parameters on power.\cr
The tool takes a minimum of 12 subjects as demanded in most BE guidances into account.\cr\cr
It should be kept in mind that this is \bold{not} a substitute for the "Sensitivity Analysis"
recommended in ICH-E9. In a real study a combination of all effects occurs simultaneously.
It is up to \bold{you} to decide on reasonable combinations and analyze their respective power.
}
\value{
Returns a list with class \code{"pwrA"} with the components
  \item{plan}{A data.frame with the result of the sample size estimation.
              See output of \code{sampleN.TOST()}.}
  \item{paCV}{A data.frame with value pairs CV, pwr for impact of
              deviations from CV.}
  \item{paGMR}{A data.frame with value pairs theta0, pwr for impact of
               deviations from theta0 (GMR).}
  \item{paN}{A data.frame with value pairs N, pwr for impact of
             deviations from planned N (drop-outs).}
  \item{method}{Method of BE decision. Here fix = "ABE".}
  \item{minpower}{Minimum acceptable power.}

The class \code{'pwrA'} has the S3 methods \code{print()} and \code{plot()}.
See \code{\link{pa.scABE}} for usage.
}
\references{
See \url{http://forum.bebac.at/mix_entry.php?id=13353}.
}
\author{
Idea and original code by H. \enc{Schütz}{Schuetz}\cr
with modifications by D. Labes to use PowerTOST infrastructure.
}
\note{
The code of deviations from planned sample size tries to keep the degree of
imbalance as low as possible between (sequence) groups.
This results in a lesser drop of power than more extreme drop-out patterns.
}

\seealso{
\code{\link{power.TOST}}, \code{\link{known.designs}}, \code{\link{pa.scABE}}
}
\examples{
# using the defaults
# design="2x2", targetpower=0.8, minpower=0.7, theta0/GMR=0.95
# BE margins from defaults of sampleN.TOST() 0.8 ... 1.25
# print & plot implicitly
pa.ABE(CV=0.2)
# print & plot
\dontrun{
res <- pa.ABE(CV=0.2)
print(res, plotit=FALSE)               # print only
plot(res, pct=FALSE, ratiolabel="GMR") # changed from defaults}
}
