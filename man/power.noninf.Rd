\encoding{utf-8}
\name{power.noninf}
\alias{power.noninf}
\title{
Power of the one-sided non-inferiority t-test
}
\description{
Function calculates of the power of the one-sided non-inferiority \emph{t}-test
for normal or log-normal distributed data.
}
\usage{
power.noninf(alpha = 0.025, logscale = TRUE, margin, theta0, CV, n, 
             design = "2x2", robust = FALSE)
}
\arguments{
  \item{alpha}{
Type I error probability, significance level. Defaults here to 0.025.
}
  \item{logscale}{
Should the data used on log-transformed or on original scale? \code{TRUE} (default) or \code{FALSE}.
}
  \item{margin}{
Non-inferiority margin.\cr
In case of \code{logscale=TRUE} it must be given as ratio, otherwise as difference to 1.\cr
Defaults to 0.8 if \code{logscale=TRUE} or to -0.2 if \code{logscale=FALSE}.
}
  \item{theta0}{
\sQuote{True} or assumed T/R ratio or difference (T--R). \cr
In case of \code{logscale=TRUE} it must be given as ratio, \cr
otherwise as difference to 1. See examples.\cr
Defaults to 0.95 if \code{logscale=TRUE} or to -0.05 if \code{logscale=FALSE}.
}
  \item{CV}{
Coefficient of variation as ratio.\cr
In case of cross-over studies this is the within-subject CV and in case of a parallel-group design the CV of the total variability.
}
  \item{n}{
Number of subjects under study.\cr
Is total number if given as scalar, else number of subjects in the (sequence) 
groups. In the latter case the length of the \code{n} vector has to be equal to the 
number of (sequence) groups.
}
  \item{design}{
Character string describing the study design.\cr
See \code{\link{known.designs}} for designs covered in this package.
}
  \item{robust}{
Defaults to \code{FALSE}. With that value the usual degrees of freedom will be used.\cr
Set to \code{TRUE} will use the degrees of freedom according to the \sQuote{robust} evaluation
(aka \enc{Senn’s}{Senn's} basic estimator). These degrees of freedom are calculated as \code{n-seq}.\cr
See \code{known.designs()$df2} for designs covered in this package.\cr
Has only effect for higher-order crossover designs.
}
}
\details{
The power is calculated exact via non-central \emph{t}-distribution.\cr

\strong{Notes on the underlying hypotheses}\cr
If the supplied margin is < 0 (\code{logscale=FALSE}) or < 1 (\code{logscale=TRUE}), 
then it is assumed higher response values are better. The hypotheses are\cr
\verb{  H0: theta0 <= margin vs. H1: theta0 > margin}\cr
where \code{theta0 = mean(test)-mean(reference)} if \code{logscale=FALSE}\cr
or\cr
\verb{  H0: log(theta0) <= log(margin) vs. H1: log(theta0) > log(margin)}\cr
where \code{theta0 = mean(test)/mean(reference)} if \code{logscale=TRUE}.\cr

If the supplied margin is > 0 (\code{logscale=FALSE}) or > 1 (\code{logscale=TRUE}), 
then it is assumed lower response values are better. The hypotheses are\cr
\verb{  H0: theta0 >= margin vs. H1: theta0 < margin}\cr
where \code{theta0 = mean(test)-mean(reference)} if \code{logscale=FALSE}\cr
or\cr
\verb{  H0: log(theta0) >= log(margin) vs. H1: log(theta0) < log(margin)}\cr
where \code{theta0 = mean(test)/mean(reference)} if \code{logscale=TRUE}.\cr
This latter case may also be considered as \sQuote{non-superiority}.
}
\value{
Value of power according to the input arguments.
}
\section{Warning }{
The function does not vectorize if design is a vector.\cr
The function vectorize properly if CV or theta0 are vectors.\cr
Other vector input is not tested yet.
}
\references{
Julious SA. \emph{Sample sizes for clinical trials with Normal data.} Stat Med. 2004;23(12):1921--86. \doi{10.1002/sim.1783}
}
\author{
D. Labes
}
\note{
This function does not rely on TOST but may be useful in planning BE studies 
if the question is not equivalence but \sQuote{non-superiority}.\cr
Hint: Evaluation of Fluctuation in the \href{https://www.ema.europa.eu/en/documents/scientific-guideline/note-guidance-modified-release-oral-transdermal-dosage-forms-section-ii-pharmacokinetic-clinical_en.pdf}{\enc{EMEA’s}{EMEA's} Note for Guidance} between a modified release 
formulation and an immediate release product.
}

\seealso{
\code{\link{known.designs}, \link{sampleN.noninf}}
}
\examples{
# using all the defaults: margin=0.8, theta0=0.95, alpha=0.025
# log-transformed, design="2x2"
# should give: 0.4916748
power.noninf(CV=0.3, n=24)
}
