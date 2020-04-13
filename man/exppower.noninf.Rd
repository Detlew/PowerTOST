\encoding{utf-8}
\name{exppower.noninf}
\alias{exppower.noninf}
\title{
Expected power of the non-inferiority test
}
\description{
Calculates the so-called expected, \emph{i.e.}, unconditional, power for a variety of
study designs used in bioequivalence studies.
}
\usage{
exppower.noninf(alpha = 0.025, logscale = TRUE, theta0, margin, CV, n,
                design = "2x2", robust = FALSE, 
                prior.type = c("CV", "theta0", "both"), prior.parm = list(), 
                method = c("exact", "approx"))
}
\arguments{
  \item{alpha}{
Significance level (one-sided). Defaults here to 0.025.
  }
  \item{logscale}{
Should the data be used on log-transformed or on original scale? \code{TRUE} (default) or \code{FALSE}.
  }  
  \item{theta0}{
Assumed \sQuote{true} (or \sQuote{observed} in case of \code{prior.type != "CV"}) ratio or difference.\cr
In case of \code{logscale=TRUE} it must be given as ratio T/R.\cr
If \code{logscale=FALSE}, the difference in means. In this case, the difference may be expressed in two ways: relative to the same (underlying) reference mean, i.e. as (T-R)/R = T/R - 1; or as difference in means T-R. Note that in the former case the units of \code{margin} and \code{CV} need also be given relative to the reference mean (specified as ratio).\cr
Defaults to 0.95 if \code{logscale=TRUE} or to -0.05 if \code{logscale=FALSE}
  }
  \item{margin}{
Non-inferiority margin.\cr
In case of \code{logscale=TRUE} it is given as ratio.\cr
If \code{logscale=FALSE}, the limit may be expressed in two ways:
difference of means relative to the same (underlying) reference mean or in units of the difference of means.
Note that in the former case the units of \code{CV} and \code{theta0} need also be given relative to the reference mean (specified as ratio).\cr
Defaults to 0.8 if \code{logscale=TRUE} or to -0.2 if \code{logscale=FALSE}.
  }
  \item{CV}{
In case of \code{logscale=TRUE} the (geometric) coefficient of variation given as ratio.\cr
If \code{logscale=FALSE} the argument refers to (residual) standard deviation of the response. In this case, standard deviation may be expressed two ways: relative to a reference mean (specified as ratio sigma/muR), i.e. again as a coefficient of variation; or untransformed, i.e. as standard deviation of the response. Note that in the former case the units of \code{theta0}, \code{theta1} and \code{theta2} need also be given relative to the reference mean (specified as ratio).\cr

In case of cross-over studies this is the within-subject CV, in case of a parallel-group design the CV of the total variability.
  }
  \item{n}{
Number of subjects under study.\cr
Is total number if given as scalar, else number of subjects in the (sequence) 
groups. In the latter case the length of n has to be equal to the 
number of (sequence) groups.
}
  \item{design}{
Character string describing the study design.
See \code{known.designs()} for designs covered in this package.
}
  \item{robust}{
Defaults to FALSE. Set to \code{TRUE} will use the degrees of freedom according 
to the \sQuote{robust} evaluation (aka \enc{Senn’s}{Senn's} basic estimator). These degrees of freedom are calculated
as \code{n-seq}.\cr 
See \code{known.designs()$df2} for designs covered in this package.
}
  \item{prior.type}{
Specifies which parameter uncertainty should be accounted for. In case of 
\code{prior.type = "CV"} (the default), only the uncertainty with respect to the
CV will be considered (i.e. the given treatment effect is assumed to be fix).
In case of \code{prior.type = "theta0" } only uncertainty with respect to the
treatment ratio/difference will be accounted for (i.e. the given CV is assumed
to be fix). In case of \code{prior.type = "both"} the power value will be
unconditional with respect to both the \code{CV} and \code{theta0}.
  }
  \item{prior.parm}{
A list of parameters expressing the prior information about the 
variability and/or treatment effect. Possible components are \code{df}, 
\code{SEM}, \code{m} and \code{design}.\cr
For \code{prior.type = "CV"} the degrees of freedom from the prior trial are
required. This information can be provided by specifying the single component 
\code{df} or the combination consisting of \code{m} and \code{design}.\cr
For \code{prior.type = "theta0"} the standard error of the treatment difference 
from the prior trial is required. This information can be provided by specifying
the single component \code{SEM} or the combination consisting of 
\code{m} and \code{design}.\cr
For \code{prior.type = "both"} the degrees of freedom and the standard error of
the treatment difference are required. This information can be provided by 
specifying the combination consisting of \code{df} and \code{SEM} or via the 
combination \code{m} and \code{design}.\cr
See 'Details' for a technical description on each component.
  }
  \item{method}{
Defaults to \code{method="exact"}.
In that case the expected power will be calculated as expected value of the
power with respect to the (prior) distribution of the respective parameter(s).\cr
Set to \code{method="approx"} the expected power according to the 
approximate formulas given in the book from Julious or in the Julious/Owen 
paper will be calculated (using non-central \emph{t}); this only affects 
\code{prior.type = "CV"}.
}
}
\details{
This function calculates the so-called expected power taking into account that
usually the parameters (CV and/or theta0) are not known but \emph{estimated} 
from a prior study with some uncertainty. The expected power is an unconditional 
power and can therefore be seen as probability for success. See references for 
further details.\cr\cr
The \code{prior.parm} argument is a list that can supply any of the following 
components: 
\describe{
\item{\code{df}}{Error degrees of freedom from the prior trial (>4, maybe non-integer). 
\code{df = Inf} is allowed and for \code{method = "exact"} the result will then
coincide with \code{power.noninf(...)}.\cr
Note: This corresponds to the df of both the CV and the difference of means.}
\item{\code{SEM}}{Standard error of the difference of means from the prior trial; 
must always be on additive scale (\emph{i.e.}, usually log-scale).}
\item{\code{m}}{Number of subjects from prior trial. Specification is analogous to
the main argument \code{n}.}
\item{\code{design}}{Study design of prior trial. Specification is analogous to the 
main argument \code{design}.}
}

For \code{prior.parm}, the combination consisting of \code{df} and \code{SEM} 
requires a somewhat advanced knowledge of the prior trial (provided in the raw
output from for example the software \acronym{SAS}, or may be obtained via
\code{emmeans::emmeans}). However, it has the advantage that if there were 
missing data the exact degrees of freedom and standard error of the difference 
can be used, the former possibly being non-integer valued (e.g. if the 
Kenward-Roger method was used).

Details on argument \code{prior.type}:
\describe{
  \item{\code{CV}}{The expectation is calculated with respect to the 
  Inverse-gamma distribution.}
  \item{\code{theta0}}{The expectation is calculated with respect to the
  conditional distribution theta0 | \eqn{\sigma^2} = s^2 
  of the posteriori distribution of (theta0, \eqn{\sigma^2}) from the prior
  trial.}
  \item{\code{both}}{The expectation is calculated with respect to the posteriori
  distribution of (theta0, \eqn{\sigma^2}) from the prior trial. Numerical calculation
  of the two-dimensional integral is performed via \code{cubature::hcubature}.}
}
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
or \cr
\verb{  H0: log(theta0) >= log(margin) vs. H1: log(theta0) < log(margin)}\cr
where \code{theta0 = mean(test)/mean(reference)} if \code{logscale=TRUE}.\cr
This latter case may also be considered as \sQuote{non-superiority}.
}

\value{
Value of expected power according to the input.
}
\references{
Grieve AP. \emph{Confidence Intervals and Sample Sizes.} Biometrics. 1991;47:1597--603. \doi{10.2307/2532411}

\enc{O’Hagan}{O'Hagan}, Stevens, JW, Campell MJ. \emph{Assurance in Clinical Trial Design.} Pharm Stat. 2005;4:187--201. \doi{10.1002/pst.175}

Julious SA, Owen RJ. \emph{Sample size calculations for clinical studies allowing for uncertainty in variance.} Pharm Stat. 2006;5:29--37. \doi{10.1002/pst.197}

Julious SA. \emph{Sample sizes for Clinical Trials.} Boca Raton: CRC Press / Chapman & Hall; 2010.

Bertsche A, Nehmitz G, Beyersmann J, Grieve AP. \emph{The predictive distribution of the residual variability in
the linear-fixed effects model for clinical cross-over trials.} Biom J. 2016;58(4):797--809. \doi{10.1002/bimj.201500245}

Box GEP, Tiao GC. \emph{Bayesian Inference in Statistical Analysis.} Boston: Addison-Wesley; 1992.

Held L, Sabanes Bove D. \emph{Applied Statistical Inference. Likelihood and Bayes.} Berlin, Heidelberg: Springer; 2014. \doi{10.1007/978-3-642-37887-4}

Senn S. \emph{Cross-over Trials in Clinical Research.} Chichester: John Wiley & Sons; 2\ifelse{html}{\out{<sup>nd</sup>}}{\eqn{^{nd}}} edition 2002.

Zierhut ML, Bycott P, Gibbs MA, Smith BP, Vicini P. \emph{Ignorance is not bliss: Statistical power is not probability of trial success.} Clin Pharmacol Ther. 2015;99:356--9. \doi{10.1002/cpt.257}
}
\author{
B. Lang, D. Labes
}
\seealso{
 \code{\link{expsampleN.noninf}, \link{power.noninf}}
}
\examples{
# Expected power for non-inferiority test for a 2x2 crossover
# with 40 subjects. CV 30\% known from a pilot 2x2 study with 
# 12 subjects 
# using all the defaults for other parameters (theta0 carved in stone)
# should give: [1] 0.6761068
exppower.noninf(CV = 0.3, n = 40, prior.parm = list(df = 12-2))
# or equivalently
exppower.noninf(CV = 0.3, n = 40, prior.parm = list(m = 12, design = "2x2"))

# May be also calculated via exppower.TOST() after setting upper acceptance limit 
# to Inf and alpha=0.025
exppower.TOST(CV = 0.3, n = 40, prior.parm = list(df = 10), theta2 = Inf, alpha=0.025)

# In contrast: Julious approximation
exppower.noninf(CV = 0.3, n = 40, prior.parm = list(df = 10), method = "approx")
# should give: [1] 0.6751358

# Compare this to the usual (conditional) power (CV known, "carved in stone")
power.noninf(CV = 0.3, n = 40)
# should give: [1] 0.7228685
# same as if setting df = Inf in function exppower.noninf()
exppower.noninf(CV = 0.3, n = 40, prior.parm = list(df = Inf))

# Expected power for a 2x2 crossover with 40 subjects
# CV 30\% and theta0 = 0.95 known from a pilot 2x2 study with 12 subjects
# using uncertainty with respect to both CV and theta0
exppower.noninf(CV = 0.3, theta0 = 0.95, n = 40, 
                prior.parm = list(m = 12, design = "2x2"), prior.type = "both")
# should give a decrease of expected power to 0.5982852
}
