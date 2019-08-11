\encoding{utf-8}
\name{exppower.TOST}
\alias{exppower.TOST}
\title{
Expected power of the TOST procedure
}
\description{
Calculates the so-called expected, \emph{i.e.}, unconditional, power for a variety of
study designs used in bioequivalence studies.
}
\usage{
exppower.TOST(alpha = 0.05, logscale = TRUE, theta0, theta1, theta2,  
              CV, n, design = "2x2", robust = FALSE, 
              prior.type = c("CV", "theta0", "both"), prior.parm = list(),
              method = c("exact", "approx"))
}
\arguments{
  \item{alpha}{
Significance level (one-sided). Commonly set to 0.05.
  }
  \item{logscale}{
Should the data be used on log-transformed or on original scale? \code{TRUE} (default) or \code{FALSE}.
  }  
  \item{theta0}{
Assumed \sQuote{true} (or \sQuote{observed} in case of \code{prior.type != "CV"}) bioequivalence 
ratio or difference.\cr
Typically set to 0.95 (default if missing) if \code{logscale=TRUE}.
Defaults to 0.05 if \code{logscale=FALSE}.
  }
  \item{theta1}{
Lower bioequivalence limit as ratio (if \code{logscale=TRUE}) or as difference.\cr
Can be missing. Defaults then to 0.8 if \code{logscale=TRUE} or to -0.2 if 
\code{logscale=FALSE}.
  }
  \item{theta2}{
Upper bioequivalence limit as ratio (if \code{logscale=TRUE}) or as difference.\cr
If not given \code{theta2} will be calculated as \code{1/theta1} if \code{logscale=TRUE},
else as \code{-theta1}.
  }
  \item{CV}{
Assumed true or observed coefficient of variation as ratio. Only values > 0 are
allowed. If \code{logscale=FALSE} CV is assumed to be the standard deviation.
  }
  \item{n}{
Number of subjects under study.\cr
Is total number if given as scalar, else number of subjects in the (sequence) 
groups. In the latter case the length of n has to be equal to the 
number of (sequence) groups.
}
  \item{design}{
Character string describing the study design.
See \link{known.designs} for designs covered in this package.
}
  \item{robust}{
Defaults to \code{FALSE}. Set to \code{TRUE} will use the degrees of freedom according 
to the \sQuote{robust} evaluation (aka Senn’s basic estimator). These df are calculated
as \code{n-seq}.\cr 
See \code{known.designs()$df2} for designs covered in this package.
}
  \item{prior.type}{
Specifies which parameter uncertainty should be accounted for. In case of 
\code{prior.type = "CV"} (the default), only the uncertainty with respect to the
CV will be considered (i.e. the given treatment effect is assumed to be fix).
In case of \code{prior.type = "theta0" } only uncertainty with respect to the
treatment ratio/difference will be accounted for (i.e. the given CV is assumed
to be fix). In case of  \code{prior.type = "both"} the power value will be
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
usually the parameters (CV and/or theta0) are not known but estimated from a
prior study with some uncertainty. The expected power is an unconditional power
and can therefore be seen as probability for success. See references for further
details.

The \code{prior.parm} argument is a list that can supply any of the following 
components: 
\describe{
\item{\code{df}}{Error degrees of freedom from the prior trial (>4, maybe non-integer). 
\code{df = Inf} is allowed and for \code{method = "exact"} the result will then
coincide with \code{power.TOST(...)}.\cr
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
\code{\link[emmeans]{emmeans}}. However, it has the advantage that if there were 
missing data the exact degrees of freedom and standard error of the difference 
can be used, the former possibly being non-integer valued (\emph{e.g.}, if the 
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
of the two-dimensional integral is performed via \code{\link[cubature]{hcubature}}.}
}
}
\value{
Value of expected power according to the input.
}
\references{
Grieve AP. \emph{Confidence Intervals and Sample Sizes.} Biometrics. 1991;47:1597--603. \doi{10.2307/2532411}

O’Hagan, Stevens, JW, Campell MJ. \emph{Assurance in Clinical Trial Design.} Pharm Stat. 2005;4:187--201. \doi{10.1002/pst.175}

Julious SA, Owen RJ. \emph{Sample size calculations for clinical studies allowing for 
uncertainty in variance.} Pharm Stat. 2006;5:29--37. \doi{10.1002/pst.197}

Julious SA. \emph{Sample sizes for Clinical Trials.} Boca Raton: CRC Press / Chapman & Hall; 2010.

Bertsche A, Nehmitz G, Beyersmann J, Grieve AP. \emph{The predictive distribution of the residual variability in
the linear-fixed effects model for clinical cross-over trials.} Biom J. 2016;58(4):797--809. \doi{10.1002/bimj.201500245}

Box GEP, Tiao GC. \emph{Bayesian Inference in Statistical Analysis.} Boston: Addison-Wesley; 1992.

Held L, Sabanes Bove D. \emph{Applied Statistical Inference. Likelihood and Bayes.} Berlin, Heidelberg: Springer; 2014. \doi{10.1007/978-3-642-37887-4}

Senn S. \emph{Cross-over Trials in Clinical Research.} Chichester: John Wiley & Sons; 2\ifelse{html}{\out{<sup>nd</sup>}}{\eqn{^{nd}}} edition 2002.

Zierhut ML, Bycott P, Gibbs MA, Smith BP, Vicini P. \emph{Ignorance is not bliss: Statistical power is not probability of trial success.} Clin Pharmacol Ther. 2015;99:356--9. \doi{10.1002/cpt.257}
}
\author{
B. Lang (thanks to G. Nehmiz for the helpful discussions), D. Labes
}
\seealso{
 \code{\link{expsampleN.TOST}, \link{power.TOST}} 
}
\examples{
# Expected power for a 2x2 crossover with 40 subjects
# CV 30\% known from a pilot 2x2 study with 12 subjects
# using all the defaults for other parameters (theta0 carved in stone)
exppower.TOST(CV = 0.3, n = 40, prior.parm = list(df = 12-2))
# should give: [1] 0.7365519
# or equivalently
exppower.TOST(CV = 0.3, n = 40, prior.parm = list(m = 12, design = "2x2"))

# In contrast: Julious approximation
exppower.TOST(CV = 0.3, n = 40, prior.parm = list(df = 10), method = "approx")
# should give: [1] 0.7359771

# Compare this to the usual (conditional) power (CV known, "carved in stone")
power.TOST(CV = 0.3, n = 40)
# should give: [1] 0.8158453
# same as if setting df = Inf in function exppower.TOST()
exppower.TOST(CV = 0.3, n = 40, prior.parm = list(df = Inf))

# Expected power for a 2x2 crossover with 40 subjects
# CV 30\% and theta0 = 0.95 known from a pilot 2x2 study with 12 subjects
# using uncertainty with respect to both CV and theta0
exppower.TOST(CV = 0.3, theta0 = 0.95, n = 40, 
              prior.parm = list(m = 12, design = "2x2"), prior.type = "both")
# should give [1] 0.5114685
}
