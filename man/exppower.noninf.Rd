\name{exppower.noninf}
\alias{exppower.noninf}
\title{
Expected power of the non-inferiority test
}
\description{
Calculates the so-called expected, i.e. unconditional, power for a variety of
study designs used in bioequivalence studies.
}
\usage{
exppower.noninf(alpha = 0.025, logscale = TRUE, theta0, margin, CV, n,
                design = "2x2", robust = FALSE, dfCV, 
                prior.type = c("CV", "theta0", "both"), prior.parm = list(), 
                method = c("exact", "approx"))
}
\arguments{
  \item{alpha}{
Significance level (one-sided). Defaults here to 0.025.
  }
  \item{logscale}{
Should the data be used on log-transformed or on original scale? TRUE or FALSE.\cr
Defaults to TRUE.
  }  
  \item{theta0}{
Assumed 'true' (or 'observed' in case of \code{prior.type != "CV"}) ratio 
or difference.\cr
Typically set to 0.95 (default if missing) if \code{logscale=TRUE}.
Defaults to -0.05 if \code{logscale=FALSE}.
  }
  \item{margin}{
Non-inferiority margin.\cr
In case of \code{logscale=TRUE} it must be given as ratio, otherwise as
difference.\cr
Defaults to 0.8 if \code{logscale=TRUE} or to -0.2 if \code{logscale=FALSE}.
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
See \code{known.designs()} for designs covered in this package.
}
  \item{robust}{
Defaults to FALSE. Set to \code{TRUE} will use the degrees of freedom according 
to the 'robust' evaluation (aka Senn's basic estimator). These df are calculated
as \code{n-seq}.\cr 
See \code{known.designs()$df2} for designs covered in this package.
}
  \item{dfCV}{
Alias of the component \code{df} in argument \code{prior.parm} for backward
compatibility.\cr
This argument is deprecated and will be removed in next versions.
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
For \code{prior.type != "CV"} the degrees of freedom and the standard error of
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
paper will be calculated (using non-central t); this only affects 
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
coincide with \code{power.noninf(...)}.\cr
Note: This corresponds to the df of both the CV and the difference of means.}
\item{\code{SEM}}{Standard error of the difference of means from the prior trial; 
must always be on additive scale (i.e. usually log-scale).}
\item{\code{m}}{Number of subjects from prior trial. Specification is analogous to
the main argument \code{n}.}
\item{\code{design}}{Study design of prior trial. Specification is analogous to the 
main argument \code{design}.}
}
For \code{prior.parm}, the combination consisting of \code{df} and \code{SEM} 
requires a somewhat advanced knowledge of the prior trial (the standard error of
the difference usually needs separate calculation or is only provided in the
raw output, for example coming from the software SAS). However, it has the advantage
that if the prior trial had missing data the exact degrees of freedom
and standard error can be used, the former being possibly non-integer valued.
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
\cr
Notes on the underlying hypotheses:\cr\cr
If the supplied margin is < 0 (\code{logscale=FALSE}) or < 1 (\code{logscale=TRUE}), 
then it is assumed higher response values are better. The hypotheses are\cr
\code{H0: theta0 <= margin vs. H1: theta0 > margin}\cr
where \code{theta0 = mean(test)-mean(reference)} if \code{logscale=FALSE}\cr
or\cr
\code{H0: log(theta0) <= log(margin) vs. H1: log(theta0) > log(margin)}\cr
where \code{theta0 = mean(test)/mean(reference)} if \code{logscale=TRUE}.
\cr\cr
If the supplied margin is > 0 (\code{logscale=FALSE}) or > 1 (\code{logscale=TRUE}), 
then it is assumed lower response values are better. The hypotheses are\cr
\code{H0: theta0 >= margin vs. H1: theta0 < margin}\cr
where \code{theta0 = mean(test)-mean(reference)} if \code{logscale=FALSE}\cr
or\cr
\code{H0: log(theta0) >= log(margin) vs. H1: log(theta0) < log(margin)}\cr
where \code{theta0 = mean(test)/mean(reference)} if \code{logscale=TRUE}.\cr
This latter case may also be considered as 'non-superiority'.
}
\value{
Value of expected power according to the input.
}
\references{
Grieve AP.\cr 
\emph{Confidence Intervals and Sample Sizes}\cr
Biometrics. 1991;47:1597--603. \doi{10.2307/2532411}

O'Hagan et al.\cr 
\emph{Assurance in Clinical Trial Design}\cr
Pharm Stat. 2005;4:187--201. \doi{10.1002/pst.175}

Julious SA, Owen RJ.\cr
\emph{Sample size calculations for clinical studies allowing for uncertainty in variance}\cr
Pharm Stat. 2006;5:29--37. \doi{10.1002/pst.197}

Julious SA.\cr 
\emph{Sample sizes for Clinical Trials}\cr
Boca Raton: CRC Press / Chapman & Hall; 2010.

Bertsche et al.\cr
\emph{The predictive distribution of the residual variability in
the linear-fixed effects model for clinical cross-over trials}\cr
Biom J. 2016;58(4):797--809. \doi{10.1002/bimj.201500245}

Box GEP, Tiao GC.\cr 
\emph{Bayesian Inference in Statistical Analysis}\cr
Boston: Addison-Wesley; 1992.

Held L, Sabanes Bove D.\cr
\emph{Applied Statistical Inference. Likelihood and Bayes}\cr
Berlin, Heidelberg: Springer; 2014. \doi{10.1007/978-3-642-37887-4}

Senn S.\cr
\emph{Cross-over Trials in Clinical Research}\cr
Chichester: Wiley; 2002.

Zierhut ML et al.\cr
\emph{Ignorance is not bliss: Statistical power is not probability of trial success}\cr
Clin Pharmacol Ther. 2015;99:356--9. \doi{10.1002/cpt.257}
}
\author{
B. Lang & D. Labes
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
