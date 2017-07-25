\name{expsampleN.TOST}
\alias{expsampleN.TOST}
\title{
Sample size based on expected power
}
\description{
Calculates the sample size based on the expected power for a variety of study 
designs used in bioequivalence studies. See known.designs() for the study 
designs covered.
}
\usage{
expsampleN.TOST(alpha = 0.05, targetpower = 0.8, logscale=TRUE, theta0, 
                theta1, theta2, CV, design = "2x2", robust = FALSE, 
                prior.type = c("CV", "theta0", "both"), prior.parm = list(),
                method = c("exact", "approx"), print = TRUE, details)
}
\arguments{
  \item{alpha}{
Significance level (one-sided). Commonly set to 0.05.
  }
  \item{targetpower}{
Power to achieve at least. Must be >0 and <1.
Typical values are 0.8 or 0.9.
  }
  \item{logscale}{
Should the data used on log-transformed or on original scale? TRUE or FALSE.\cr
Defaults to TRUE.
  }  
  \item{theta0}{
Assumed \sQuote{true} (or \sQuote{observed} in case of \code{prior.type != "CV"}) bioequivalence 
ratio or difference.\cr
Typically set to 0.95 (default if missing) if \code{logscale=TRUE}.
Defaults to 0.05 if \code{logscale=FALSE}.
  }
  \item{theta1}{
Lower bioequivalence limit as ratio if \code{logscale=TRUE} or as difference.\cr
Can be missing. Defaults then to 0.8 if \code{logscale=TRUE} or to -0.2 if 
\code{logscale=FALSE}.
  }
  \item{theta2}{
Upper bioequivalence limit as ratio if \code{logscale=TRUE} or as difference.\cr
If not given \code{theta2} will be calculated as \code{1/theta1} if \code{logscale=TRUE},
else as \code{-theta1}.
  }
  \item{CV}{
Assumed true or observed coefficient of variation as ratio. Only values > 0 are allowed.\cr
If \code{prior.type="CV"} may be given as vector: The CVs are then pooled (as a 
weighted mean with their degrees of freedoms as weights).
  }
  \item{design}{
Character string describing the study design.\cr
See \code{known.designs()} for designs covered in this package.
  }
  \item{robust}{
Defaults to FALSE. With that value the usual degrees of freedom will be used.\cr
Set to \code{TRUE} will use the degrees of freedom according to the \sQuote{robust} evaluation
(aka Senn's basic estimator). These df are calculated as \code{n-seq}.\cr
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
\code{SEM}, \code{m}, \code{design}.\cr
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
paper will be calculated (using non-central \emph{t}); this only affects 
\code{prior.type = "CV"}.
}
  \item{print}{
If TRUE (default) the function prints its results.\cr
If FALSE only a data.frame with the results will be returned.
  }
  \item{details}{
If TRUE the design characteristics and the steps during
sample size calculations will be shown.\cr 
If not specified, the default value is FALSE for \code{prior.type != "both"}
and TRUE otherwise.
  }
}
\details{
The sample size is calculated based on iterative evaluation of expected power.

The starting value of the sample size search is taken from a large sample 
approximation if \code{prior.type = "CV"}. Else an empirical start value is
obtained.

Note that in case of \code{prior.type = "both"} the calculation may still take several seconds.

Note also that the expected power is always bounded above by the so-called
probability of technical success (\acronym{PTS}) which may be a value less than 1. 
Therefore, it may be possible that it is either not possible to calculate the 
required sample size at all or that the sample size gets very large 
if the given targetpower is less but close to the \acronym{PTS}.
}
\value{
A data.frame with the input values and the result of the sample size estimation.\cr
The "Sample size" column contains the \bold{total} sample size in case of all
designs implemented.
}
\references{
Grieve AP.\cr 
\emph{Confidence Intervals and Sample Sizes}\cr
Biometrics. 1991;47:1597--603. \doi{10.2307/2532411}

O'Hagan et al.\cr 
\emph{Assurance in Clinical Trial Design}\cr 
Pharm Stat. 2005;4:187--201. \doi{10.1002/pst.175}

Julious SA, Owen RJ.\cr 
\emph{Sample size calculations for clinical studies allowing for 
uncertainty in variance} \cr
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
 \code{\link{exppower.TOST}, \link{known.designs}, \link{sampleN.TOST}}
}
\examples{
# Classical 2x2 cross-over, target power = 80\%,
# BE limits 80 ... 125\%, assumed true BE ratio = 95\%,
# intra-subject CV=30\% estimated from prior 2x2 trial 
# with m = 30 subjects
expsampleN.TOST(CV=0.3, prior.parm = list(m = 30, design = "2x2"))
# -> gives n = 42 with achieved expected power 0.806262
# Compare this to the usual sample size with CV assumed known ('carved in stone')
sampleN.TOST(CV=0.3)
# -> gives n = 40 subjects
# Compare this to the case where uncertainty is accounted for CV and theta0
# Not run due to timing policy of CRAN - may run several seconds
\dontrun{
expsampleN.TOST(CV=0.3, prior.parm = list(m = 30, design = "2x2"), 
                prior.type = "both")}
# -> gives n = 72 subjects

# More than one CV with corresponding degrees of freedom 
# other settings as above in first example
CVs <- c(0.25, 0.3)
dfs <- c(22, 10)
expsampleN.TOST(CV=CVs, prior.parm = list(df = dfs))
# -> gives a pooled CV=0.2664927 with df=32
# and a sample size n=34 with achieved expected power 0.812653 exact
# achieved expected power 0.815019 approximate acc. Julious
}
