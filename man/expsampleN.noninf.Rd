\encoding{utf-8}
\name{expsampleN.noninf}
\alias{expsampleN.noninf}
\title{
Sample size based on expected power for the non-inferiority test
}

\description{
  Calculates the sample size based on the expected power for a variety
  of designs used in bioequivalence studies. See \link{known.designs}
  for the study designs covered.
}

\usage{
expsampleN.noninf(alpha = 0.025, targetpower = 0.8, logscale = TRUE, 
                  theta0, margin, CV, design = "2x2", robust = FALSE, 
                  prior.type = c("CV", "theta0", "both"), prior.parm = list(),
                  method = c("exact", "approx"), print = TRUE, details)
}

\arguments{
  \item{alpha}{
    Significance level (one-sided). Defaults here to 0.025.
  }
  \item{targetpower}{
    Power to achieve at least. Must be \verb{>0} and \verb{<1}.
    Typical values are \verb{0.8} or \verb{0.9}.
  }
  \item{logscale}{
    Should the data used on log-transformed or on original scale?
    \code{TRUE} or \code{FALSE}.\cr
    Defaults to \code{TRUE}.
  }  
  \item{theta0}{
    Assumed \sQuote{true} (or \sQuote{observed} in case of
    \code{prior.type != "CV"}) ratio or difference.\cr
    Typically set to \verb{0.95} (default if missing) if \code{logscale=TRUE}.
    Defaults to \verb{-0.05} if \code{logscale=FALSE}.
  }
  \item{margin}{
    Non-inferiority margin.\cr
    In case of \code{logscale=TRUE} it must be given as a ratio, otherwise
    as a difference.\cr
    Defaults to \verb{0.8} if \code{logscale=TRUE} or to \verb{-0.2} if
    \code{logscale=FALSE}.
  }
  \item{CV}{
    Assumed true or observed coefficient of variation as ratio (not
    percent). Only values \verb{> 0} are allowed.
    If \code{logscale=FALSE} CV is assumed to be the standard deviation.
    If \code{prior.type="CV"} may be given as vector: The CVs are then
    pooled (as a weighted mean with their degrees of freedom as weights).
  }
  \item{design}{
    Character string describing the study design.\cr
    See \code{known.designs()} for designs covered in this package.
  }
  \item{robust}{
    Defaults to \code{FALSE}. With that value the usual degrees of
    freedom will be used.\cr
    Setting to \code{TRUE} will use the degrees of freedom according to
    the \sQuote{robust} evaluation (aka Senn\enc{’}{'}s basic estimator).
    These df are calculated as \code{n-seq}.\cr
    See \code{known.designs()} for designs covered in this package.
  }
  \item{prior.type}{
    Specifies which parameter uncertainty should be accounted for.
    In case of \code{prior.type="CV"} (the default), only the
    uncertainty with respect to the CV will be considered (\emph{i.e.},
    the given treatment effect is assumed to be fix). In case of
    \code{prior.type="theta0"} only uncertainty with respect to the
    treatment ratio/difference will be accounted for (\emph{i.e.}, the
    given CV is assumed to be fix). In case of \code{prior.type="both"}
    the power value will be unconditional with respect to both the
    \code{CV} and \code{theta0}.
  }
  \item{prior.parm}{
    A list of parameters expressing the prior information about the 
    variability and/or treatment effect. Possible components are \code{df}, 
    \code{SEM}, \code{m}, \code{design}.\cr
    For \code{prior.type="CV"} the degrees of freedom from the prior trial
    are required. This information can be provided by specifying the
    single component\code{df} or the combination consisting of \code{m}
    and \code{design}.\cr
    For \code{prior.type = "theta0"} the standard error of the treatment difference 
    from the prior trial is required. This information can be provided by specifying
    the single component \code{SEM} or the combination consisting of 
    \code{m} and \code{design}.\cr
    For \code{prior.type = "both"} the degrees of freedom and the standard error of
    the treatment difference are required. This information can be provided by 
    specifying the combination consisting of \code{df} and \code{SEM} or via the 
    combination \code{m} and \code{design}.\cr
    See section \sQuote{Details} for a technical description of each
    component.
  }
  \item{method}{
    Defaults to \code{method="exact"}. In that case the expected power
    will be calculated as expected value of the power with respect to the
    (prior) distribution of the respective parameter(s).\cr
    Set to \code{method="approx"} the expected power according to the 
    approximate formulas given by Julious or Julious & Owen will be
    calculated (using the non-central \emph{t}); this only affects 
    \code{prior.type = "CV"}.
  }
  \item{print}{
    If \code{TRUE} (default) the function prints its results.\cr
    If \code{FALSE} only a data.frame with the results will be returned.
  }
  \item{details}{
    If \code{TRUE} the design characteristics and the steps during
    sample size calculations will be shown.\cr 
    If not specified, the default value is \code{FALSE} for
    \code{prior.type != "both"} and \code{TRUE} otherwise.
  }
}
\details{
  The sample size is calculated based on iterative evaluation of
  expected power. The starting value of the sample size search is
  taken from a large sample approximation if \code{prior.type="CV"}.
  Else an empirical start value is obtained. Note that in case of
  \code{prior.type="both"} the calculation may still take several seconds.
  
  Note also that the expected power is always bounded above by the
  so-called probability of technical success (\acronym{PTS}) which
  may be a value less than 1.Therefore, it may be possible that it
  is either not possible to calculate the required sample size at
  all or that the sample size gets very large if the given targetpower
  is less but close to the \acronym{PTS}.\cr
  
  \strong{Notes on the underlying hypotheses}\cr
  If the supplied margin is \verb{< 0} (\code{logscale=FALSE}) or
  \verb{< 1} (\code{logscale=TRUE}), then it is assumed \emph{higher} response
  values are better. The hypotheses are\cr
  \verb{  H0: theta0 <= margin}\cr
  \verb{  H1: theta0 >  margin}\cr
  where \code{theta0 = mean(test)-mean(reference)} if \code{logscale=FALSE}\cr
  or\cr
  \verb{  H0: log(theta0) <= log(margin)}\cr
  \verb{  H1: log(theta0) >  log(margin)}\cr
  where \code{theta0 = mean(test)/mean(reference)} if \code{logscale=TRUE}.\cr
  
  If the supplied margin is \verb{> 0} (\code{logscale=FALSE}) or
  \verb{> 1} (\code{logscale=TRUE}), then it is assumed \emph{lower} response
  values are better. The hypotheses are\cr
  \verb{  H0: theta0 >= margin}\cr
  \verb{  H1: theta0 <  margin}\cr
  where \code{theta0 = mean(test)-mean(reference)} if \code{logscale=FALSE}\cr
  or\cr
  \verb{  H0: log(theta0) >= log(margin)}\cr
  \verb{  H1: log(theta0) <  log(margin)}\cr
  where \code{theta0 = mean(test)/mean(reference)} if \code{logscale=TRUE}.\cr
  This latter case may also be considered as \sQuote{non-superiority}.
}

\value{
  A data.frame with the input values and the result of the sample
  size estimation.\cr
  The \code{Sample size} column contains the \emph{total} sample
  size in case of all designs implemented.
}

\author{
B. Lang, D. Labes
}

\references{
Grieve AP. \emph{Confidence Intervals and Sample Sizes.} Biometrics. 1991;47:1597--603. \doi{10.2307/2532411}

O’Hagan, Stevens, JW, Campell MJ. \emph{Assurance in Clinical Trial Design.} Pharm Stat. 2005;4:187--201. \doi{10.1002/pst.175}

Julious SA, Owen RJ. \emph{Sample size calculations for clinical studies allowing for uncertainty in variance.} Pharm Stat. 2006;5:29--37. \doi{10.1002/pst.197}

Julious SA. \emph{Sample sizes for Clinical Trials.} Boca Raton: CRC Press; 2010.

Bertsche A, Nehmitz G, Beyersmann J, Grieve AP. \emph{The predictive distribution of the residual variability in
the linear-fixed effects model for clinical cross-over trials.} Biom J. 2016;58(4):797--809. \doi{10.1002/bimj.201500245}

Box GEP, Tiao GC. \emph{Bayesian Inference in Statistical Analysis.} Boston: Addison-Wesley; 1992.

Held L, Sabanes Bove D. \emph{Applied Statistical Inference. Likelihood and Bayes.} Berlin, Heidelberg: Springer; 2014. \doi{10.1007/978-3-642-37887-4}

Senn S. \emph{Cross-over Trials in Clinical Research.} Chichester: John Wiley & Sons; 2\ifelse{html}{\out{<sup>nd</sup>}}{\eqn{^{nd}}} edition 2002.

Zierhut ML, Bycott P, Gibbs MA, Smith BP, Vicini P. \emph{Ignorance is not bliss: Statistical power is not probability of trial success.} Clin Pharmacol Ther. 2015;99:356--9. \doi{10.1002/cpt.257}
}

\seealso{
 \code{\link{exppower.noninf}, \link{known.designs}, \link{sampleN.noninf}}
}
\examples{
# Classical 2x2 cross-over, target power = 80\%,
# assumed true ratio = 95\%, margin = 0.8,
# intra-subject CV=30\% estimated from prior 2x2 trial 
# with m = 12 subjects
expsampleN.noninf(theta0 = 0.95, margin = 0.8, CV = 0.3, design = "2x2",
                  prior.parm = list(m = 12, design = "2x2"))
# gives n = 58 with achieved expected power 0.809148 
# Compare this to the usual sample size with CV assumed
# as 'carved in stone'
sampleN.noninf(theta0 = 0.95, margin = 0.8, CV = 0.3)

# Perform 'non-superiority' (lower is better) with assumed
# true ratio = 105\% and margin 125\%
expsampleN.noninf(theta0 = 1.05, margin = 1.25, CV = 0.3, design = "2x2",
                  prior.parm = list(m = 12, design = "2x2"))
# should give n = 56 with achieved expected power 0.806862

# More than one CV with corresponding degrees of freedom 
# other settings as above in first example
CVs <- c(0.25, 0.3)
dfs <- c(22, 10)
expsampleN.noninf(theta0 = 0.95, margin = 0.8, CV = CVs, 
                  prior.parm = list(df = dfs))
# should give a pooled CV=0.2664927 with 32 df and a sample
# size n=42 with achieved expected power 0.814073 exact
# achieved expected power 0.816163 approximate acc. to Julious

# Uncertainty is accounted for CV and theta0
\dontrun{
expsampleN.noninf(CV = 0.3, prior.type = "both",
                  prior.parm = list(m = 12, design = "2x2"))
# gives a dramatic increase in sample size (n = 194)
# due to small pilot trial}
}
