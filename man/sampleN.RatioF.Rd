\encoding{utf-8}
\name{sampleN.RatioF}
\alias{sampleN.RatioF}
\title{
Sample size for equivalence of the ratio of two means with normality on original scale
}
\description{
Estimates the necessary sample size to have at least a given power based on 
Fieller’s confidence (\sQuote{fiducial}) interval.
}
\usage{
sampleN.RatioF(alpha = 0.025, targetpower = 0.8, theta1 = 0.8, theta2, 
               theta0 = 0.95, CV, CVb, design = "2x2", print = TRUE,
               details = FALSE, imax=100, setseed=TRUE)
}
\arguments{
  \item{alpha}{
Type I error probability. \cr
Defaults here to 0.025 because this function is intended for studies
with clinical endpoints. 
  }
  \item{targetpower}{
Power to achieve at least. Must be >0 and <1. Typical values are 0.8 or 0.9.
  }
  \item{theta1}{
Lower bioequivalence limit. Typically 0.8 (default).
  }
  \item{theta2}{
Upper bioequivalence limit. Typically 1.25.\cr
Is set to \code{1/theta1} if missing.
  }
  \item{theta0}{
\sQuote{True} or assumed T/R ratio. Typically set to 0.95.
  }
  \item{CV}{
Coefficient of variation as ratio. In case of \code{design="parallel"} this is
the CV of the total variability, in case of \code{design="2x2"} the intra-subject
CV (CVw in the reference).
  }
  \item{CVb}{
CV of the between-subject variability. Only necessary for \code{design="2x2"}.
  }
  \item{design}{
A character string describing the study design.\cr
\code{design="parallel"} or \code{design="2x2"} allowed for a two-parallel 
group design or a classical TR|RT crossover design.
  }
  \item{print}{
If \code{TRUE} (default) the function prints its results. If \code{FALSE} only
a data.frame with the results will be returned. 
  }
  \item{details}{
If \code{TRUE} the steps during sample size calculations will be shown.\cr 
Defaults to \code{FALSE}.
  }
  \item{imax}{
Maximum number of steps in sample size search.\cr 
Defaults to 100. Adaption only in rare cases needed.
  }
  \item{setseed}{
If set to \code{TRUE} the dependence of the power from the state of the random number
generator is avoided.
  }
}
\details{
The sample size is based on exact power calculated using the bivariate 
non-central \emph{t}-distribution via function \code{\link[mvtnorm]{pmvt}} of the package \code{mvtnorm}.\cr
Due to the calculation method used in package mvtnorm these
probabilities are dependent from the state of the random number generator
within the precision of the power.\cr\cr
The CV(within) and CVb(etween) in case of \code{design="2x2"} are obtained 
via an appropriate ANOVA from the error term and from the difference 
\code{(MS(subject within sequence)-MS(error))/2}.\cr\cr
The estimated sample size gives always the \emph{total} number of subjects (not subject/sequence in crossovers or subjects/group in parallel designs -- like in some other software packages).
}
\value{
A data.frame with the input values and results will be returned.\cr
The sample size n returned is the total sample size for both designs.\cr
}
\references{
Fieller EC. \emph{Some Problems in Interval Estimation.} J Royal Stat Soc B. 1954;16(2):175--85. \doi{10.1111/j.2517-6161.1954.tb00159.x}

Sasabuchi S. \emph{A test of a multivariate normal mean with composite hypotheses determined by linear inequalities.} Biometrika. 1980;67(2):429--39. \doi{10.1093/biomet/67.2.429}

Hauschke D, Kieser M, Diletti E, Burke M. \emph{Sample size determination for proving equivalence based on the ratio of two means for normally distributed data.} Stat Med. 1999;18(1):93--105.

Hauschke D, Steinijans V, Pigeot I. \emph{Bioequivalence Studies in Drug Development.} Chichester: Wiley; 2007. Chapter 10.

European Agency for the Evaluation of Medicinal Products, CPMP. \emph{Points to Consider on Switching between Superiority and Non-Inferiority.} London, 27 July 2000. \href{https://www.ema.europa.eu/en/documents/scientific-guideline/points-consider-switching-between-superiority-non-inferiority_en.pdf}{CPMP/EWP/482/99}

}
\author{D. Labes}

\note{
This function is intended for studies with clinical endpoints.\cr In such studies
the 95\% confidence intervals are usually used for equivalence testing.\cr
Therefore, alpha defaults here to 0.025 (see EMEA 2000).
}

\seealso{
\code{\link{power.RatioF}}
}

\examples{
# sample size for a 2x2 cross-over study
# with CVw=0.2, CVb=0.4
# alpha=0.025 (95\% CIs), target power = 80\%
# 'true' ratio = 95\%, BE acceptance limits 80-125\%
# using all the defaults:
sampleN.RatioF(CV = 0.2, CVb = 0.4)
# gives n=28 with an achieved power of 0.807774
# see Hauschke et.al. (2007) Table 10.3a

# sample size for a 2-group parallel study
# with CV=0.4 (total variability) 
# alpha=0.025 (95\% CIs), target power = 90\%
# 'true' ratio = 90\%, BE acceptance limits 75-133.33\%
sampleN.RatioF(targetpower = 0.9, theta1 = 0.75,
               theta0 = 0.90, CV = 0.4, design = "parallel")
# gives n=236 with an achieved power of 0.900685
# see Hauschke et.al. (2007) Table 10.2

# a rather strange setting of ratio0! have a look at n.
# it would be better this is not the sample size but your account balance ;-).
sampleN.RatioF(theta0 = 0.801, CV = 0.2, CVb = 0.4)

}
