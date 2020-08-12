\encoding{utf-8}
\name{CI.RatioF}
\alias{CI.RatioF}
\title{
1--2*alpha Fieller CI given point estimate, CV (, CVb) and n
}

\description{
  Utility function to calculate the
  \ifelse{html}{\out{1&ndash;2&alpha;}}{\eqn{1-2\alpha}{1-2*alpha}}
  Fieller confidence interval given the point estimate, 
  CV (, CVb), and n for the parallel group and
  \ifelse{html}{\out{2&times;2}}{\eqn{2 \times 2}{2x2}}
  crossover.
}

\usage{
CI.RatioF(alpha = 0.025, pe, CV, CVb, n, design = c("2x2", "parallel"))
}

\arguments{
  \item{alpha}{
    Type I error probability, aka significance level.\cr
    Defaults here to 0.025 because this function is intended for studies
    with clinical endpoints. 
  }
  \item{pe}{
    Point estimate of T/R ratio.
  }
  \item{CV}{
    Coefficient of variation as ratio (not percent). In case of \code{design="parallel"}
    this is the CV of the total variability, in case of \code{design="2x2"}
    the intra-subject CV.
  }
  \item{CVb}{
    CV of the between-subject variability. Only necessary for
    \code{design="2x2"}.
  }
  \item{n}{
    Total number of subjects if a scalar is given.\cr
    Number of subjects in (sequence) groups if given as vector.
  }
  \item{design}{
    A character string describing the study design.\cr 
    \code{design="parallel"} or \code{design="2x2"} allowed for a
    parallel two-group design or a classical TR|RT crossover design.
  }
}

\details{
  The CV(within) and CVb(etween) in case of \code{design="2x2"} are
  obtainedvia an appropriate ANOVA from the error term and from the
  difference \code{(MS(subject within sequence)-MS(error))/2}.
}

\value{
  Returns the
  \ifelse{html}{\out{1&ndash;2&alpha;}}{\eqn{1-2\alpha}{1-2*alpha}}
  confidence interval.
}

\note{
  The function assumes an evaluation using un-transformed data.\cr
  The function assumes equal variances in case of \code{design="parallel"}.\cr
  The formula implemented covers balanced and unbalanced designs.\cr\cr
  Note that when the mean of the denominator of the ratio is close to zero, 
  confidence intervals might be degenerated and are returned as \code{NA}. In such
  a case a warning is issued.\cr\cr
  Whether the function vectorizes properly is not thoroughly tested.\cr\cr
  This function is intended for studies with clinical endpoints. In such
  studies the 95\% confidence intervals are usually used for equivalence 
  testing. Therefore, alpha defaults here to 0.025 (see \acronym{EMA} 2000).
}

\author{
  D. Labes
}

\references{
Locke CS. \emph{An exact confidence interval from untransformed data for the ratio of 
two formulation means.} J Pharmacokin Biopharm. 1984;12(6):649--55. \doi{10.1007/BF01059558}

Hauschke D, Steinijans VW, Pigeot I. \emph{Bioequivalence Studies in Drug Development.} Chichester: John Wiley; 2007. Chapter 10. \doi{10.1002/9780470094778.fmatter}

European Medicines Agency, Committee for Proprietary Medicinal Products. \emph{Points to consider on switching between superiority and non-inferiority.} London, 27 July 2000. \href{https://www.ema.europa.eu/documents/scientific-guideline/points-consider-switching-between-superiority-non-inferiority_en.pdf}{CPMP/EWP/482/99}
}

\seealso{
\code{\link{CI.BE}}, \code{\link{power.RatioF}}
}

\examples{
# 95\% Fieller CI for the 2x2 crossover
CI.RatioF(pe = 0.95, CV = 0.3, CVb = 0.6, n = 32)
}
