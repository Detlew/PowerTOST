\encoding{utf-8}
\name{CI.BE}
\alias{CI.BE}
\title{
1--2*alpha confidence interval given point estimate, CV, and n
}

\description{
  Utility function to calculate the
  \ifelse{html}{\out{1&ndash;2&alpha;}}{\eqn{1-2\alpha}{1-2*alpha}}
  CI given point estimate, CV, and n for the various designs
  covered in this package.
}

\usage{
CI.BE(alpha = 0.05, pe, CV, n, design = "2x2", robust = FALSE)
}

\arguments{
  \item{alpha}{
    Type I error probability, significance level. Defaults to 0.05.
  }
  \item{pe}{
    Point estimate (GMR).
  }
  \item{CV}{
    Coefficient of variation as ratio (not percent). 
  }
  \item{n}{
    Total number of subjects if a scalar is given.\cr
    Number of subjects in (sequence) groups if given as vector.
  }
  \item{design}{
    Character string describing the \enc{study’s}{study's} design.\cr
    See \code{known.designs()} for designs covered in this package.
  }
  \item{robust}{
    Defaults to \code{FALSE}.\cr
    Setting to \code{TRUE} will use the degrees of freedom according
    to the \sQuote{robust} evaluation (aka \enc{Senn’s}{Senn's} basic
    estimator). These degrees of freedom are calculated as \code{n-seq}.\cr
    See \code{known.designs()$df2} for designs covered in this package.
  }
}

\value{
  Returns the \ifelse{html}{\out{1&ndash;2&alpha;}}{\eqn{1-2\alpha}{1-2*alpha}}
  confidence interval.\cr
  Returns a vector with named elements \code{lower}, \code{upper} if
  arguments \code{pe} and \code{CV} are scalars, else a matrix with
  columns \code{lower}, \code{upper} is returned.
}

\note{
  The function assumes an evaluation using log-transformed data.\cr
  The function assumes equal variances in case of \code{design="parallel"}
  and the higher order crossover designs.\cr
  The implemented formula covers balanced and unbalanced designs.\cr\cr
  Whether the function vectorizes properly is not thoroughly tested.
}

\author{
D. Labes
}

\examples{
# 90\% confidence interval for the 2x2 crossover
# n(total) = 24
CI.BE(pe = 0.95, CV = 0.3, n = 24)
# should give
#     lower     upper 
# 0.8213465 1.0988055 
# same total number but unequal sequences
CI.BE(pe = 0.95, CV = 0.3, n = c(13, 11))
#     lower     upper
# 0.8209294 1.0993637   
}
