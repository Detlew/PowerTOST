\encoding{utf-8}
\name{CVwRfromU}
\alias{CVwRfromU}
\alias{U2CVwR}
\title{CVwR from the upper expanded limit (ABEL)}

\description{
  Calculates the intra-subject CV (coefficient of variation) of the reference from
  the upper expanded limit of a BE study (replicate design for ABEL).
  Useful if no
  \ifelse{html}{\out{<em>CV</em><sub>wR</sub>}}{\eqn{CV_{wR}}{CVwR}} but the
  expanded limits were given.
}

\usage{
CVwRfromU(U, regulator = "EMA")
U2CVwR(U, regulator = "EMA")
}

\arguments{
  \item{U}{
    Upper expanded limit.\cr
    Must be within \{1.2500, 1.4319\} if \code{regulator="EMA"} and within \{1.2500, 1.5000\}
    if \code{regulator="HC"}.
  }
  \item{regulator}{
    Regulatory body\enc{’}{'}s settings for expanding the BE acceptance limits,
    given as a string from the choices \code{"EMA"} or \code{"HC"}. 
    Defaults to \code{regulator="EMA"}.
  }
}

\details{
  Only the upper expanded limit is supported since it offers one more 
  significant digit than the lower expanded limit.
}

\value{
  Numeric value of the \code{CVwR} as ratio, where
  \code{CVwR = sqrt(exp((log(U)/r_const)^2)-1)}.
}

\note{
  \code{U2CVwR()} is simply an alias to \code{CVwRfromU()}.
}

\author{
  H. \enc{Schütz}{Schuetz}
}

\examples{
# Given the upper expanded limit and using the defaults
CVwRfromU(U = 1.38)
# should give [1] 0.44355, i.e., a CVwR ~ 44\%
# Upper limit from a study according the Health Canada’s rules
CVwRfromU(U = 1.48, regulator = "HC")
# should give [1] 0.55214
}