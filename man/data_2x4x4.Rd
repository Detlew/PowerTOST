\encoding{utf-8}
\name{ct9.6.4+ct9.6.8}
\alias{data2x4x4}
\alias{ct9.6.4}
\alias{ct9.6.8}
\title{Sample Size Tables for the 2x4x4 Replicate Crossover Design}
\description{These data.frames give sample size tables calculated with 
\code{sampleN.TOST()} for the \enc{2×4×4}{2x4x4} replicate crossover design
(2-treatment 4-sequence 4-period design).}

\details{
The data.frames can be accessed by their names.\cr

  \tabular{ll}{
    data.frame \tab Description\cr
    ct9.6.4 \tab Additive model, theta1=--0.2, theta2=+0.2 (BE limits 0.80 -- 1.20)\cr
    \tab approximate power via shifted non-central \emph{t}-distribution\cr
    ct9.6.8 \tab Multiplicative model, theta1=0.8, theta2=1.25 (1/theta1)\cr
    \tab approximate power via shifted non-central \emph{t}-distribution
  }
Attention! CV is se (standard error) of residuals.
}

\source{
  \tabular{lcl}{
    data.frame \tab Origin \tab Details\cr
    ct9.6.4 \tab Chow & Liu \tab Table 9.6.4 (p 294)\cr
    ct9.6.8 \tab Chow & Liu \tab Table 9.6.8 (p 298)
  }
}

\references{
Chow SC, Liu JP. \emph{Design and Analysis of Bioavailability and Bioequivalence Studies.} Boca Raton: CRC Press; 3\ifelse{html}{\out{<sup>rd</sup>}}{\eqn{^{rd}}} edition 2009.
}

\note{Scripts for creation of these data.frames can be found in the \code{/test}
sub-directory of the package.\cr
Comparing the results of these scripts to the corresponding data.frames can 
be used for validation purposes.}

\author{PowerTOST}

\examples{
ct9.6.4
ct9.6.8
}
