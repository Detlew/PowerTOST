\encoding{utf-8}
\name{ct5.1+ct5.2+ct5.3+ct5.4.1}
\alias{data2x2}
\alias{ct5.1}
\alias{ct5.2}
\alias{ct5.3}
\alias{ct5.4.1}
\title{Sample Size Tables for the Classical 2x2 Crossover Design}
\description{These data.frames give sample size tables calculated with 
\code{sampleN.TOST()} for the 2×2 design.}

\details{
The data.frames can be accessed by their names.\cr

  \tabular{cl}{
    data.frame \tab Description\cr
    ct5.1 \tab Multiplicative model, theta1=0.8, theta2=1.25 (1/theta1), exact\cr
    ct5.2 \tab Multiplicative model, theta1=0.75, theta2=1.3333 (1/theta1), exact\cr
    ct5.3 \tab Multiplicative model, theta1=0.9, theta2=1.1111 (1/theta1), exact\cr
    ct5.4.1 \tab Additive model, theta1=--0.2, theta2=+0.2 (BE limits 0.80 -- 1.20), exact
  }
}

\source{
  \tabular{cll}{
    data.frame \tab Origin \tab Details\cr
    ct5.1 \tab Hauschke \emph{et al.} \tab Table 5.1 (p 113--114)\cr
    ct5.2 \tab Hauschke \emph{et al.} \tab Table 5.2 (p 115--116) \cr
    ct5.3 \tab Hauschke \emph{et al.} \tab Table 5.3 (p 118)\cr
    ct5.4.1 \tab Chow & Liu \tab Table 5.4.1 (p 158)
  }
}

\references{
Hauschke D, Steinijans VW, Pigeot I. \emph{Bioequivalence Studies in Drug Development.} Chichester: John Wiley; 2007.

Chow SC, Liu JP. \emph{Design and Analysis of Bioavailability and Bioequivalence Studies.} Boca Raton: CRC Press; 3\ifelse{html}{\out{<sup>rd</sup>}}{\eqn{^{rd}}} edition 2009.
}

\note{Scripts for creation of these data.frames can be found in the \code{/test}
sub-directory of the package.\cr
Comparing the results of these scripts to the corresponding data.frames can 
be used for validation purposes.}

\author{PowerTOST}

\examples{
ct5.1
ct5.2
ct5.3
ct5.4.1
}
