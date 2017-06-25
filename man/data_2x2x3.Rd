\encoding{utf-8}
\name{ct9.6.2+ct9.6.6}
\alias{data2x2x3}
\alias{ct9.6.2}
\alias{ct9.6.6}
\title{Sample size tables for the 2×2×3 replicate crossover}
\description{These data.frames give sample size tables calculated with 
\code{sampleN.TOST()} for the 2×2×3 replicate crossover design
(2-treatment 2-sequence 3-period design.}

\details{
After \code{data("data2x2x3")} the data.frames can be accessed by their names.\cr

ct9.6.2 is Table 9.6.2 from\cr
Chow SC, Liu JP. \emph{Design and Analysis of Bioavailability and Bioequivalence Studies.}\cr
Third edition, CRC Press, Chapman & Hall, Boca Raton (2009)\cr
Additive model, theta1=--0.2, theta2=+0.2 (BE limits 0.80 -- 1.20),\cr 
approximate power via shifted non-central \emph{t}-distribution.\cr

ct9.6.6 is Table 9.6.6 from the same reference.\cr
Multiplicative model, theta1=0.8, theta2=1.25 (1/theta1), 
power via shifted non-central \emph{t}-distribution.\cr
Attention! Chow and Liu’s CV is se (standard error) of residuals.
}
\note{Scripts for creation of these data.frames can be found in the \\test
sub-directory of the package.\cr
Comparing the results of that scripts to the corresponding data.frames can 
be used for validation purposes.}

\author{PowerTOST}

\examples{
data("data2x2x3")
ct9.6.2
ct9.6.6
}
