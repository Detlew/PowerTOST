\name{ct9.6.4+ct9.6.8}
\alias{data2x4x4}
\alias{ct9.6.4}
\alias{ct9.6.8}
\title{Sample size tables for the 2x4x4 replicate crossover}
\description{These data.frames give sample size tables calculated with 
\code{sampleN.TOST()} for the 2x4x4 replicate crossover design
(2-treatment-4-sequence-4-period design).}

\details{
The data.frame's can be accessed by their names or by \code{data("name")}.\cr

ct9.6.4 is Table 9.6.4 from\cr
Chow S.C., Liu J.P.\cr
"Design and Analysis of Bioavailability and Bioequivalence Studies",\cr
Third edition, CRC Press, Chapman & Hall, Boca Raton (2009)\cr
Additive model, theta1=-0.2, theta2=+0.2 (BE limits 0.80 - 1.20),\cr 
approximate power via shifted non-central t-distribution.\cr

ct9.6.8 is Table 9.6.8 from the same reference.\cr
Multiplicative model, theta1=0.8, theta2=1.25 (1/theta1), 
power via shifted non-central t-distribution.\cr
Attention! Chow and Liu's CV in case of multiplicative model is se 
(standard error) of residuals.
}
\note{Scripts for creation of these data.frame's can be found in the \\test
sub-directory of the package.\cr
Comparing the results of that scripts to the corresponding data.frames can 
be used for validation purposes.}

\author{PowerTOST}
