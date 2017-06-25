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

ct5.1 is Table 5.1 from\cr
Hauschke D, Steinijans V, Pigeot I. \emph{Bioequivalence studies in Drug Development.}\cr
John Wiley & Sons, Chichester (2007)\cr
Multiplicative model, theta1=0.8, theta2=1.25 (1/theta1), exact\cr

ct5.2 is Table 5.2 from the same source\cr
Multiplicative model, theta1=0.75, theta2=1.3333 (1/theta1), exact\cr

ct5.3 is Table 5.3 from the same source\cr
Multiplicative model, theta1=0.9, theta2=1.1111 (1/theta1), exact\cr

ct5.4.1 is Table 5.4.1 from\cr
Chow SC, Liu JP. \emph{Design and Analysis of Bioavailability and Bioequivalence Studies.}\cr
Third edition, CRC Press, Chapman & Hall, Boca Raton (2009)\cr
Additive model, theta1=--0.2, theta2=+0.2 (BE limits 0.80 -- 1.20), exact\cr
}
\note{Scripts for creation of these data.frames can be found in the \\test
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
