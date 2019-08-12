\encoding{utf-8}
\name{ctSJ.VIII.10+ctSJ.VIII.20+ctCW.III}
\alias{data_parallel}
\alias{ctSJ.VIII.10}
\alias{ctSJ.VIII.20}
\alias{ctCW.III}
\title{Sample Size Tables for the Parallel Group Design}
\description{These data.frames give sample size tables calculated with 
\code{sampleN.TOST()} for the parallel group design (2 groups).}
\details{
The data.frames can be accessed by their names.\cr

  \tabular{ll}{
    data.frame \tab Description\cr
    ctSJ.VIII.10 \tab Multiplicative model, theta1=0.9, theta2=1.1111 (1/theta1), target power=90\%\cr
    \tab approximate power via non-central \emph{t}-distribution\cr
    ctSJ.VIII.20 \tab Multiplicative model, theta1=0.8, theta2=1.25 (1/theta1), target power=90\%\cr
    \tab approximate power via non-central \emph{t}-distribution\cr
    ctCW.III \tab Additive model, theta1=--0.2, theta2=+0.2 (BE limits 0.80 -- 1.20), exact
  }
Attention! Julious gives sample size per group.
}

\source{
  \tabular{lcl}{
    data.frame \tab Origin \tab Details\cr
    ctSJ.VIII.10 \tab Julious \tab Table VIII (p. 1972), column \sQuote{Level of bioequivalence 10\%}\cr
    ctSJ.VIII.20 \tab Julious \tab Table VIII (p. 1972), column \sQuote{Level of bioequivalence 20\%}\cr
    ctCW.III \tab Chow & Wang \tab Table III (p. 164)
  }
Seems the last reference is not very reliable (compare to the table in the paper).
}

\references{
Julious SA. \emph{Tutorial in Biostatistics. Sample sizes for clinical trials with Normal data.} Stat Med. 2004;23(12):1921--86. \doi{10.1002/sim.1783}

Chow SC, Wang H. \emph{On Sample Size Calculation in Bioequivalence Trials.} J Pharmacokinet Pharmacodyn. 2001;28(2):155--69. \doi{10.1023/A:1011503032353}
}

\note{Scripts for creation of these data.frames can be found in the \code{/test}
sub-directory of the package.\cr
Comparing the results of these scripts to the corresponding data.frames can 
be used for validation purposes.}

\author{PowerTOST}

\examples{
ctSJ.VIII.10
ctSJ.VIII.20
ctCW.III
}
