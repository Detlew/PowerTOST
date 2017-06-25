\encoding{utf-8}
\name{ctSJ.VIII.10+ctSJ.VIII.20+ctCW.III}
\alias{data_parallel}
\alias{ctSJ.VIII.10}
\alias{ctSJ.VIII.20}
\alias{ctCW.III}
\title{Sample size tables for the parallel group design}
\description{These data.frames give sample size tables calculated with 
\code{sampleN.TOST()} for the parallel group design (2 groups).}
\details{
After \code{data("data_parallel")} the data.frames can be accessed by their names.\cr

ctSJ.VIII.10 is Table VIII, column 'level of bioequivalence 10\%' from\cr
Julious SA. \emph{Tutorial in Biostatistics. Sample sizes for clinical trials with Normal data.}\cr
Stat Med. 2004;23(12):1921--86. \doi{10.1002/sim.1783}

Multiplicative model, theta1=0.9, theta2=1.1111 (1/theta1), target power=90\%,\cr
power approximate via non-central \emph{t}-distribution.\cr
Attention! Julious gives sample size per group.

ctSJ.VIII.20 is Table VIII from the same source\cr
column 'level of bioequivalence 20\%'\cr
Multiplicative model, theta1=0.8, theta2=1.25 (1/theta1), target power=90\%,\cr
power approximate via non-central \emph{t}.\cr

ctCW.III is Table III from\cr
Chow SC, Wang H. \emph{On Sample Size Calculation in Bioequivalence Trials}\cr
J Pharmacokinet Pharmacodyn. 2001. 28(2):155--69.\cr
Additive model, theta1=--0.2, theta2=+0.2 (BE limits 0.80 -- 1.20), exact.\cr

Seems the last reference is not very reliable (compare to the Table in the paper).
}
\note{Scripts for creation of these data.frames can be found in the \\test
sub-directory of the package.\cr
Comparing the results of that scripts to the corresponding data.frames can 
be used for validation purposes.}

\author{PowerTOST}

\examples{
data("data_parallel")
ctSJ.VIII.10
ctSJ.VIII.20
ctCW.III
}
