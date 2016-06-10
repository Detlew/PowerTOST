\name{ctSJ.VIII.10+ctSJ.VIII.20+ctCW.III}
\alias{data_parallel}
\alias{ctSJ.VIII.10}
\alias{ctSJ.VIII.20}
\alias{ctCW.III}
\title{Sample size tables for the parallel group design}
\description{These data.frames give sample size tables calculated with 
\code{sampleN.TOST()} for the parallel group design (2 groups).}
\details{
The data.frame's can be accessed by their names or by \code{data("name")}.\cr

ctSJ.VIII.10 is Table VIII, column 'level of bioequivalence 10\%' from\cr
S.A.Julious\cr
"Tutorial in Biostatistics\cr
Sample sizes for clinical trials with Normal data"\cr
'Statistics in Medicine, Vol. 23, 1921-1986 (2004)\cr
Multiplicative model, theta1=0.9, theta2=1.1111 (1/theta1), target power=90\%,\cr
power approximate via non-central t-distribution.\cr
Attention! Julious gives sample size per group.

ctSJ.VIII.20 is Table VIII from the same source\cr
column 'level of bioequivalence 20\%'\cr
Multiplicative model, theta1=0.8, theta2=1.25 (1/theta1), target power=90\%,\cr
power approximate via non-central t.\cr

ctCW.III is Table III from\cr
Chow and Wang\cr
"On Sample Size Calculation in Bioequivalence Trials"\cr
J. Pharmacokin. Biopharm. Vol. 28(2), 155-169 (2001)\cr
Additive model, theta1=-0.2, theta2=+0.2 (BE limits 0.80 - 1.20), exact.\cr

Seems the last reference is not very reliable (compare to the Table in the paper).
}
\note{Scripts for creation of these data.frame's can be found in the \\test
sub-directory of the package.\cr
Comparing the results of that scripts to the corresponding data.frames can 
be used for validation purposes.}
\author{PowerTOST}
