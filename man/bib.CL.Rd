\name{bib.CL}
\alias{bib.CL}
\title{
Design matrices of period balanced incomplete block designs
}
\description{
This function returns the 'design' matrix of incomplete block designs described
in Chow & Liu's book. The design matrices where recoded 1=R, 2=T1, 3=T2 ...

}
\usage{
bib.CL(trt, p)
}
\arguments{
  \item{trt}{
Number of treatments (3 to 5). 
}
  \item{p}{
Number of periods (2 to trt-1).
}
}
\value{
Matrix containing the sequences in rows and periods in columns.\cr 
The entry (i,j) of the matrix corresponds to the treatment or dose (index) a 
subject within i-th sequence gets in the j-th period.
}
\references{
Chow, Liu\cr
"Design and Analysis of Bioavailability and Bioequivalence Studies"\cr
Third edition, Chapman and Hall/CRC 2008\cr
Chapter 2.6
}
\author{
D. Labes
}
\examples{
# 4 treatments/doses, 3 periods
bib.CL(4,3)
# gives 4 sequences
# to see this in Chow & Liu's coding
tmt <- c("R", "T1", "T2", "T3")
matrix(tmt[bib.CL(4,3)], ncol=3)
}
