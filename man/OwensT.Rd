\encoding{utf-8}
\name{OwensT}
\alias{OwensT}
\title{
Owen's T-function
}
\description{
Calculates the definite integral from \code{0} to \code{a} of\cr
\verb{  exp(-0.5*h^2*(1+x^2))/(1+x^2)/(2*pi)}.
}
\usage{
OwensT(h, a)
}
\arguments{
  \item{h}{parameter h}
  \item{a}{upper limit of integration}
}
\details{
The function is an R port of FORTRAN code given in the references and MATLAB 
code given by John Burkardt under the GNU LGPL license.\cr\cr
The arguments of \code{OwensT()} have to be scalars because the implementation 
doesn’t vectorize.
}
\value{
Numerical value of the definite integral.
}
\author{
MATLAB code by J. Burkardt, R port by D. Labes
}
\note{
This function is only needed as auxiliary in \code{\link{OwensQOwen}}.\cr
But may be useful for others.
}
\references{
Goedhart PW, Jansen MJW. \emph{Remark AS R89: A Remark on Algorithm AS 76: An Integral Useful in Calculating Central t and Bivariate Normal Probabilities.} J Royal Stat Soc C. 1992;41(2):496--7. \doi{10.2307/2347586}

Boys R. \emph{Algorithm AS R80: A Remark on Algorithm AS 76: An Integral Useful in Calculating Noncentral t and Bivariate Normal Probabilities.} J Royal Stat Soc C. 1989;38(3):580--2. \doi{10.2307/2347755}

Thomas GE. \emph{Remark ASR 65: A Remark on Algorithm AS76: An Integral Useful in Calculating Non-Central t and Bivariate Normal Probabilities.} J Royal Stat Soc C. 1986;35(3):310--2. \doi{10.2307/2348031}

Chou Y-M. \emph{Remark AS R55: A Remark on Algorithm AS 76: An Integral Useful in Calculating Noncentral T and Bivariate Normal Probabilities.} J Royal Stat Soc C. 1985;34(1):100--1. \doi{10.2307/2347894}

Thomas GE. \emph{Remark AS R30: A Remark on Algorithm AS 76: An Integral Useful in Calculating Non-Central t and Bivariate Normal Probabilities.} J Royal Stat Soc C. 1979;28(1):113. \doi{10.2307/2346833}

Young JC, Minder C. \emph{Algorithm AS 76: An Integral Useful in Calculating Non-Central t and Bivariate Normal Probabilities.} J Royal Stat Soc C. 1974;23(3):455--7. \doi{10.2307/2347148}

Burkardt J. \emph{ASA076. Owen's T Function.} \url{https://people.math.sc.edu/Burkardt/f_src/asa076/asa076.html}

Owen DB. \emph{Tables for Computing Bivariate Normal Probabilities.} Ann Math Stat. 1956;27(4):1075--90. \doi{10.1214/aoms/1177728074}
}
\seealso{
\code{\link{OwensQOwen}, \link{OwensQ}}
}
\examples{
OwensT(2.5, 0.75)
# should give [1]  0.002986697
# value from Owen's tables is 0.002987
OwensT(2.5, -0.75)
# should give [1] -0.002986697
}
