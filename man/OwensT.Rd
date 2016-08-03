\name{OwensT}
\alias{OwensT}
\title{
Owen's T-function
}
\description{
Calculates the definite integral from \code{0} to \code{a} of 
\code{exp(-0.5*h^2*(1+x^2))/(1+x^2)/(2*pi)}.
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
code given on \url{http://people.sc.fsu.edu/~jburkardt/m_src/asa076/asa076.html}
by John Burkardt under the GNU LGPL license.\cr\cr
The arguments of \code{OwensT()} have to be scalars because the implementation 
doesn't vectorize.
}
\value{
Numerical value of the definite integral.
}
\author{
MATLAB code by John Burkardt\cr
R port by D. Labes
}
\note{
This function is only needed as auxiliary in \code{OwensQOwen()}.\cr
But may be useful for others.
}
\references{
Goedhart PW, Jansen MJW\cr
Remark AS R89: A Remark on Algorithm AS76: An Integral Useful in Calculating 
Noncentral T and Bivariate Normal Probabilities\cr
J. of the Royal Statistical Society. Series C (Applied Statistics), Vol. 41 (2), 1992, pages 496-497\cr
\href{http://dx.doi.org/10.2307/2347586}{doi 10.2307/2347586}

Boys R \cr
Remark AS R80: A Remark on Algorithm AS76: An Integral Useful in Calculating 
Noncentral T and Bivariate Normal Probabilities \cr
J. of the Royal Statistical Society. Series C (Applied Statistics), Vol. 38 (3), 1989, pp. 580-582\cr
\href{http://dx.doi.org/10.2307/2347755}{doi 10.2307/2347755}

Thomas GE\cr
Remark ASR 65: A Remark on Algorithm AS76: An Integral Useful in Calculating 
NonCentral t and Bivariate Normal Probabilities\cr
J. of the Royal Statistical Society. Series C (Applied Statistics), Vol. 35 (3), 1986, pp. 310-312\cr
\href{http://dx.doi.org/10.2307/2348031}{doi 10.2307/2348031}

Chou Y-M \cr
Remark AS R55: A Remark on Algorithm AS76: An Integral Useful in Calculating 
Noncentral T and Bivariate Normal Probabilities\cr
J. of the Royal Statistical Society. Series C (Applied Statistics), Vol. 34 (1), 1985, pp. 100-101\cr
\href{http://dx.doi.org/10.2307/2347894}{doi 10.2307/2347894}

Porter MA, Winstanley DJ\cr
Remark AS R30: A Remark on Algorithm AS76: An Integral Useful in Calculating 
Noncentral T and Bivariate Normal Probabilities\cr
J. of the Royal Statistical Society. Series C (Applied Statistics), Vol. 28 (1), 1979, p. 113\cr
\href{http://dx.doi.org/10.2307/2346833}{doi 10.2307/2346833}

Young JC, Minder C\cr
Algorithm AS 76: An Algorithm Useful in Calculating Non-Central T and Bivariate 
Normal Distributions\cr
J. of the Royal Statistical Society. Series C (Applied Statistics), Vol. 23 (3), 1974, pp. 455-457\cr
\href{http://dx.doi.org/10.2307/2347148}{doi 10.2307/2347148}

Owen DB\cr
Tables for Computing Bivariate Normal Probabilities\cr
Annals of Mathematical Statistics, Vol. 27 (4), December 1956, pp 1075-1090\cr
\href{http://dx.doi.org/10.1214/aoms/1177728074}{doi 10.1214/aoms/1177728074}
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
