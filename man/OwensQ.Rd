\name{OwensQ}
\alias{OwensQ}
\title{Owen's Q-function}
\description{Calculates Owen's Q function.}
\usage{OwensQ(nu, t, delta, a, b)}
\arguments{
  \item{nu}{degree of Owen's Q}
  \item{t}{parameter t}
  \item{delta}{parameter delta}
  \item{a}{lower integration limit}
  \item{b}{upper integration limit}
}
\details{
Uses \code{integrate()} from package stats to perform the numerical evaluation 
of the definite integral in Owen's Q function.\cr
See \url{../doc/BE_power_sample_size_excerpt.pdf} in the package sub-directory 
/doc for the definition of Owen's Q and implmentation details. \cr
In case of high delta and/or high upper integration limit b where the 
implementation via R's function \code{integrate()} may fail the function 
\code{OwensQOwen()} is used. \cr\cr
The arguments to the function must be scalars. No vectors allowed.
}
\value{
Numeric value of Owen's Q-function at given input arguments.
}
\note{
This function is intended for internal use in the power calculations.\cr
But may be useful for others.
}
\section{Warning}{
Since for really large values of nu and the upper integration limit b the integrand is a function which is zero over nearly all its range, the 
\code{integrate()} function may fail (see \code{?integrate}) and 
\code{OwensQ()} then returns erroneously 0.\cr
The function now tries to return a value via non-central t-approximation in 
such cases. This approximation is up to 6 decimals correct as far as tested.\cr 
\code{OwensQ()} issues a warning if the nct-approximation is used.
}
\author{D. Labes}
\references{
Owen, D. B. (1965) \cr
"A Special Case of a Bivariate Non-central t-Distribution"\cr
Biometrika, 52, 437-446.
}
\seealso{
\code{\link{OwensQOwen}}
}
\examples{
# This function is mainly intended for internal use.
OwensQ(10,2.5,5,0,2)
#should give [1] 9.388137e-06 
OwensQ(10,-2.5,-5,0,2)
#should give [1] 0.05264363 
}

