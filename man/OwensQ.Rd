\name{OwensQ}
\alias{OwensQ}
\title{Owen's Q-function}
\description{Calculates Owen's Q function.}
\usage{OwensQ(nu, t, delta, a=0, b)}
\arguments{
  \item{nu}{degree of Owen's Q}
  \item{t}{parameter t}
  \item{delta}{parameter delta}
  \item{a}{lower integration limit, only a=0 implemented}
  \item{b}{upper integration limit}
}
\details{
Uses now the relationship to non-central t-distribution (see Chou YM)\cr\cr
   \eqn{OwensQ = pt(t, df=nu, ncp=delta) - Integal_b_Inf(Q_integrand)}\cr\cr
The definite integral is numerically evaluated using \code{integrate()} from 
package \code{stats} after a variables transformation resulting in the 
integration range from 0 to 1 instead of the semi-infinite original range. 
This may result in higher precision and better numerical stability.\cr\cr
The arguments to the function must be scalars. No vectors allowed.
}
\value{
Numeric value of Owen's Q-function at given input arguments.
}
\note{
This function is intended for internal use in the power calculations.\cr
But may be useful for others.
}
\author{D. Labes}
\references{
Owen, DB (1965) \cr
"A Special Case of a Bivariate Non-central t-Distribution"\cr
Biometrika, 52, 437-446.\cr
\doi{10.2307/2333696}

Chou YM (1992)\cr
"A bivariate noncentral T-distibution with applications"\cr
Communications in Statistics - Theory and Methods, 21:12, 3427-3462\cr 
\doi{10.1080/03610929208830988}

}
\seealso{
\code{\link{OwensQOwen}}
}
\examples{
# This function is mainly intended for internal use.
OwensQ(10, 2.5, 5, 0, 2)
#should give [1] 9.388137e-06 
OwensQ(10, -2.5, -5, 0, 2)
#should give [1] 0.05264363 
}
