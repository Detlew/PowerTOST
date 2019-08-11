\encoding{utf-8}
\name{OwensQOwen}
\alias{OwensQOwen}
\title{
Owen's Q-function via repeated integration by parts
}
\description{
This is an implementation of the algorithm given by Owen via repeated integration by parts. 
}
\usage{
OwensQOwen(nu, t, delta, a=0, b)
}
\arguments{
  \item{nu}{degree of Owen’s Q}
  \item{t}{parameter t}
  \item{delta}{parameter delta}
  \item{a}{lower integration limit.\cr
           Only \code{a=0} implemented, other values give an error.}
  \item{b}{upper integration limit}
}
\value{
Numeric value of Owen’s Q function.
}
\note{
The argument \code{a=0} could be dropped but is retained for sake of completeness. 
}
\references{
Owen DB. \emph{A special case of a bivariate non-central t-distribution.} Biometrika. 1965;52(3/4):437--46. \doi{10.2307/2333696}
}
\author{
D. Labes
}
\note{
This function is mainly for comparative / validation purposes.\cr
The function requireds \code{\link{OwensT}} function.
}

\seealso{
\code{\link{OwensQ}}, \code{\link{OwensT}}
}
\examples{
# comparison of the results of both implementations
# both should give [1] 0.0731726
OwensQ(2, 2.92, 4.2135, 0, 2.0407)
OwensQOwen(2, 2.92, 4.2135, 0, 2.0407)
}

