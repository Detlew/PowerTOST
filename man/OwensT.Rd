\name{OwensT}
\alias{OwensT}
\title{
Owen's T-function
}
\description{
Calculates the definite integral from 0 to a of exp(-0.5*h^2*(1+x^2))/(1+x^2)/(2*pi).
}
\usage{
OwensT(h, a)
}
\arguments{
  \item{h}{parameter h}
  \item{a}{upper limit of integration}
}
\details{
The function is simply implemented via stats function \code{integrate()}.
}
\value{
Numeric value of the definite integral.
}
\author{
D. Labes
}
\note{
This function is only needed in \code{OwensQOwen()}.\cr
But may be useful for others.
}
\seealso{
\code{\link{OwensQOwen}, \link{OwensQ}}
}
\examples{
OwensT(2.5,0.75)
# should give [1] 0.002986697
}

