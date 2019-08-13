\name{CV2se+se2CV+CV2mse+mse2CV}
\alias{CV2se}
\alias{se2CV}
\alias{CV2mse}
\alias{mse2CV}
\title{Helper functions}
\description{
  Calculates the standard error or the mean squared error from a 
  given CV and vice versa for log-normal data.}
\usage{
CV2se(CV)
se2CV(se)
CV2mse(CV)
mse2CV(mse)
}
\arguments{
  \item{CV}{coefficient of variatio as ratio (not percent)}
  \item{se}{standard error}
  \item{mse}{mean squared error (aka residual variance)}
}

\value{
  Returns\cr
  \verb{  se  = sqrt(log(CV^2+1))}\cr
  \verb{  CV  = sqrt(exp(se^2)-1)}\cr
  \verb{  mse = log(CV^2+1)}\cr
  \verb{  CV  = sqrt(exp(mse)-1)}
}

\note{
These functions were originally intended for internal use only but may be useful for others.
}

\author{
  D. Labes
}

\examples{
# these functions are one liners:
CV2se <- function(CV) return(sqrt(log(1.0 + CV^2)))
se2CV <- function(se) return(sqrt(exp(se^2)-1))

CV2se(0.3)
# should give:
# [1] 0.2935604 

se2CV(0.2935604)
# [1] 0.3 
}

