\name{scABEL}
\alias{scABEL}
\title{
Scaled (widened) BE Acceptance Limits
}
\description{
The (widened) scaled BE acceptance limits are calculated according to the 
regulatory settings of EMA, HC, FDA or via user defined regulatory settings.
}
\usage{
scABEL(CV, regulator)
}
\arguments{
  \item{CV}{
Coefficient of variation (of the Reference) as ratio.
}
  \item{regulator}{
Regulatory body settings for the widening of the BE acceptance limits.\cr
May be given as character from the choices "EMA", "HC", "FDA" or as an object of
class 'regSet' (see \code{\link{reg_const}}).\cr
Defaults to \code{regulator="EMA"} if missing.\cr\cr
The former \code{regulator="ANVISA"} is no longer allowed. Use "EMA" since ANVISA 
now recommends the use of EMA regulatory settings.\cr\cr
The former \code{regulator="USER"} is no longer accepted but can be handled now
via function \code{reg_const()} to define an object with class 'regSet'.
}
}
\value{
Returns a vector of lenghth 2 if one CV is given or a matrix if CV is given as vector 
with named components "lower" and "upper" of the scaled acceptance limits.
}
\details{
The widened BE acceptance limits will be calculated by the formula\cr
   \code{[lBEL,uBEL] =exp(-+ r_const*sWR)}\cr
with r_const the regulatory constant and sWR the standard deviation of the within
subjects variability of the Reference.\cr 
r_const=0.76 (~log(1.25)/0.29356) is used in case of \code{regulator="EMA"} or 
\code{regulator="HC"} and in case of \code{regulator="FDA"} 
r_const=0.89257...(=log(1.25)/0.25).\cr
If the CVwR of the Reference is < CVswitch=0.3 the conventional ABE limits 
apply (mixed procedure).\cr\cr 
In case of \code{regulator="EMA"} a cap is placed on the widened limits if 
CVwr>0.5, i.e. the widened limits are held at value calculated for CVwR=0.5.
In case of \code{regulator="HC"} the capping is done such that the acceptance
limits are 0.6666 ... 1.5 at maximum, i.e. CVcap=0.57382. 
Literally it was given rounded to 3 significant digits as 57.4\%.\cr
See \url{http://www.hc-sc.gc.ca/dhp-mps/prodpharma/activit/announce-annonce/notice-avis-be-hvdp-nb-pphv-eng.php} \cr\cr

}
\note{
The scaled acceptance limits are not used directly in the BE evaluation 
for highly variable drugs recommended by the FDA. They are included here
for comparative purposes. Moreover there are controversies where to locate the 
'implied' BE acceptance limits. See reference below.
}
\references{
Davit et al.\cr
\emph{Implementation of a Reference-Scaled Average Bioequivalence Approach for 
Highly Variable Generic Drug Products by the US Food and Drug Administration}\cr
AAPS J. 2012;14(4):915--24. \doi{10.1208/s12248-012-9406-x}
}
\author{
D. Labes
}
\seealso{
\code{\link{power.scABEL}, \link{sampleN.scABEL}, \link{reg_const}}
}
\examples{
scABEL(CV=0.3, regulator="EMA")
# should give the usual BE limits:
# lower  upper
#  0.80   1.25
scABEL(CV=0.4, regulator="EMA")
# should give the widened limits:
#   lower     upper
# 0.746177 1.340165
#
# define old ANVISA settings via reg_const()
rc <- reg_const("USER", r_const=0.76, CVswitch=0.4, CVcap=0.5)
rc$name <- "ANVISAold"
scABEL(CV=0.4, regulator=rc)
# should give the not widened limits:
# lower upper 
#  0.80  1.25 
}
