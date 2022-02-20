\encoding{utf-8}
\name{defunct}
\alias{defunct}
\title{Removed functions in PowerTOST}
\description{
Details about removed functions in PowerTOST.
}
\details{
\code{power.scABEL2} and \code{sampleN.scABEL2} were deprecated in version 1.4.1. In \code{power.scABEL} and \code{sampleN.scABEL} the regulator component \code{"est_method"} is used for switching between simulations based on the EMA’s ANOVA evaluation or ISC evaluation, respectively.\cr
\code{power.scABEL2} and \code{sampleN.scABEL2} were removed in version 1.4.3.

\code{power2.TOST} was deprecated in version 1.2.6 since \code{power.TOST} was modified in order to handle unbalanced sequences. \code{power2.TOST} was removed in version in version 1.2.7.
}

\seealso{
\code{\link{power.TOST}}, \code{\link{power.scABEL}}, \code{\link{sampleN.scABEL}} for the new functions.
}
