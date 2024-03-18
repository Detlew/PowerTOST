\encoding{utf-8}
\name{defunct}
\alias{defunct}
\title{Removed functions in PowerTOST}
\description{
Details about removed functions in PowerTOST.
}
\details{
\code{type1error.2TOST} was designed to calculate the type I error (\emph{TIE}) rate of two simultaneous TOST procedures with some specified correlation of the parameters. It suffered from poor precision to obtain the \emph{TIE} via simulations.\cr
Due to the intersection-union principle the \emph{TIE} is always upper bounded to \ifelse{latex}{\out{$\alpha$}}{\ifelse{html}{\out{<em>&alpha;</em>}}{alpha}} by theory and hence, \code{type1error.2TOST} was removed in version 1.4-7.

\code{power.scABEL2} and \code{sampleN.scABEL2} were deprecated in version 1.4-1. In \code{power.scABEL} and \code{sampleN.scABEL} the regulator component \code{"est_method"} is used for switching between simulations based on the EMA’s ANOVA evaluation or ISC evaluation, respectively.\cr
\code{power.scABEL2} and \code{sampleN.scABEL2} were removed in version 1.4-3.

\code{power2.TOST} was deprecated in version 1.2-6 since \code{power.TOST} was modified in order to handle unbalanced sequences. \code{power2.TOST} was removed in version in version 1.2-7.

The functions \code{power.NTIDFDA}, \code{sampleN.NTIDFDA} and \code{pa.NTIDFDA} were deprecated in 
version 1.5-4 and were removed in version 1.5-6.

}

\references{
Berger RL, Hsu JC. \emph{Bioequivalence Trials, Intersection-Union Tests and Equivalence Confidence Sets.} Stat Sci. 1996;11(4):283--302. \href{https://www.jstor.org/stable/2246021}{JSTOR:2246021}
}

\seealso{
\code{\link{power.TOST}}, \code{\link{power.scABEL}}, \code{\link{sampleN.scABEL}}, \code{\link{power.NTID}}, \code{\link{sampleN.NTID}}, and \code{\link{pa.NTID}} for the new functions.
}


