\encoding{utf-8}
\name{known.designs}
\alias{known.designs}
\title{Show the 'known' designs}
\description{
Returns the known study designs for which power and sample size can be calculated within this package.
}
\usage{
known.designs()
}
\details{
This function is for informal purposes and used internally for obtaining characteristics of the designs used in calculation formulas.
}
\value{
Returns a data.frame with \cr
  \item{\code{no}}{number of the design}
  \item{\code{design}}{character string for identifying the design}
  \item{\code{df}}{degrees of freedom of the design}
  \item{\code{df2}}{\sQuote{robust} degrees of freedom of the design}
  \item{\code{steps}}{step width in the iterative sample size estimation}
  \item{\code{bk}}{so-called design constant in terms of total n}
  \item{\code{bkni}}{design constant in terms of number of subjects in (sequence) groups}

The design character string has to be used in the functions calls for power and sample size.
}
\note{
The design string for higher order crossover designs is named as:\cr
\code{treatments x sequences x periods} in case of replicate designs and\cr
\code{treatments x periods} in case of crossover designs for more then 2 treatments 
with number of sequences equal number of treatments.\cr\cr
The df for the replicate crossover designs are those without carry-over in the model.\cr
Chen \emph{et al.} used models with carry-over, \emph{i.e.}, one df lower than here.\cr\cr
The design constant bk in case of design 2x2x4 is here bk=1.\cr
Chen \emph{et al.} used bk=1.1 due to carry-over in the model.\cr\cr
n is the \bold{total} number of subjects for all designs implemented.\cr
df2 = degrees of freedom for the so-called \sQuote{robust} analysis (aka \enc{Senn’s}{Senn's} basic estimator).\cr
These degrees of freedom are often also more appropriate in case of evaluation via a \sQuote{true} 
mixed model (\emph{e.g.} the FDA\enc{’}{'} for replicate designs).\cr\cr
The design \code{2x2x2r} is the 2-treatment-2-sequence-2-period design with 
2 repeated targets determined in each period (sequences TT|RR or RR|TT) described 
by Liu. Implemented are the characteristics of this design for the evaluation 
via assuming no S×F interaction and equal variabilities of Test and Reference.
}

\references{
Chen KW, Chow SC, Liu G. \emph{A Note on Sample Size Determination for Bioequivalence Studies with 
Higher-order Crossover Designs.} J Pharmacokin Biopharm. 1997;25(6):753--65.

Senn S. \emph{Cross-over Trials in Clinical Research.} Chichester: John Wiley & Sons; 2\ifelse{html}{\out{<sup>nd</sup>}}{\eqn{^{nd}}} edition 2002.

U.S. Department of Health and Human Services, Food and Drug Administration, Center for Drug Evaluation and Research (CDER). \emph{Guidance for Industry. Statistical Approaches to Establishing Bioequivalence.} January 2001. \href{https://www.fda.gov/media/70958/download}{download}

Liu J-p. \emph{Use of the Repeated Crossover design in Assessing Bioequivalence.} Stat Med. 1995;14(9-10):1067--78.
}
\author{
D. Labes
}
\examples{
known.designs()
}
