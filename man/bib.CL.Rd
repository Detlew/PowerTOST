\encoding{utf-8}
\name{bib.CL}
\alias{bib.CL}
\title{
Design matrices of period balanced incomplete block designs
}
\description{
  This function returns the \sQuote{design} matrix of incomplete
  block designs described by Chow & Liu. The design matrices were
  recoded \verb{1=R}, \verb{2=T1}, \verb{3=T2}, \ldots
}

\usage{
bib.CL(trt, p)
}

\arguments{
  \item{trt}{
    Number of treatments (\verb{3} to \verb{5}). 
  }
  \item{p}{
    Number of periods (\verb{2} to \verb{trt-1}).
  }
}

\value{
  Matrix containing the sequences in rows and periods in columns.\cr 
  The entry \verb{(i, j)} of the matrix corresponds to the treatment or
  dose (index) a subject within \emph{i}-th sequence gets in the
  \emph{j}-th period.
}

\author{
D. Labes
}

\references{
Chow SC, Liu JP. \emph{Design and Analysis of Bioavailability and Bioequivalence Studies.} Boca Raton: CRC Press; 3\ifelse{html}{\out{<sup>rd</sup>}}{\eqn{^{rd}}} edition 2009. Chapter 2.6.
}

\examples{
# 4 treatments/doses, 3 periods
bib.CL(4, 3)
# gives 4 sequences
# to see this in Chow & Liu's coding
tmt <- c("R", "T1", "T2", "T3")
matrix(tmt[bib.CL(4, 3)], ncol=3)
}
