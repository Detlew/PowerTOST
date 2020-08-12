\name{power.dp}
\alias{power.dp}
\title{
Power of dose-proportionality studies evaluated via Power model
}
\description{
Calculates the power of dose-proportionality studies using the power model for 
crossover (Latin square) or parallel group designs via a confidence interval
equivalence criterion.
}
\usage{
power.dp(alpha = 0.05, CV, doses, n, beta0, theta1 = 0.8, theta2 = 1/theta1, 
         design = c("crossover", "parallel", "IBD"), dm = NULL, CVb)
}
\arguments{
  \item{alpha}{
Type 1 error. Commonly set to 0.05.
}
  \item{CV}{
Coefficient of variation for intra-subject variability if \code{design="crossover"} 
or CV of total variability in case of \code{design="parallel"}.
}
  \item{doses}{
Vector of dose levels. At least two doses have to be given.
}
  \item{n}{
Number of subjects. Is total number if given as scalar, else number of subjects
in the (sequence) groups. In the latter case the length of n vector has to be 
the same as length of vector doses.\cr
n has to be >2.
}
  \item{beta0}{
\sQuote{True} slope of power model. If missing defaults to \code{1+log(0.95)/log(rd)}
where \code{rd} is the ratio of highest to lowest dose.
}
  \item{theta1}{
Lower acceptance limit for the ratio of dose normalized means (Rdmn).\cr
Transformes into slope acceptance range as described under item \code{beta0}.
}
  \item{theta2}{
Upper acceptance limit for the ratio of dose normalized means (Rdmn).
}
  \item{design}{
Crossover design (default), parallel group design or incomplete block design (IBD).\cr
Crossover design means Latin square design with number of doses as dimension.
}
  \item{dm}{
'Design matrix' of the incomplete block design (IBD) if \code{design="IBD"}.\cr
This matrix contains the sequences in rows and periods in columns. 
The entry (\emph{i}, \emph{j}) of the design matrix corresponds to the dose (index) a subject 
with \emph{i}-th sequence gets in the \emph{j}-th period.
Can be obtained f.i. via functions of package \code{crossdes} or via function 
\code{bib.CL()}.
}
  \item{CVb}{
Coefficient of variation of the between-subject variability.\cr
Only necessary if \code{design="IBD"}. Will be set to 2*CV if missing.
This is only a crude rule of thumb. Better obtain an estimate of CVb from a 
previous crossover study.\cr\cr
Set \code{CVb=0} if an all-effects-fixed model shall be used. This model gives
higher power than the random subject effects model.
}
}
\details{
The power calculations are based on TOST for testing equivalence of the slope 
of the power model with alternativ hypothesis slope = 1.\cr
Power is calculated via non-central t-approximation only.\cr
The calculations are based on mixed effects model (random intercept aka
random subject effect). For \code{design="cossover"} or \code{design="parallel"}
the results coincide with all-effects-fixed model.
}
\value{
Value of power according to the input arguments.
}
\references{
Patterson S, Jones B. \emph{Bioequivalence and Statistics in Clinical Pharmacology.} Boca Raton: Chapman & Hall/CRC: 2006. p. 239.\cr
(contains presumably a bug)

Sethuraman VS, Leonov S, Squassante L, Mitchell TR, Hale MD. \emph{Sample size calculation for the Power Model for dose proportionality studies.} Pharm Stat. 2007;6(1):35--41. \doi{10.1002/pst.241}

Hummel J, McKendrick S, Brindley C, French R. \emph{Exploratory assessment of dose proportionality: review of current approaches and proposal for a practical criterion.} Pharm. Stat. 2009;8(1):38--49. \doi{10.1002/pst.326}
}
\author{
D. Labes
}
\section{Warning }{
This function is \sQuote{experimental} only since it is not thorougly tested yet. 
Especially for \code{design="IBD"} reliable test cases are missing.
}

\seealso{
\code{\link{sampleN.dp}}, \code{\link{bib.CL}}
}
\examples{
# using all the defaults, i.e. latin square crossover design, alpha=0.05, 
# beta0=1+log(0.95)/log(rd), theta1=0.8, theta2=1.25
power.dp(CV = 0.2, doses = c(1,2,8), n = 15)
#
# period balanced IBD with 3 doses, 2 periods and 3 sequences,
ibd <- matrix(c(1, 2, 3, 2, 3, 1), nrow = 3, ncol = 2)
power.dp(CV = 0.2, doses = c(1,2,8), n = 12, design = "IBD", dm = ibd)
# considerably lower than 3x3 Latin square
}
