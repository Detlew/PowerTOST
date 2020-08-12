\encoding{utf-8}
\name{sampleN.dp}
\alias{sampleN.dp}
\title{
Sample size estimation of dose-proportionality studies evaluated via the power model
}
\description{
Performes a sample size estimation for dose-proportionality studies using the 
power model for cossover (Latin square), parallel group designs or incomplete
block designs via a confidence interval equivalence criterion.
}
\usage{
sampleN.dp(alpha = 0.05, CV, doses, targetpower = 0.8, beta0, theta1 = 0.8, 
           theta2 = 1/theta1, design = c("crossover", "parallel", "IBD"), 
           dm=NULL, CVb, print = TRUE, details = FALSE, imax = 100)
}
\arguments{
  \item{alpha}{
Type 1 error. Usually set to 0.05.
}
  \item{CV}{
Coefficient of variation. Is intra-subject CV for \code{design="crossover"} and
CV of total variability in case of \code{design="parallel"}
}
  \item{doses}{
Vector of dose values under study. At least two doses have to be given.
}
  \item{targetpower}{
Power to achieve at least. Must be >0 and <1.\cr
Typical values are 0.8 or 0.9.
}
  \item{beta0}{
\sQuote{True} or assumed slope of the power model. If missing defaults to \code{1+log(0.95)/log(rd)}
where \code{rd} is the ratio is the ratio of the highest to the lowest dose.\cr 
Has to be within slope acceptance range according to
\code{1+log(theta1)/log(rd)} and \code{1+log(theta2)/log(rd)}. Otherwise, the 
function issues an error.
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
Can be obtained f.i. via functions of package \code{crossdes}. See examples.\cr
Function \code{\link{bib.CL}} returns some IBDs described by Chow & Liu.
}
  \item{CVb}{
Coefficient of variation of the between-subject variability.\cr
Only necessary if \code{design="IBD"}. Will be set to 2*CV if missing.
This is only a crude rule of thumb. Better obtain an estimate of CVb from a previous crossover study.\cr\cr
Set \code{CVb=0} if all-effects-fixed model shall be used. This model gives
lower sample sizes than the mixed model with random subject effects 
(random intercept).  
}
  \item{print}{
If \code{TRUE} (default) the function prints its results.\cr 
If set to \code{FALSE} only the data.frame with the results will be returned. 
}
  \item{details}{
If \code{details=TRUE} the steps during sample size search will be shown. Defaults to \code{FALSE}.
}
  \item{imax}{
Maximum number of steps in sample size search. \cr 
Defaults to 100. Adaption only in rare cases needed, if any.
}
}
\details{
The sample size is estimated via iterative evaluation of \code{power.dp()}.\cr
Start value for the sample size search is taken from a large sample approximation.\cr
The sample size is bound to number of dose or sequence groups as minimum.\cr
Balanced designs are used although this is not absolutely necessary.\cr\cr
The estimated sample size gives always the \emph{total} number of subjects (not subject/sequence in crossovers or subjects/group in parallel designs -- like in some other software packages).
}
\value{
A data.frame with the input and results will be returned.\cr
The \code{Sample size} column contains the total sample size.
}
\references{
Chow SC, Liu JP. \emph{Design and Analysis of Bioavailability and Bioequivalence Studies.} Boca Raton: CRC Press; 3\ifelse{html}{\out{<sup>rd</sup>}}{\eqn{^{rd}}} edition 2009.

Patterson S, Jones B. \emph{Bioequivalence and Statistics in Clinical Pharmacology.} Boca Raton: Chapman & Hall/CRC: 2006. p. 239.\cr
(contains presumably a bug)

Sethuraman VS, Leonov S, Squassante L, Mitchell TR, Hale MD. \emph{Sample size calculation for the Power Model for dose proportionality studies.} Pharm Stat. 2007;6(1):35--41. \doi{10.1002/pst.241}

Hummel J, McKendrick S, Brindley C, French R. \emph{Exploratory assessment of dose proportionality: review of current approaches and proposal for a practical criterion.} Pharm. Stat. 2009;8(1):38--49. \doi{10.1002/pst.326}
}
\author{
D. Labes
}
\section{Warning }{
This function is \sQuote{experimental} only, since it is not thorougly tested yet. 
Especially for \code{design="IBD"} reliable test cases are missing.
}
\seealso{
\code{\link{power.dp}}, \code{\link{bib.CL}}
}
\examples{
# using all the defaults, i.e. crossover design, alpha=0.05
# theta1=0.8, theta2=1.25 but true slope slightly off 1
sampleN.dp(CV = 0.2, doses = c(1, 2, 8), beta0 = 1.02)
# should give n=18, power=0.854528

\donttest{
# incomplete block design with 5 doses, 3 periods
# from library(crossdes)
doses  <- c(5, 25, 50, 100, 200)
CVb    <- mse2CV(0.8)
levels <- length(doses)
per    <- 3
block  <- levels*(levels-1)/(per-1)
# IBD based on balanced minimal repeated measurements design
# gives n=30 and 10 sequences
ibd    <- crossdes::balmin.RMD(levels, block, per)
sampleN.dp(CV = 0.2, doses = doses, beta0 = 1, design = "IBD", dm = ibd, 
           CVb = CVb, targetpower=0.9)}
}
