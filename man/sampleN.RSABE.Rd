\encoding{utf-8}
\name{sampleN.RSABE}
\alias{sampleN.RSABE}
\title{
Sample size estimation for BE decision via linearized scaled ABE criterion
}
\description{
This function performs the Sample size estimation for the BE decision via 
linearized scaled ABE criterion based on simulations.
}
\usage{
sampleN.RSABE(alpha = 0.05, targetpower = 0.8, theta0, theta1, theta2, CV,
              design = c("2x3x3", "2x2x4", "2x2x3"), regulator = c("FDA", "EMA"),
              nsims = 1e+05, nstart, imax=100,
              print = TRUE, details = TRUE, setseed=TRUE)
}
\arguments{
  \item{alpha}{
Type I error probability. Per convention mostly set to 0.05.
}
  \item{targetpower}{
Power to achieve at least. Must be >0 and <1.\cr
Typical values are 0.8 or 0.9.
}
  \item{theta0}{
'True' or assumed bioequivalence ratio. \cr
Defaults to 0.90 according to the two Laszlo's if not given explicitly.
}
  \item{theta1}{
Conventional lower ABE limit to be applied in the mixed procedure if 
CVsWR <= CVswitch.\cr
Also Lower limit for the point estimate constraint.\cr
Defaults to 0.8 if not given explicitly.
}
  \item{theta2}{
Conventional upper ABE limit to be applied in the mixed procedure if 
CVsWR <= CVswitch. Also upper limit for the point estimate constraint.\cr
Defaults to 1.25 if not given explicitly.
}
  \item{CV}{
Coefficient(s) of variation as ratio.\cr
If \code{length(CV) = 1} the same CV is assumed for Test and Reference.\cr
If \code{length(CV) = 2} the CV for Test must be given in CV[1] and for 
Reference in CV[2].
}
  \item{design}{
Design of the study to be planned.\cr
2x3x3 is the partial replicate design (TRR|RTR|RRT).\cr
2x2x4 is the full replicate design with 2 sequences and 4 periods.\cr
2x2x3 is the 3-period design with sequences (TRT|RTR).\cr
Defaults to \code{design="2x3x3"}
}
  \item{regulator}{
Regulatory body settings for the scaled ABE criterion.\cr
Defaults to \code{design="FDA"}.\cr
Also the scaled ABE criterion is usually calculated with the FDA constant
r_const=log(1.25)/0.25 you can override this behavior to use the EMA setting
r_const=0.76 to avoid the discontinuity at CV=30\% and be more stringent.
}
  \item{nsims}{
Number of simulations to be performed to obtain the (empirical) power.
}
  \item{nstart}{
Set this to a start for the sample size search if a previous run failed.\cr
After reworking the start n in version 1.1-05 seldom needed.
}
  \item{imax}{
Maximum number of steps in sample size search. Defaults to 100.
}
  \item{print}{
If \code{TRUE} (default) the function prints its results.\cr 
If \code{FALSE} only the result data.frame will be returned. 
}
  \item{details}{
If set to \code{TRUE}, the default, the steps during sample size search are shown.
}
  \item{setseed}{
Simulations are dependent on the starting point of the (pseudo) random number 
generator. To avoid differences in power for different runs a \code{set.seed(123456)} 
is issued if \code{setseed=TRUE}, the default. 
}
}
\details{
The linearized scaled ABE criterion is calculated according to the SAS code
given in the FDA progesterone guidance.\cr\cr
The simulations are done via the distributional properties of the statistical
quantities necessary for deciding BE based on scaled ABE. For more details see
a document "Implementation_scaledABE_simsVx.yy.pdf" in the doc subdirectory of 
the package.
}
\value{
Returns a data.frame with the input and sample size results.\cr
The "Sample size" column contains the total sample size.\cr
The "nlast" column contains the last n value. May be useful for restarting.
}
\references{
FDA "Draft Guidance on Progesterone"\cr
Recommended Apr 2010; Revised Feb 2011\cr
http://www.fda.gov/downloads/Drugs/GuidanceComplianceRegulatoryInformation/
Guidances/UCM209294.pdf

\enc{Lászlo}{Laszlo} \enc{Tóthfalusi}{Tothfalusi} and 
\enc{Lászlo}{Laszlo} \enc{Endrényi}{Endrenyi}\cr
"Sample Sizes for Designing Bioequivalence Studies for Highly Variable Drugs"\cr
J. Pharm. Pharmaceut. Sci. (www.cspsCanada.org) 15(1) 73 - 84, 2011\cr
\href{http://ejournals.library.ualberta.ca/index.php/JPPS/article/download/11612/9489}{free download}

\enc{Tóthfalusi}{Tothfalusi} L., \enc{Endrényi}{Endrenyi} L. and A. Garcia Arieta\cr
"Evaluation of Bioequivalence for Highly Variable Drugs with Scaled Average 
Bioequivalence" \cr
Clin. Pharmacokin. 48/11, 725-743 (2009)
}
\author{
D. Labes
}
\note{
The sample size estimation is done only for balanced designs since the 
break down of the total subject number in case of unbalanced sequence groups
is not unique. Moreover the formulas used are only for balanced designs.\cr
The minimum sample size is n=6, even if the power is higher than the intended 
targetpower.
}
\section{Warning }{
The sample size estimation for theta0 >1.2 and <0.85 may be very time consuming 
and will eventually also fail since the start values chosen are not really 
reasonable in that ranges. This is especially true in the range about CV = 0.3 
and regulatory constant according to FDA.\cr
If you really need sample sizes in that range be prepared to restart the sample 
size estimation via the argument nstart.\cr
Since the dependence of power from n is very flat in the mentioned region you 
may also consider to adapt the number of simulations not to tap in the simulation 
error trap.
}
\seealso{
\code{\link{power.RSABE}}, \code{\link{power.scABEL}}
}
\examples{
# using all the defaults:
# design=2x3x3 (partial replicate design), theta0=0.90, 
# ABE limits, PE constraint 0.8 - 1.25
# targetpower=80\%, alpha=0.05, 1E5 sims
sampleN.RSABE(CV=0.3)
# should result in a sample size n=45, power=0.80344
}
