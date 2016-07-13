\encoding{utf-8}
\name{power.scABEL}
\alias{power.scABEL}
\alias{power.scABEL2}
\title{
(Empirical) Power of BE decision via scaled (widened) BE acceptance limits
}
\description{
These functions performs the power calculation of the BE decision via 
scaled (widened) BE acceptance limits by simulations.
}
\usage{
power.scABEL(alpha = 0.05, theta1, theta2, theta0, CV, n, 
             design = c("2x3x3", "2x2x4", "2x2x3"), regulator, 
             nsims = 1e+05, details = FALSE, setseed = TRUE)

power.scABEL2(alpha = 0.05, theta1, theta2, theta0, CV, n, 
              design = c("2x3x3", "2x2x4", "2x2x3"), regulator, 
              nsims = 1e+05, details = FALSE, setseed = TRUE)
}
\arguments{
  \item{alpha}{
Type I error probability, significance level. Conventionally mostly set to 0.05.
}
  \item{theta1}{
Conventional lower ABE limit to be applied in the mixed procedure if 
CVsWR <= CVswitch. Also lower limit for the point estimate constraint.\cr
Defaults to 0.8 if not given explicitly.
}
  \item{theta2}{
Conventional upper ABE limit to be applied in the mixed procedure if 
CVsWR <= CVswitch. Also upper limit for the point estimate constraint.\cr
Defaults to 1.25 if not given explicitly.
}
  \item{theta0}{
'True' or assumed bioequivalence ratio. \cr
Defaults to 0.90 according to the two Laszlo's if not given explicitly.
}
  \item{CV}{
Coefficient(s) of variation as ratio.\cr
If \code{length(CV) = 1} the same CV is assumed for Test and Reference.\cr
If \code{length(CV) = 2} the CV for Test must be given in CV[1] and for 
Reference in CV[2].
}
  \item{n}{
Number of subjects under study.\cr
May be given as vector. In that case it is assumed that n contains the number
of subjects in the sequence groups.\cr
If n is given as single number (total sample size) and this number is not 
divisible by the number of sequences of the design an unbalanced design is 
assumed. A corresponding message is thrown showing the numbers of subjects 
in sequence groups.\cr
Attention! In case of the 2x2x3 (TRT|RTR) design the order of n's is important 
if given as vector. n[1] is for sequence group 'TRT' and n[2] is for 
sequence group 'RTR'. 
}
  \item{design}{
Design of the study to be planned.\cr
2x3x3 is the partial replicate design (TRR/RTR/RRT).\cr
2x2x4 is the full replicate design with 2 sequences and 4 periods.\cr
2x2x3 is the 3-period design with sequences TRT|RTR.\cr
Defaults to \code{design="2x3x3"}.
}
  \item{regulator}{
Regulatory settings for the widening of the BE acceptance limits.\cr
May be given as character from the choices "EMA", "HC" "FDA" or as an object of
class 'regSet' (see \code{\link{reg_const}}).\cr
Defaults to \code{regulator="EMA"} if missing.\cr
This argument may be given also in lower case if given as character.\cr\cr
The former \code{regulator="ANVISA"} is also allowed but is internally changed 
to "EMA" since ANVISA now recommends the use of EMA regulatory settings. 
Don't use "ANVISA" any longer since it will be removed in next versions.
A corresponding message informs about this.
}
  \item{nsims}{
Number of simulations to be performed to obtain the empirical power.
Defaults to 100 000 = 1e+05.\cr
If simulations are aimed for empirical alpha \code{nsims=1e+06} is recommended.
}
  \item{details}{
If set to \code{TRUE} the computational time is shown as well as the components
for the BE decision.\cr
p(BE-wABEL) is the probability that the CI is within (widened) limits.\cr
p(BE-PE) is the probability that the point estimate is within theta1 ... theta2.\cr
p(BE-ABE) is the simulated probability for the conventional ABE test.
}
  \item{setseed}{
Simulations are dependent on the starting point of the (pseudo) random number 
generator. To avoid differences in power for different runs a \code{set.seed()} 
is issued if \code{setseed=TRUE}, the default. 
}
}
\details{
The methods rely on the analysis of log-transformed data, i.e. assume a 
log-normal distribution on the original scale.\cr\cr
The widened BE acceptance limits will be calculated by the formula\cr
   \code{[lBEL,uBEL] =exp(-+ r_const*sWR)}\cr
with r_const the regulatory constant and sWR the standard deviation of the within
subjects variability of the Reference. r_const=0.76 (~log(1.25)/0.29356) is used 
in case of \code{regulator="EMA"} or \code{regulator="HC"} and in case of 
\code{regulator="FDA"} r_const=0.89257...(=log(1.25)/0.25).
If the CVwR of the Reference is < CVswitch=0.3 the conventional ABE limits 
apply (mixed procedure).\cr\cr 
In case of \code{regulator="EMA"} a cap is placed on the widened limits if 
CVwr>0.5, i.e. the widened limits are held at value calculated for CVwR=0.5.
In case of \code{regulator="HC"} the capping is done such that the acceptance
limits are 0.6666 ... 1.5 at maximum.\cr\cr
The former inofficial regulatory settings for \code{regulator="ANVISA"} are now
covered by \code{regulator="EMA"}.\cr\cr
The simulations are done via the distributional properties of the statistical
quantities necessary for deciding BE based on widened ABEL.\cr
For more details see a document "Implementation_scaledABE_simsVx.yy.pdf" in the 
/doc subdirectory of the package.\cr\cr
Function \code{power.scABEL()} implements the simulation via distributional 
characteristics of the 'key' statistics obtained from the EMA recommended 
evaluation via ANOVA if \code{regulator="EMA"} or if the regulator component 
\code{est_method} is set to "ANOVA" if regulator is an object of class 'regSet'.\cr
Otherwise the simulations are based on the distributional characteristis of the 
'key' statistics obtained from evaluation via intra-subject contrasts (ISC), 
as recommended by the FDA.\cr\cr
Function \code{power.scABEL2()} is soleyly based on simulation using the 
distributional characteristics of the 'key' statistics obtained from evaluation 
via ISC. This function is deprecated.
}
\value{
Returns the value of the (empirical) power if argument \code{details=FALSE}.\cr\cr
Returns a named vector if argument \code{details=TRUE}.\cr
p(BE) is the power, p(BE-wABEL) is the power of the widened ABEL criterion alone 
and p(BE-pe) is the power of the criterion 'point estimat within acceptance 
range' alone. p(BE-ABE) is the power of the conventional ABE test given for 
comparative purposes.
}
\references{
\enc{Tóthfalusi}{Tothfalusi}, L, \enc{Endrényi}{Endrenyi}, L.\cr
Sample Sizes for Designing Bioequivalence Studies for Highly Variable Drugs\cr
J Pharm Pharmaceut Sci. 2011;15(1):73-84. \cr
\href{http://ejournals.library.ualberta.ca/index.php/JPPS/article/download/11612/9489}{free download}
}
\author{
D. Labes
}
\note{
In case of \code{regulator="FDA"} the (empirical) power is only approximate since
the BE decision method is not exactly what is expected by the FDA. But the "Two
Laszlos" state that the scABEL method should be 'operational equivalent' to the
FDA method.\cr
To get the power for the FDA favored method via linearized scaled ABE criterion
use function \code{\link{power.RSABE}}.\cr\cr
In case of \code{regulator="HC"} in function \code{power.scABEL()}, based on
ISC, power is also only approximative since the Health Canada recommends an 
evaluation via mixed model approach. This could only implemented via 
subject data simulations which are very time consuming. But ISC may be a good 
substitute.
}
\section{Warning }{
Cross-validation of the simulations as implemented here and via the 'classical' 
subject data simulation have shown somewhat unsatisfactory results for the 
2x3x3 design if the variabilities for Test and Reference are different.\cr
The function \code{power.scABEL()} therefore gives a warning if calculations 
with different CVwT and CVwR are requested for the 2x3x3 partial replicate design.
For more details see the above mentioned document "Implementation_scaledABE_simsVy.xx.pdf"\cr\cr 
\code{power.scABEL2()} is deprecated and will be removed in future. A corresponding
warning is thrown if this function is used.
}
\seealso{
\code{\link{sampleN.scABEL}, \link{power.RSABE}, \link{reg_const}}
}
\examples{
# using all the defaults:
# design="2x3x3", EMA regulatory settings
# PE constraint 0.8-1.25, cap on widening if CV>0.5
# true ratio =0.90, 1E+6 simulations
power.scABEL(CV=0.4, n=29)
# should give:
# Unbalanced design. n(i)=10/10/9 assumed.
# [1] 0.66113
#
# with details=TRUE to view the computational time and components
power.scABEL(CV=0.5, n=54, theta0=1.15, details=TRUE)
# should give (times may differ depending on your machine):
# 1e+05sims. Time elapsed (sec): 0.07
# 
#      p(BE) p(BE-wABEL)    p(BE-pe)   p(BE-ABE) 
#    0.81727     0.82078     0.85385     0.27542
#
# exploring pure ABEL with the EMA regulatory constant
# (without mixed method, without capping, without pe constraint)
rs <- reg_const("EMA")
rs$CVswitch  <- 0
rs$CVcap     <- Inf
rs$pe_constr <- FALSE
power.scABEL(CV=0.5, n=54, theta0=1.15, regulator=rs)
# should give
# [1] 0.8519
}
