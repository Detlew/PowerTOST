\encoding{utf-8}
\name{power.scABEL.sdsims}
\alias{power.scABEL.sdsims}
\title{
(Empirical) Power of BE decision via scaled (widened) BE acceptance limits
}
\description{
These function performs the power calculation of the BE decision via 
scaled (widened) BE acceptance limits based on \bold{subject data} simulations.
}
\usage{
power.scABEL.sdsims(alpha = 0.05, theta1, theta2, theta0, CV, n, 
                    design = c("2x3x3", "2x2x4", "2x2x3"), design_dta=NULL,
                    regulator, nsims = 1e+05, details = FALSE, setseed = TRUE, 
                    progress)
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
2x3x3 is the partial replicate design (TRR|RTR|RRT).\cr
2x2x4 is the full replicate design with 2 sequences and 4 periods (TRTR|RTRT).\cr
2x2x3 is the 3-period design with sequences TRT|RTR.\cr
Defaults to \code{design="2x3x3"}.
}
  \item{design_dta}{
Alternatively to using the arguments \code{design} and \code{n} the design may 
be defined via a data.frame with columns \code{subject, sequence, period} and
\code{tmt}. This feature is experimental in the sense that the data.frame is
not checked for complying with the assumed structure.\cr
If you use the argument \code{design_dta} you don't need to specify the arguments
\code{design} and \code{n}.\cr
The default \code{design_dta = NULL} means that \code{design} and \code{n} are used
for the internal construction of the design data.frame.
}
  \item{regulator}{
Regulatory settings for the widening of the BE acceptance limits.\cr
May be given as character == "EMA" or as an object of
class 'regSet' (see \code{\link{reg_const}}).\cr
Defaults to \code{regulator="EMA"} if missing.\cr
This argument may be given also in lower case if given as character.\cr\cr
If given as object of class 'regSet' the component \code{est_method} can not be "ISC". 
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
  \item{progress}{
Should a progressbar be shown? Defaults to TRUE if missing and nsims >5E5.
}
}
\details{
The methods rely on the analysis of log-transformed data, i.e. assume a 
log-normal distribution on the original scale.\cr\cr
The widened BE acceptance limits will be calculated by the formula\cr
   \code{[lBEL,uBEL] =exp(-+ r_const*sWR)}\cr
with r_const the regulatory constant and sWR the standard deviation of the within
subjects variability of the Reference. r_const=0.76 (~log(1.25)/0.29356) is used 
in case of \code{regulator="EMA"}.
If the CVwR of the Reference is < CVswitch=0.3 the conventional ABE limits 
apply (mixed procedure).\cr 
In case of \code{regulator="EMA"} a cap is placed on the widened limits if 
CVwr>0.5, i.e. the widened limits are held at value calculated for CVwR=0.5.\cr\cr
The simulations are done by simulating subject data (all effects fixed except the
residuals) and evaluating these data via ANOVA of all data to get the point esimate
of T vs. R along with its 90\% CI and an ANOVA of the data under R(eference) only
to get an estimate of s2wR.\cr
The data.frame with columns  \code{subject, sequence, period} and \code{tmt} 
necessary for evalution of simulated subject data is constructed internally from 
the arguments \code{design} and \code{n} or may be given user defined via the argument 
\code{design_dta}. The last option is usefull if missing data have to be considered
or if designs have to be evaluated which are not in the list of argument 
\code{design}.\cr
This feature is experimental in the sense that the data.frame is not checked 
for complying with the assumed structure.
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
\enc{Tóthfalusi}{Tothfalusi} L, \enc{Endrényi}{Endrenyi} L.\cr
\emph{Sample Sizes for Designing Bioequivalence Studies for Highly Variable Drugs}\cr
J Pharm Pharmaceut Sci. 2011;15(1):73--84. \href{http://ejournals.library.ualberta.ca/index.php/JPPS/article/download/11612/9489}{free download}
}
\author{
D. Labes & B. Lang
}
\note{
The function is mainly intended for crosscheck of \code{power.scABEL()} results.\cr
But may be mandatory for cases where \code{power.scABEL()} results are inaccurate
(low samplesizes and/or heteroscedasticity).\cr
It is relatively slow. The run-time of this function doing 1 Mio sims is between 
~ 7-8 sec for n=12 and ~ 3-4 min for n=120 on a machine with an Intel core i7 processor.\cr
Thus be patient and go for a cup of coffee if you use this function with high 
sample sizes!
}
\seealso{
\code{\link{power.scABEL}, \link{reg_const}}
}
\examples{
# using all the defaults:
# design="2x3x3", EMA regulatory settings
# PE constraint 0.8-1.25, cap on widening if CV>0.5
# true ratio =0.90, 1E+5 simulations
\dontrun{power.scABEL.sdsims(CV=0.4, n=36)}
# should give:
# [1] 0.74321
}
