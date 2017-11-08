\encoding{utf-8}
\name{power.RSABE2L.sdsims}
\alias{power.RSABE2L.sdsims}
\alias{power.RSABE2L.sds}
\title{
(Empirical) Power of BE Decision via Reference Scaled ABE
}
\description{
These function performs the power calculation of the BE decision via 
the reference scaled ABE based on \bold{subject data} simulations.
Implemented are the methods ABEL, Hyslop and 'exact' as described in the 
references.\cr
The estimation method of the key statistics needed to perform the RSABE decision 
is the usual ANOVA.
}
\usage{
power.RSABE2L.sdsims(alpha = 0.05, theta1, theta2, theta0, CV, n, 
                     design = c("2x3x3", "2x2x4", "2x2x3"), design_dta=NULL,
                     SABE_test="exact", regulator, nsims = 1e+05, 
                     details = FALSE, setseed = TRUE, progress)
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
  \item{SABE_test}{
This argument specifies the test method to be used for the reference scaled
ABE decision.\cr
Default is the "exact" ncTOST method of the 2 Laszlo's. Other choices are "ABEL",
"hyslop" and "fda". See Details.
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
p(BE-RSABE) is the probability of a positive outcome of the SABE test.\cr
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
The data.frame with columns  \code{subject, sequence, period} and \code{tmt} 
necessary for evalution of simulated subject data is constructed internally from 
the arguments \code{design} and \code{n} or may be given user defined via the argument 
\code{design_dta}. The last option is usefull if missing data have to be considered
or if designs have to be evaluated which are not in the list of argument 
\code{design}.\cr\cr
The estimation method for obtaining the statistics necessary to perform the
reference scaled ABE decision is the usual ANOVA with effects treatment, period,
sequence and subject within sequence for the evaluation of all data and period,
sequence and subject within sequence for the evaluation of the Reference formulation
data only.\cr\cr
The SABE tests implemented are:\cr
\tabular{ll}{
\code{"exact"} \tab - 'exact' based method of the 2 Laszlo's (see references, 
called there ncTOST)\cr
\code{"ABEL"} \tab - Average bioequivalence with expanding limits\cr
\code{"hyslop"} \tab - BE decision via the linearized RSABE criterion and its 
upper 95\% CI\cr
\code{"fda"} \tab - Hyslop with an additional bias correction term as implemented 
in the SAS code of the FDA Guidance on Progesterone.
}
}
\value{
Returns the value of the (empirical) power if argument \code{details=FALSE}.\cr\cr
Returns a named vector if argument \code{details=TRUE}.\cr
p(BE) is the power, p(BE-RSABE) is the power of using the reference scaled ABE 
alone and p(BE-pe) is the power of the criterion 'point estimat within acceptance 
range' alone. p(BE-ABE) is the power of the conventional ABE test given for 
comparative purposes.
}
\references{
FDA \emph{Draft Guidance on Progesterone}\cr
Recommended Apr 2010. Revised Feb 2011. 
\href{http://www.fda.gov/downloads/Drugs/GuidanceComplianceRegulatoryInformation/Guidances/UCM209294.pdf}{download}

\enc{Tóthfalusi}{Tothfalusi} L, \enc{Endrényi}{Endrenyi} L.\cr
\emph{An Exact Procedure for the Evaluation of Reference-Scaled Average
Bioequivalence}\cr
The AAPS J., 2016 Mar;18(2):476-89\cr
Epub 2016 Jan 29. \doi{10.1208/s12248-016-9873-6}.

\enc{Tóthfalusi}{Tothfalusi} L, \enc{Endrényi}{Endrenyi} L.\cr
\emph{Algorithms for evaluating reference scaled average bioequivalence: power, 
bias, and consumer risk}\cr
Stat Med. 2017 Nov 30;36(27):4378-4390\cr 
Epub 2017 Aug 29. \doi{10.1002/sim.7440}
}
\author{
D. Labes
}
\note{
The function is relatively slow. The run-time of doing 1 Mio sims 
is between ~ 1 upto 6 min for n=12 or n=120 and 1 Mio sims (see the call under 
examples) on a machine with an Intel core i7 processor.\cr
Thus be patient and go for a cup of coffee if you use this function with higher 
sample sizes and aim for calculating type 1 error!
}
\seealso{
\code{\link{power.RSABE}, \link{reg_const}}
}
\examples{
# Not run due to timing policy of CRAN
\dontrun{
# pure EMA settings without mixed procedure, cap on widening and PE constraint
# as in the reference from 2017
reg <- reg_const("EMA")
reg$CVswitch  <- 0
reg$CVcap     <- Inf
reg$pe_constr <- FALSE
reg$name      <- "EMA pure"
power.RSABE2L.sds(CV=0.4, n=12, theta0=exp(0.05), design="2x2x4", regulator=reg, nsims=50000)
# should give:
# [1] 0.46504 (compared to 47.1\% in the 2017 reference)}
}
