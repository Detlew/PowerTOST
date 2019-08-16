\encoding{utf-8}
\name{sampleN.RSABE2L.sdsims}
\alias{sampleN.RSABE2L.sdsims}
\alias{sampleN.RSABE2L.sds}
\title{
Sample size estimation for BE Decision via Reference Scaled ABE
}
\description{
These functions performs the sample size estimation of the BE decision via 
the reference scaled ABE based on \bold{subject data} simulations.
Implemented are the methods ABEL, Hyslop and \sQuote{exact} (see the References in \code{\link{power.RSABE2L.sdsims}}).\cr
The estimation method of the key statistics needed to perform the \acronym{RSABE} decision 
is the usual \acronym{ANOVA}.\cr
This function has an alias sampleN.RSABE2L.sds().
}
\usage{
sampleN.RSABE2L.sdsims(alpha = 0.05, targetpower = 0.8, theta0,
                       theta1, theta2, CV,
                       design = c("2x3x3", "2x2x4", "2x2x3"),
                       SABE_test = "exact", regulator, nsims=1e5,
                       nstart, imax = 100, print = TRUE,
                       details = TRUE, setseed = TRUE, progress)
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
\sQuote{True} or assumed T/R ratio. \cr
Defaults to 0.90 according to the two \enc{Laszlós}{Laszlos} if not given explicitly.
}
  \item{theta1}{
Conventional lower ABE limit to be applied in the mixed procedure if 
\code{CVsWR <= CVswitch}.\cr
Also Lower limit for the point estimate constraint.\cr
Defaults to 0.8 if not given explicitly.
}
  \item{theta2}{
Conventional upper ABE limit to be applied in the mixed procedure if 
\code{CVsWR <= CVswitch}. Also upper limit for the point estimate constraint.\cr
Defaults to 1.25 if not given explicitly.
}
  \item{CV}{
Intra-subject coefficient(s) of variation as ratio (not percent).
  \itemize{
    \item If given as a scalar (\code{length(CV)==1}) the \emph{same} CV of Test
          and Reference is assumed (homoscedasticity, \code{CVwT==CVwR}).
    \item If given as a vector (\code{length(CV)==2}), \emph{i.e.}, assuming
          heteroscedasticity, the CV of the Test \strong{must} be given in \code{CV[1]} and the one of the Reference in the \code{CV[2]}.
  }
}
  \item{design}{
Design of the study to be planned.\cr
\code{"2x3x3"} is the partial replicate design.\cr
\code{"2x2x4"} is a full replicate design with 2 sequences and 4 periods.\cr
\code{"2x2x3"} is a full replicate design with 2 sequences and 3 periods.\cr
Defaults to \code{design="2x3x3"}. Details are given the section about Designs.
}
  \item{SABE_test}{
This argument specifies the test method to be used for the reference scaled
\acronym{ABE} decision.\cr
Default is the \code{"exact"} \sQuote{ncTOST} method of the two \enc{Laszlós}{Laszlos}. Other choices are \code{"ABEL"}, \code{"Hyslop"} and \code{"FDA"}. See Details.\cr
This argument may be given also in lower case.
}
  \item{regulator}{
Regulatory settings for the widening of the BE acceptance limits.\cr
May be given as character \code{"EMA"} or as an object of
class 'regSet' (see \code{\link{reg_const}}).\cr
Defaults to \code{regulator="EMA"} if missing.\cr
This argument may be given also in lower case if given as character.\cr
If given as object of class 'regSet' the component \code{est_method} can not be \code{"ISC"}. 
}
  \item{nsims}{
Number of simulations to be performed to obtain the (empirical) power.
The default value 100,000 = 1e+5 is usually sufficient. Consider to rise 
this value if theta0<=0.85 or >=1.25. But see the warning section.
}
  \item{nstart}{
Set this to a start for the sample size search if a previous run failed.\cr
After reworking the start n in version 1.1-05 rarely needed.
}
  \item{imax}{
Maximum number of steps in sample size search. Defaults to 100.
}
  \item{print}{
If \verb{TRUE} (default) the function prints its results. If \verb{FALSE} only the result data.frame will be returned.
}
  \item{details}{
If set to \code{TRUE} (default), the steps during sample size search are shown.
}
  \item{setseed}{
Simulations are dependent on the starting point of the (pseudo) random number 
generator. To avoid differences in power for different runs a \code{set.seed(123456)} 
is issued if \code{setseed=TRUE}, the default. 
}
  \item{progress}{
Should a progressbar be shown? Defaults to \code{TRUE} if missing and nsims >5E5.
}
}
\details{
The methods rely on the analysis of log-transformed data, \emph{i.e.}, assume a 
log-normal distribution on the original scale.\cr\cr
The widened BE acceptance limits will be calculated by the formula\cr
\verb{  [L, U] = exp(± r_const * sWR)}\cr
with \code{r_const} the regulatory constant and \code{sWR} the standard deviation of the within
subjects variability of the Reference. \code{r_const = 0.76} (~log(1.25)/0.29356) is used 
in case of \code{regulator="EMA"}.
If the CVwR of the Reference is < CVswitch=0.3 the conventional ABE limits 
apply (mixed procedure).\cr 
In case of \code{regulator="EMA"} a cap is placed on the widened limits if 
CVwr>0.5, \emph{i.e.}, the widened limits are held at value calculated for CVwR=0.5.\cr\cr
The simulations are done by simulating subject data (all effects fixed except the
residuals) and evaluating these data via ANOVA of all data to get the point estimate
of T vs. R along with its 90\% CI and an ANOVA of the data under R(eference) only
to get an estimate of s2wR.
}
\value{
Returns a data.frame with the input settings and sample size results.\cr
The \code{Sample size} column contains the total sample size.\cr
The \code{nlast} column contains the last \code{n} value. May be useful for restarting.
}
\author{
H. \enc{Schütz}{Schuetz}
}
\note{
We are doing the sample size estimation only for balanced designs since the 
break down of the total subject number in case of unbalanced sequence groups
is not unique. Moreover the formulas used are only for balanced designs.\cr
The minimum sample size is 6, even if the power is higher than the intended 
targetpower.\cr\cr
Subject simulations are relatively slow. Thus be patient and go for a cup of coffee if you use this function with high sample sizes!
}
\section{Designs}{
Although some designs are more \sQuote{popular} than others, sample size estimations are valid for \emph{all} of the following designs:
  \tabular{ll}{
    \code{"2x2x4"} \tab TRTR | RTRT\cr
    \tab TRRT | RTTR\cr
    \tab TTRR | RRTT\cr
    \code{"2x2x3"} \tab TRT | RTR\cr
    \tab TRR | RTT\cr
    \code{"2x3x3"} \tab TRR | RTR | RRT
  }
}
\section{Warning }{
The sample size estimation for very extreme theta0 (<0.83 or >1.21) may be very 
time consuming and will eventually also fail since the start values chosen are 
not really reasonable in that ranges.\cr
If you really need sample sizes in that range be prepared to restart the sample 
size estimation via the argument nstart.\cr
Since the dependence of power from n is very flat in the mentioned region you may 
also consider to adapt the number of simulations not to tap in the simulation 
error trap.
}
\seealso{
\code{\link{power.RSABE2L.sdsims}}, \code{\link{sampleN.scABEL}}, \code{\link{reg_const}}
}
\examples{
# using the defaults:
# partial replicate design, targetpower=80\%,
# true assumed ratio = 0.90, 1E+5 simulated studies
# ABE limits, PE constraint 0.8 - 1.25
# EMA regulatory settings
# compare results
\dontrun{
CV     <- 0.4
method <- c("exact", "abel", "hyslop", "fda")
res    <- data.frame(SABE_test = c("ncTOST", "ABEL",
                                   "Hyslop", "FDA"),
                     n = NA, power = NA)
for (i in 1:nrow(res)) {
  res[i, 2:3] <- sampleN.RSABE2L.sdsims(CV = CV,
                                        SABE_test = method[i],
                                        details = FALSE,
                                        print = FALSE)[8:9]
}
print(res, digits = 4, row.names = FALSE)}
# should result in a sample size n=48 with all methods,
# power=0.8197 (ncTOST), 0.8411 (ABEL), 0.8089 (Hyslop), 0.8113 (FDA)
}
