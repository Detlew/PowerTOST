\encoding{utf-8}
\name{power.TOST.sds}
\alias{power.TOST.sds}
\title{
Power calculation of the BE decision with models incorporating groups
}
\description{
The power is calculated via subject data simulations.\cr
Three models are implemented:\cr
gmodel==1 is full FDA model but only for testing group-by-treatment interaction\cr
gmodel==2 is full FDA model but without group-by-treatment interaction\cr
gmodel==3 is model with pooled groups, i.e. without any group term
}
\usage{
power.TOST.sds(alpha = 0.05, theta1, theta2, theta0, CV, n, 
               design = c("2x2", "2x2x2", 2x3x3", "2x2x4", "2x2x3"),
               design_dta = NULL, grps = 2, ngrp = NULL, gmodel = 2,
               nsims = 1e+05, details = FALSE, setseed = TRUE, progress)
}
\arguments{
  \item{alpha}{
Type I error probability, significance level. Conventionally mostly set to 0.05.
}
  \item{theta1}{
Lower BE limit. Defaults to 0.8 if not given explicitely.
}
  \item{theta2}{
Upper BE limit. Defaults to 1.25 if not given explicitely.
}
  \item{theta0}{
\sQuote{True} or assumed T/R ratio. \cr
Defaults to 0.95 if not given explicitly.
}
  \item{CV}{
Intra-subject coefficient(s) of variation as ratio (not percent).
  \itemize{
    \item If given as a scalar (\code{length(CV)==1}) the \emph{same} CV of Test
          and Reference is assumed (homoscedasticity, \code{CVwT==CVwR}).
    \item If given as a vector (\code{length(CV)==2}), \emph{i.e.}, assuming
          heteroscedasticity, the CV of the Test \strong{must} be given in 
          \code{CV[1]} and the one of the Reference in the \code{CV[2]}.
}
}
  \item{n}{
Number of subjects under study.\cr
May be given as vector. In that case it is assumed that \code{n} contains the number
of subjects in the sequence groups.\cr
If \code{n} is given as single number (total sample size) and this number is not 
divisible by the number of sequences of the design an unbalanced design is 
assumed. A corresponding message is thrown showing the numbers of subjects 
in sequence groups.\cr
Attention! In case of the \code{"2x2x3"} (TRT|RTR) design the order of sample sizes is important 
if given as vector. \code{n[1]} is for sequence group 'TRT' and \code{n[2]} is for 
sequence group 'RTR'. 
}
  \item{design}{
Design of the study to be planned.\cr
\code{"2x2"} or \code{"2x2x2"} is the conventional cross-over design.\cr
\code{"2x3x3"} is the partial replicate design (TRR|RTR|RRT).\cr
\code{"2x2x4"} is the full replicate design with 2 sequences and 4 periods.\cr
\code{"2x2x3"} is the 3-period design with sequences TRT|RTR.\cr
Defaults to \code{design="2x2"}.
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
  \item{grps}{
Number of (logistical) groups. Defaults to 2.
}
  \item{ngrp}{
Vector of number of subjects in groups.
}
  \item{gmodel}{
\code{gmodel=1} is full FDA model for testing grp x treatment interaction
followed by an evaluation with gmodel=2 if the interaction is not significant 
at p=0.1. If significant gmodel=3 with the data from the largest group is used
for the BE decision.\cr
\code{gmodel=2} is full FDA model but without group by treatment interaction,
i.e. poolability of groups assumed but groups considered in the model.\cr
\code{gmodel=3} is model with pooled groups, i.e. without any group term.\cr
Defaults to \code{gmodel=2}.
}
  \item{nsims}{
Number of simulations to be performed to obtain the empirical power.
Defaults to 100,000 = 1e+05.\cr
If simulations are aimed for empirical alpha \code{nsims=1e+06} is recommended.
}
  \item{details}{
If set to \code{TRUE} the computational time is shown.
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
The power is calculated via subject data sims.\cr
The evaluation of BE is done via 1-2*alpha confidence interval using classical ANOVA
for the models with group effects.\cr
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
Returns the value of the (empirical) power
}
\references{
Schütz H. \emph{Multi-Group Studies in Bioequivalence. To pool or not to pool?}\cr
Presentation at BioBriges 2018, Prague. \url{https://bebac.at/lectures/Prague2018.pdf}
}
\author{
D. Labes
}
\note{
The run time of the function may be relatively long.\cr
Take a cup of coffee and be patient.
}
\examples{
# power for gmodel=2, 2x2 crossover, grps=3 with even number of subjects
power.TOST.sds(CV=0.2, n=18, grps=3)
# gives [1] 0.78404
# without considering groups
power.TOST.sds(CV=0.2, n=18, gmodel=3)
# gives [1] 0.7887
}
