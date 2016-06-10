\name{CI.RatioF}
\alias{CI.RatioF}
\title{
1-2*alpha Fieller confidence interval given point est., CV (,CVb) and n
}
\description{
Utility function to calculate the 1-2*alpha Fieller CI's given point est., 
CV (, CVb) and n for the parallel group and 2x2 crossover.
}
\usage{
CI.RatioF(alpha = 0.025, pe, CV, CVb, n, design = c("2x2", "parallel"))
}
\arguments{
  \item{alpha}{
Type I error probability, aka significance level. \cr
Defaults here to 0.025 because this function is intended for studies
with clinical endpoints. 
}
  \item{pe}{
point estimate of ratio T/R.
}
  \item{CV}{
Coefficient of variation as ratio. In case of \code{design="parallel"} this is
the CV of the total variability, in case of \code{design="2x2"} the 
intra-subject CV.
}
  \item{CVb}{
CV of the between-subject variability. Only necessary for \code{design="2x2"}.
}
  \item{n}{
Total number of subjects if a scalar is given.\cr
Number of subjects in (sequence) groups if given as vector.
}
  \item{design}{
A character string describing the study design.\cr 
\code{design="parallel"} or \code{design="2x2"} allowed for a two-parallel 
group design or a classical TR/RT crossover design.
}
}
\details{
The CV(within) and CVb(etween) in case of \code{design="2x2"} are obtained 
via an appropriate ANOVA from the error term and from the difference 
\code{(MS(subject within sequence)-MS(error))/2}.
}
\value{
Returns the 1-2*alpha confidence interval.
}
\note{
The function assumes an evaluation using un-transformed data.\cr
The function assumes equal variances in case of \code{design="parallel"}.\cr
The formula implemented covers balanced and unbalanced designs.\cr\cr
Note that when the mean of the denominator of the ratio is close to zero, 
confidence intervals might be degenerated and are returned as NA.
In that case a warning is issued.\cr\cr
Whether the function vectorizes properly is not thoroughly tested.\cr\cr
This function is intended for studies with clinical endpoints. 
In such studies the 95\% confidence intervals are usually used for equivalence 
testing. Therefore alpha defaults here to 0.025.\cr
See CPMP/EWP/482/99 "Points to consider on switching between superiority and 
non-inferiority" EMEA, London (2000).
}
\author{
D. Labes
}
\references{
Locke C.S.\cr
"An exact confidence interval from untransformed data for the ratio of 
two formulation means."\cr
J Pharmacokinet Biopharm. 12(6):649-55 (1984)

Hauschke D., Steinijans V. and Pigeot I.\cr
"Bioequivalence Studies in Drug Development"\cr
Chapter 10., John Wiley & Sons, Chichester (2007)
}
\seealso{
\code{\link{CI.BE}}, \code{\link{power.RatioF}}
}
\examples{
# 95\% Fieller CI for the 2x2 crossover
CI.RatioF(pe=1.05,CV=0.3,CVb=0.6, n=24)
}

