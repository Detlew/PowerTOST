\encoding{utf-8}
\name{reg_const}
\alias{reg_const}
\title{
Constructor of an object with class 'regSet' containing the regulatory settings for ABEL
}
\description{
This function may be used to define regulatory settings not implemented yet
in PowerTOST.
}
\usage{
reg_const(regulator, r_const, CVswitch, CVcap, pe_constr)
}
\arguments{
  \item{regulator}{
Name of the regulatory body as a string. Implemented settings are for \code{"EMA"},
\code{"FDA"}, \code{"HC"}, and \code{"GCC"}.\cr
The former (inofficial) settings for \code{"ANVISA"} are covered by the EMA settings.\cr
In case of \code{regulator="USER"} the other arguments must be given.
Otherwise, they may be missing.
}
  \item{r_const}{
Regulatory constant.
}
  \item{CVswitch}{
CV to switch to the widened acceptance limits.
}
  \item{CVcap}{
CV for capping the widening of the acceptance limits.
}
  \item{pe_constr}{
Logical. Shall pe constraint be applied? Defaults to \verb{TRUE}.
}
}
\value{
Returns an object of class 'regSet', a list with components
  \item{name }{Name of the settings}
  \item{CVswitch }{see arguments}
  \item{r_const }{Regulatory constant}
  \item{CVcap }{see arguments}
  \item{pe_constr}{see arguments}
  \item{est_method}{\code{"ANOVA"} or \code{"ISC"}}
Class 'regSet' has a S3 print method.\cr\cr
The component \verb{est_method} is automatically set to \code{"ANOVA"}, except for
\code{regulator="FDA"} or \code{regulator="HC"} where \code{"ISC"} is used.
}
\author{
D. Labes
}
\note{
The former inofficial regulatory settings for \code{regulator="ANVISA"} are covered by \code{regulator="EMA"} (BEBA Forum).\cr\cr
The settings for CVcap of Health Canada (\code{regulator="HC"}) were chosen in such a way
that the limits of the acceptance range are capped nearly exact to 1/1.5
up to 1.5. Literally it is given rounded to 3 significant digits.
}

\references{
BEBA Forum. May 2016. \href{http://forum.bebac.at/mix_entry.php?id=16319}{online}

Health Canada, Therapeutic Products Directorate. \emph{Comparative Bioavailability Standards: Formulations Used for Systemic Effects, 2.1.1.8 Highly variable drug products} Ottawa, 08 June 2018. \href{https://www.canada.ca/en/health-canada/services/drugs-health-products/drug-products/applications-submissions/guidance-documents/bioavailability-bioequivalence/comparative-bioavailability-standards-formulations-used-systemic-effects.html#a2.1.1.8}{online}
}

\examples{
# to retrieve the EMA settings
reg_const("EMA")
# to define the old ANVISA settings
reg      <- reg_const("USER", r_const = 0.76, CVswitch = 0.4, CVcap = 0.5)
reg$name <- "Old ANVISA"
# Use reg as argument in the scaled ABEL power / sample size functions for
}
