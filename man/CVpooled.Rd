\name{CVpooled}
\alias{CVpooled}
\alias{print.CVp}
\title{Pooled CV from several studies}
\description{
This function pools CVs of several studies.
}
\usage{
CVpooled(CVdata, alpha = 0.2, logscale=TRUE, robust = FALSE)
\method{print}{CVp}(x, digits=4, verbose=FALSE, \dots)
}
\arguments{
  \item{CVdata}{
A data.frame that must contain the columns \code{CV}, \code{n} and \code{design} 
where CV are the error CVs from the studies, n the number of subjects 
and design is a character string describing the study design.\cr
See \code{known.designs()} for designs covered in this package.\cr
If the design column is missing the classical 2x2 crossover is assumed for each study.\cr
A message is displayed under that circumstances.\cr\cr
A data.frame that contains the columns \code{CV} and giving the degrees of freedom 
\code{df} directly is also accepted as \code{CVdata}.
  }
  \item{alpha}{
Error probability for calculating an upper confidence limit of the pooled CV.\cr
Recommended 0.2-0.25 for use in subsequent sample size estimation. \cr
See f.i one of H. Schuetz lectures \url{http://bebac.at/lectures/MU2010-CD2.pdf}
  }
  \item{logscale}{
Defaults to TRUE.\cr
Should the calculations be done for log-transformed data?  
  }
  \item{robust}{
Defaults to FALSE.\cr
Set to \code{TRUE} will use the degrees of freedom according to the 'robust' evaluation
(aka Senn's basic estimator). These df's are calculated as \code{n-seq}.\cr
They are also often more appropriate if the CV comes from a 'true' mixed model 
evaluation (FDA model for average bioequivalence).\cr
See \code{known.designs()$df2} for the designs covered in this package.
  }
  \item{x}{An object of class \code{"CVp"}.}
  \item{digits}{Number of digits for CV.}
  \item{verbose}{Defaults to \code{FALSE}. Prints only the pooled CV and the df.\cr
  If set to \code{TRUE} the upper confidence limit is also printed.
  }
  \item{\dots}{More args to print(). None used.}
}
\details{
The pooled CV is obtained from the weighted average of the error variances obtained from 
the CVs of the single studies, weights are the df (degrees of freedom).\cr
If only \code{n} is given in the input \code{CVdata}, the df's are calculated via
the formulas given in \code{known.designs()}. If both \code{n} and \code{df} are given
the \code{df} column precedes.\cr\cr
If \code{logscale=TRUE} the error variances are obtained via function \code{CV2se()}.
Otherwise the pooled CV is obtained via pooling the CV^2.
}
\value{
A list of class \code{"CVp"} with components
  \item{CV}{value of the pooled CV}
  \item{df}{pooled degrees of freedom}
  \item{CVupper}{upper confidence interval of the pooled CV}
  \item{alpha}{input value}
The class \code{"CVp"} has a S3 methods \code{print.CVp}.
}
\references{
H. Schuetz lectures about sample size challenges \cr
at \url{http://bebac.at/lectures.htm}.\cr\cr
Patterson, Jones\cr
"Bioequivalence and Statistics in Clinical Pharmacology"\cr
Chapter 5.7 "Determining Trial Size"\cr
Chapman & Hall/CRC, Boca Raton 2006
}
\author{
D. Labes
}
\note{
The calculations for \code{logscale=FALSE} are not described in the references.
They are implemented by analogy to the case via log-transformed data.\cr
The calculations are based on a common variance of Test and Reference formulations\cr
in replicate crossover studies or parallel group study, respectively.
}
\section{Warning}{
Pooling of CVs from parallel-group and cross-over designs does not make any sense.\cr
Also the function does not throw an error if you do so.
}
\seealso{
  \code{\link{known.designs}, \link{CVfromCI} }
}
\examples{
# some data:
# the values for AUC, study 1 and study 2 are Example 3 of H. Schuetz lecture
CVs <- ("
 PKmetric | CV   | n  |design|source
    AUC   | 0.20 | 24 | 2x2  | study 1
    Cmax  | 0.25 | 24 | 2x2  | study 1
    AUC   | 0.30 | 12 | 2x2  | study 2
    Cmax  | 0.31 | 12 | 2x2  | study 2
    AUC   | 0.25 | 12 | 2x2x4| study 3 (replicate)
")
txtcon <- textConnection(CVs)
CVdata <- read.table(txtcon, header=TRUE, sep="|", strip.white=TRUE, as.is=TRUE)
close(txtcon)

# evaluation of the AUC CVs
CVsAUC <- subset(CVdata, PKmetric=="AUC")
CVpooled(CVsAUC, alpha=0.2, logscale=TRUE)
# df of the 'robust' evaluation
CVpooled(CVsAUC, alpha=0.2, logscale=TRUE, robust=TRUE)
#print also the upper CL, data example 3
CVsAUC3 <- subset(CVsAUC,design != "2x2x4")
print(CVpooled(CVsAUC3, alpha=0.2, robust=TRUE), digits=3, verbose=TRUE)
# will give the output:
#Pooled CV = 0.235 with 32 degrees of freedom (robust df's)
#Upper 80\% confidence limit of CV = 0.266
}