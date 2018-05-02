#---------------------------------------------------------------------------
# unified function
# chooses the power function according to regulator$est_method
#
# author dlabes
#---------------------------------------------------------------------------
power.scABEL <- function(alpha=0.05, theta1, theta2, theta0, CV, n,   
                         design=c("2x3x3", "2x2x4", "2x2x3"), regulator,
                         nsims=1E5, details=FALSE, setseed=TRUE)
{
  # design must be checked outside
  desi <- match.arg(design)
  # check regulator
  if (missing(regulator)) regulator <- "EMA"
  reg  <- reg_check(regulator)
  pwrfun <- "power.scABEL1"
  if (reg$est_method=="ISC") pwrfun <- "power.scABEL2"
  # next doesn't function if one or more theta's missing
  # r <- do.call(pwrfun,
  #              list(alpha, theta1, theta2, theta0, CV, n, design=desi, reg, 
  #                   nsims, details, setseed))
  if (reg$est_method!="ISC"){
    # power via key statistic sims with empirical adaptions to obtain
    # better agreement to sims based on subject data
    r <- power.scABEL1(alpha, theta1, theta2, theta0, CV, n, design=desi, reg, 
                       nsims, details, setseed)
  } else {
    r <- power.scABEL2(alpha, theta1, theta2, theta0, CV, n, design=desi, reg, 
                       nsims, details, setseed)
  } 
  r
}  # end function
