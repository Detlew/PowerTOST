#---------------------------------------------------------------------------
# unified function
# chooses the power function according to regulator$est_method
#
# author dlabes
#---------------------------------------------------------------------------
power.scABEL <- function(alpha=0.05, theta1, theta2, theta0, CV, n,   
                         design=c("2x3x3", "2x2x4", "2x2x3"), regulator,
                         sdsims = FALSE, nsims=1E5, details=FALSE, setseed=TRUE)
{
  # design must be checked outside
  desi <- match.arg(design)
  # check regulator
  if (missing(regulator)) regulator <- "EMA"
  reg  <- reg_check(regulator)

  if (reg$est_method!="ISC"){
    if (sdsims) {
      # power via subject data sims
      # additional to the args of power.scABEL1 / power.scABEL2 
      # are the argumets design_dta, progress
      r <- power.scABEL.sdsims(alpha, theta1, theta2, theta0, CV, n, design=desi,  
                               design_dta=NULL, regulator=reg, nsims, details, 
                               setseed)
    } else {
      # power via key statistic sims with empirical adaptions to obtain
      # better agreement to sims based on subject data
      r <- power.scABEL1(alpha, theta1, theta2, theta0, CV, n, design=desi, reg, 
                         nsims, details, setseed)
    }
  } else {
    # power via key statistic sims
    r <- power.scABEL2(alpha, theta1, theta2, theta0, CV, n, design=desi, reg, 
                       nsims, details, setseed)
  } 
  r
}  # end function
