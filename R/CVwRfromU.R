# ----------------------------------------------------------------
# function to calculate the CVwR from upper expanded limit (ABEL)
# 
# Author: hschuetz
# ----------------------------------------------------------------

CVwRfromU <- function(U, regulator = "EMA")
{
  if (missing(U))
    stop("Upper expanded limit \'U\' must be given!")
  if (!regulator %in% c("EMA", "HC"))
    stop("Regulator must be \'EMA\' or \'HC\'.")
  reg <- reg_const(regulator)
  # min and max upper limits
  EL <- c(scABEL(CV=reg$CVswitch, regulator=regulator)[["upper"]],
          scABEL(CV=reg$CVcap, regulator=regulator)[["upper"]])
  if (U <= EL[1] || U >= EL[2]) {
    stop(sprintf("Calculation only possible if 1.2500 > U < %.4f!",
                 scABEL(CV=reg$CVcap, regulator=regulator)["upper"]))
  }
  # same precision as U
  CVwR <- signif(sqrt(exp((log(U)/reg$r_const)^2)-1), 5)
  return(CVwR)
}
# ----------------------------------------------------------------
# alias to CVwRfromUL
U2CVwR <- function(U, regulator = "EMA")
{
  CVwRfromU(U, regulator = "EMA")
}
