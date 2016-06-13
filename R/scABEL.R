# ---------------------------------------------------------------------------
# helper function for the scABEL regulatory settings
# constructor of an object of class 'regSet'
reg_const <- function(regulator, r_const, CVswitch, CVcap, pe_constr)
{
  regulator <- toupper(regulator)
  if (regulator=="USER"){
    if (missing(CVswitch) | missing(CVcap) | missing(r_const)){
      stop("r_const, CVswitch and CVcap must be given.")
    }
    if (missing(pe_constr)) pe_constr <- TRUE
    r <- list(name="USER", CVswitch=CVswitch, r_const=r_const, CVcap=CVcap,
              pe_constr=pe_constr)
  }
  else if (regulator=="FDA"){
    r <- list(name="FDA", CVswitch=0.3, r_const=log(1.25)/0.25, CVcap=Inf,
              est_method="ISC")
  } 
  # else if (regulator=="ANVISA"){
  #   # same regulatory const. as EMA but
  #   # switch to widened limits if CVRef>=40% (inofficial)
  #   # now no longer valid, same as EMA settings
  #   r <- list(name="ANVISA", CVswitch=0.4, r_const=log(1.25)/CV2se(0.3), CVcap=0.5)
  #   r <- list(name="ANVISA", CVswitch=0.3, r_const=0.76, CVcap=0.5)
  # }
  else if (regulator=="EMA"){
    # r_const taken literally from BE guideline
    r <- list(name="EMA", CVswitch=0.3, r_const=0.76, CVcap=0.5)
  }
  else if (regulator=="HC"){
    # r_const taken literally from guideline, 
    # CVcap chosen to obtain nearly exact 1.50 as upper limit, 
    # literally it was given as 57.4%
    r <- list(name="HC", CVswitch=0.3, r_const=0.76, 
              CVcap=0.57382,  # se2CV(log(1.5)/0.76) = 0.57381995
              est_method="ISC") 
  } else {
    stop("Unknown regulator.")
  }
  class(r) <- "regSet"
  # default is with pe constraint
  if (is.null(r$pe_constr)) r$pe_constr <- TRUE
  if (is.null(r$est_method)) r$est_method <- "ANOVA"
  r
} 

# ---------------------------------------------------------------------------
# helper function to check regulatory settings
# returns the regulatory settings as an object of class 'regSet'
# regulator may be a character string or an object of class 'regSet'
reg_check <- function(regulator, choices=c("EMA", "HC", "FDA", "ANVISA"))
{
  #browser()
  if (class(regulator)=="character"){
    reg <- toupper(regulator)
    reg <- match.arg(reg, choices)
    if (reg=="ANVISA"){
      message("Former (inofficial) ANVISA regulatory settings changed to EMA settings.")
      message("Don't use 'ANVISA' any longer since it will be removed in next versions.\n")
      reg <- "EMA"
    }
    reg <- reg_const(reg)
  } else if (class(regulator)=="regSet") {
    reg <- regulator
    # more checks?
  } else {
    stop("Arg. regulator has to be character or an object of class 'regSet'.")
    reg <- NULL
  }
  reg
}  
# -------------------------------------------------------------------------
# function to nicely print the ABEL regulatory settings
print.regSet <- function(x, ...)
{
  if(x$name=="USER"){
    cat(x$name, "defined regulatory settings\n")
  } else {
    cat(x$name, "regulatory settings\n")
  }
  cat("- CVswitch            =", x$CVswitch,"\n")
  if (!is.null(x$CVcap)) {
    if (is.finite(x$CVcap)){
      cat("- cap on scABEL if CVw(R) > ", x$CVcap,"\n",sep="")
    } else {
      cat("- no cap on scABEL\n", sep="")
    }  
  }
  cat("- regulatory constant =", x$r_const,"\n")
  if (is.null(x$pe_constr)) x$pe_constr <- TRUE
  if (x$pe_constr) {
    cat("- pe constraint applied")
  } else {
    cat("- no pe constraint")
  }
  cat("\n")
}

# -------------------------------------------------------------------------
# function to calculate the (widened) ABE acceptance limits
# -------------------------------------------------------------------------
scABEL <- function(CV, regulator)
{
  if (missing(regulator)) regulator <-"EMA"
  
  rc <- reg_check(regulator)
  CVcap    <- rc$CVcap
  CVswitch <- rc$CVswitch
  r_const  <- rc$r_const
  
  # upper acceptance limit, 1e-10 to assure 1.25 for CV nearly equal CVswitch
  ret <- ifelse(CV <= (CVswitch + 1e-10), 1.25, exp(r_const*CV2se(CV)))
  # cap on widening
  ret <- ifelse(CV>CVcap, exp(r_const*CV2se(CVcap)), ret)
  # lower acceptance limit is set to 1/upper
  if (length(CV)>1){
    ret <- cbind(1/ret, ret)
    colnames(ret) <- c("lower", "upper")
  } else {
    ret <- c(1/ret, ret)
    names(ret) <- c("lower", "upper")
  }
  ret
}

# --------------------------------------------------------------------------
# function to calculate "leveling-off" ABEL according to Karalis et al. 2011
# Eur J Pharm Sci, 44, 497-505
# --------------------------------------------------------------------------
scABEL_LO <- function(CV)
{
  gamma <- 0.0336 # Karalis et al.
  sw0   <- 0.3853

  gamma <- 0.03361 # Own fit with stepsize 0.01 in CVwR
  sw0   <- 0.38535
  
  beta  <- scABEL(CV=0.5)["upper"]
  # sigmoidal
  uppr <- 1.25 + (beta - 1.25)/(1 + exp(-(CV2se(CV)-sw0)/gamma))
  # Weibull
  
  if (length(CV)>1){
    ret <- cbind(1/uppr, uppr)
    colnames(ret) <- c("lower", "upper")
  } else {
    ret <- c(1/uppr, uppr)
    names(ret) <- c("lower", "upper")
  }
  ret
}