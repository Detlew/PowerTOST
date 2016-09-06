#------------------------------------------------------------------------------
# Authors: dlabes, Benjamin Lang
#------------------------------------------------------------------------------

# Approximate "expected" power according to Julious book ----------------------
# taking into account the uncertainty of an estimated se with 
# dfse degrees of freedom
# Only for log-transformed data.
# Raw function: see the call in .exppower.TOST()
.approx.exppower.TOST <- function(alpha=0.05, ltheta1, ltheta2, diffm, sem, 
                                  dfse, df, pts = FALSE) 
{
  if (pts) {
    # Need to integrate dinvgamma() from 0 to Inf
    # cf. .exact.exppower.TOST with prior.type = "CV"
    # Result will be 1 since dinvgamma is density
    return(1)
  }
  tval <- qt(1 - alpha, df, lower.tail = TRUE)
  d1   <- sqrt((diffm-ltheta1)^2/sem^2)
  d2   <- sqrt((diffm-ltheta2)^2/sem^2)
  # in case of diffm=ltheta1 or =ltheta2 and se=0
  # d1 or d2 have then value NaN (0/0)
  d1[is.nan(d1)] <- 0
  d2[is.nan(d2)] <- 0
  pow <- pt(d1, dfse, tval) + pt(d2, dfse, tval) - 1
  # approx. may have led to negative expected power
  pow[pow < 0] <- 0
  return(pow)
}

# Exact implementation of expected power --------------------------------------
# Author(s) B. Lang & D. Labes
.exact.exppower.TOST <- function(alpha=0.05, ltheta1, ltheta2, ldiff, se,
                                 sefac_n, df_n, df_m, sem_m, prior.type,
                                 pts = FALSE) {
  # Infinite df_m should give expected power identical to conditional power   
  if (is.infinite(df_m)) {
    return(.power.TOST(alpha, ltheta1, ltheta2, ldiff, sefac_n*se, df_n))
  }
  if (prior.type == "CV") {
    # Define expected power integrand
    # For prior densitiy, see for example Bertsche et al or 
    # Held and Sabanes Bove Example 6.25
    # NB: prior density is the posterior distribution from the prior trial
    d <- function(v) {
      dinvgamma(v, shape = df_m/2, scale = df_m/2 * se^2)
    }
    # If pts = TRUE we only want to integrate the density function
    f <- function(v) if (pts) 1 else
      .power.TOST(alpha, ltheta1, ltheta2, ldiff, sefac_n*sqrt(v), df_n)
    
    # Numerical integration from 0 to Inf is likely to result in wrong result
    # because the function d (and thus f*d as well) will be zero over nearly
    # all its range. Apply Chebyshev's inequality with k=10 to avoid this
    k <- 10
    # Mean of inverse gamma with alpha=dfCV/2, beta=se^2*dfCV/2
    minvg <- (df_m/2 * se^2) / (df_m/2 - 1)
    # Variance of inverse gamma
    vinvg <- (df_m/2 * se^2)^2/(df_m/2 - 1)^2/(df_m/2 - 2)
    lwr <- if (pts) 0 else max(0, minvg - k*sqrt(vinvg))
    # Modify a bit: heavier tail to the right, use 2*k
    upr <- minvg + 2*k*sqrt(vinvg)
    pwr <- integrate(function(v) f(v) * d(v), lwr, upr, rel.tol = 1e-05, 
                     stop.on.error = FALSE)
    if (pwr$message != "OK")
      warning(pwr$message)
    return(pwr$value)
  } else if (prior.type == "theta0") {
    # Define expected power integrand
    # Use conditional density of posteriori distribution, i.e.
    # theta0 | sigma^2=sigma^2_hat
    # See e.g. Held and Sabanes Bove (6.23) and Example 6.25
    # NB: We do not use the marginal distribution/density of theta0 (i.e. the
    # non-standardized t-distribution)
    # Define lambda parameter, follows from (6.28) & (6.30) Held + Bove
    lambda <- (se / sem_m)^2  # No missing data => lambda = 1 / sefac_m^2
    # R CMD check grumbles if here d and f are used as function names:
    #* checking R code for possible problems ... NOTE
    #  .exact.exppower.TOST: multiple local function definitions for 'd' with
    #  different formal arguments
    #  .exact.exppower.TOST: multiple local function definitions for 'f' with
    # different formal arguments
    dt0 <- function(t) {
      dnorm(t, mean = ldiff, sd = se / sqrt(lambda))
      #dt_ls(t, df_m, ldiff, sem_m)  # non-standardized t-distr.
    }
    ft0 <- function(t) if (pts) 1 else 
      .power.TOST(alpha, ltheta1, ltheta2, t, sefac_n*se, df_n)
    
    s <- se / sqrt(lambda)
    k <- 5
    # If PTS = TRUE need to integrate density from ltheta1 to ltheta2
    lwr <- if (pts) ltheta1 else ldiff - k*s
    upr <- if (pts) ltheta2 else ldiff + k*s
    pwr <- integrate(function(t) ft0(t) * dt0(t), lwr, upr, rel.tol = 1e-05, 
                     stop.on.error = FALSE)
    if (pwr$message != "OK") 
      warning(pwr$message)
    return(pwr$value)
  } else if (prior.type == "both") {
    # Define expected power integrand
    # Use 2-dimensional posterior density
    # See e.g. Held and Sabanes Bove Example 6.26 for density and parameters
    lambda <- (se / sem_m)^2
    # R CMD check grumbles if here d and f are used as function names:
    #* checking R code for possible problems ... NOTE
    #  .exact.exppower.TOST: multiple local function definitions for 'd' with
    #  different formal arguments
    #  .exact.exppower.TOST: multiple local function definitions for 'f' with
    # different formal arguments
    d2 <- function(v, t) {
      dninvgamma(m = t, v = v, mu = ldiff, lambda = lambda, 
                 alpha = df_m/2, beta = df_m/2 * se^2)
    }
    f2 <- function(v, t) if (pts) 1 else
      .power.TOST(alpha, ltheta1, ltheta2, t, sefac_n*sqrt(v), df_n)
    
    k <- 10
    minvg <- (df_m/2 * se^2) / (df_m/2 - 1)
    vinvg <- (df_m/2 * se^2)^2/(df_m/2 - 1)^2/(df_m/2 - 2)
    lwr1 <- if (pts) 0 else max(0, minvg - k*sqrt(vinvg))
    upr1 <- minvg + 2*k*sqrt(vinvg)
    s <- se / sqrt(lambda)
    lwr2 <- if (pts) ltheta1 else ldiff - (k/2)*s
    upr2 <- if (pts) ltheta2 else ldiff + (k/2)*s
    # Perform 2D-integration using pracma package function
    pwr <- pracma::quad2d(function(v, t) f2(v, t) * d2(v, t), lwr1, upr1, 
                          lwr2, upr2, n = 42)
    return(pwr)
  } else {
    return(NA)
  }
}

# Working horse for exppower.TOST() and expsampleN.TOST() ---------------------
.exppower.TOST <- function(alpha=0.05, ltheta1, ltheta2, ldiff, se, sefac_n,
                           df_n, df_m, sem_m, method="exact", prior.type,
                           pts = FALSE) {
  if (method == "approx" && prior.type == "CV") {
    return(.approx.exppower.TOST(alpha, ltheta1, ltheta2, ldiff, sefac_n*se,
                                 df_m, df_n, pts))
  }
  if (method == "approx" && prior.type != "CV")
    warning(paste0("Argument method = \"approx\" only affects caluclations ",
                   "if prior.type = \"CV\"."), call. = FALSE)
  return(.exact.exppower.TOST(alpha, ltheta1, ltheta2, ldiff, se, sefac_n, 
                              df_n, df_m, sem_m, prior.type, pts))
}

# Main function for expected power --------------------------------------------
exppower.TOST <- function(alpha = 0.05, logscale = TRUE, theta0, theta1, theta2, 
                          CV, n, design = "2x2", robust = FALSE, 
                          dfCV, # to be removed in next versions,
                          prior.type = c("CV", "theta0", "both"),
                          prior.parm = list(),
                          method = c("exact", "approx")) {
  # Error handling
  stopifnot(is.numeric(alpha), alpha >= 0, alpha <= 1, is.logical(logscale),
            is.numeric(CV), CV > 0, is.numeric(n), is.character(design),
            is.logical(robust))
  if (length(CV) > 1) {
    CV <- CV[[1]]
    warning("CV has to be a scalar here. Only CV[1] used.", call. = FALSE)
  }
  
  # Data log-normal or normal?
  if (logscale) {
    if (missing(theta0)) theta0 <- 0.95
    if (missing(theta1)) theta1 <- 0.8 
    if (missing(theta2)) theta2 <- 1/theta1
    ltheta1 <- log(theta1)
    ltheta2 <- log(theta2)
    ldiff <- log(theta0)
    se <- CV2se(CV)
  } else {
    if (missing(theta0)) theta0 <- 0.05
    if (missing(theta1)) theta1 <- -0.2 
    if (missing(theta2)) theta2 <- -theta1
    ltheta1 <- theta1
    ltheta2 <- theta2
    ldiff <- theta0
    se <- CV
  }
  
  # Check method and uncertainty type
  method <- tolower(method)
  method <- match.arg(method)
  prior.type <- match.arg(prior.type)
  # Derive df and prefactor of SEM for study to be conducted
  ds_n <- get_df_sefac(n = n, design = design, robust = robust)
  if (ds_n$df < 1)
    stop("n too low. Degrees of freedom <1!", call. = FALSE)
  
  # Check how prior.parm was specified
  names(prior.parm) <- tolower(names(prior.parm))  # also allow "sem"
  if (length(prior.parm$df) > 1 || length(prior.parm$sem) > 1)
    warning("df and SEM must be scalar, only first entry used.", call. = FALSE) 
  # Check if other components have been specified
  unpar_nms <- c("df", "sem", "m", "design")  # correct possibilities
  if (sum(no_nms <- !(names(prior.parm) %in% unpar_nms)) > 0)
    warning("Unknown names in prior.parm: ", 
            paste(names(prior.parm)[no_nms], collapse = ", "), 
            call. = FALSE)
  nms_match <- unpar_nms %in% names(prior.parm)  # check which parms are given
  # Based on information in nms_match derive degrees of freedom and 
  # standard error of prior trial
  df_m <- sem_m <- NA
  if (!missing(dfCV)) {  # Temporary code for dfCV
    warning(paste0("Argument dfCV has been moved to component df of argument ", 
                   "prior.parm and will not function in the next versions."), 
            call. = FALSE)
    df_m <- dfCV[[1]]
  }
  if (sum(nms_match[3:4]) == 2) {  # m and design given
    if (prior.parm$design == "parallel" && design != "parallel")
      stop(paste0("CV in case of parallel design is total variability. This ",
                  "cannot be used to plan a future trial with ",
                  "intra-individual comparison."), call. = FALSE)
    if (prior.parm$design != "parallel" && design == "parallel")
      warning(paste0("The meaning of a CV in an intra-individual design ",
                     "is not the same as in a parallel group design.", 
                     " The result may not be meaningful."), call. = FALSE)
    ds_m <- get_df_sefac(n = prior.parm$m, 
                         design = prior.parm$design, robust = robust)
    df_m <- ds_m$df
    sem_m <- ds_m$sefac * se[[1]]
  }
  if (prior.type == "CV") {
    if (nms_match[1]) {  # df given
      df_m <- prior.parm$df[[1]]  # possibly overwrite df_m
    }
    if (is.na(df_m))
      stop("df or combination m & design must be supplied to prior.parm!", 
           call. = FALSE)
  }
  if (sum(nms_match[1:2]) == 2) {  # df and SEM given
    df_m <- prior.parm$df[[1]]  # possibly overwrite df_m
    sem_m <- prior.parm$sem[[1]]  # and sem_m
  }
  if (is.na(df_m) && is.na(sem_m))
    stop("Combination df & SEM or m & design must be supplied to prior.parm!", 
         call. = FALSE)
  # No check for !is.na(df_m) needed as it should be different by now
  if (!is.numeric(df_m) || df_m <= 4)
    stop("Degrees of freedom need to be numeric and > 4.", call. = FALSE)
  # For sem_m however this check is needed
  if (prior.type != "CV") {
    # in contrast to !is.na(sem_m) this check avoids calc. and possible warnings
    # for a case with prior.type = "CV" and prior.parm$sem specification
    if (!is.numeric(sem_m) || sem_m < 0)
      stop("SEM needs to be numeric and >= 0.", call. = FALSE)
    # Rough plausibility check for relationship between df and SEM
    # (should also give a warning if SEM > 1)
    semc <- sqrt(2 / (df_m + 2)) * CV2se(CV)
    if (sem_m > 0 && abs(semc - sem_m) / sem_m > 0.5)
      warning(paste0("Input values df and SEM do not appear to be consistent. ", 
                     "While the formal calculations may be correct, ", 
                     "the resulting power may not be meaningful."), 
              call. = FALSE)
  }
  
  # Call working horse
  .exppower.TOST(alpha = alpha, ltheta1 = ltheta1, ltheta2 = ltheta2, 
                 ldiff = ldiff, se = se, sefac_n = ds_n$sefac, df_n = ds_n$df, 
                 df_m = df_m, sem_m = sem_m, method = method, 
                 prior.type = prior.type, pts = FALSE)
}
