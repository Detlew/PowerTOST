exppower.noninf <- function(alpha = 0.025, logscale = TRUE, theta0, margin, 
                            CV, n, design = "2x2", robust = FALSE,
                            prior.type = c("CV", "theta0", "both"),
                            prior.parm = list(),
                            method=c("exact", "approx")) {
  # Error handling
  stopifnot(is.numeric(alpha), alpha >= 0, alpha <= 1, is.logical(logscale),
            is.numeric(CV), CV > 0, is.numeric(n), is.character(design),
            is.logical(robust))
  if (length(CV) > 1) {
    CV <- CV[[1]]
    warning("CV has to be a scalar here. Only CV[1] used.", call. = FALSE)
  }
  
  # Data log-normal or normal?
  theta2 <- Inf  # for non-inferiority
  if (logscale) {
    if (missing(theta0)) theta0 <- 0.95
    if (theta0 <= 0) 
      stop("theta0 must be > 0.", call. = FALSE)
    # We use the same notation as for TOST (for ease of recycling code)
    if (missing(margin)) {
      theta1 <- 0.8
    } else {
      theta1 <- margin
    }
    if (theta1 < 0) stop("margin must be >= 0!", call. = FALSE)
    if ((theta0 <= theta1) && (theta1 < 1))  # non-inferiority error
      stop("Ratio ",theta0," must be above margin ", theta1, "!", call. = FALSE)
    if ((theta0 < theta1) && (theta1 > 1)) {
      # non-superiority
      # reduce this case to non-inferiority case
      theta1 <- 1/theta1
      theta0 <- 1/theta0
    }
    if ((theta0 >= theta1) && (theta1 > 1))  # non-superiority error
      stop("Ratio ",theta0," must be below margin ", theta1, "!", call. = FALSE)
    ltheta1 <- log(theta1)
    ltheta2 <- log(theta2)
    ldiff <- log(theta0)
    se <- CV2se(CV)
  } else {
    if (missing(theta0)) theta0 <- -0.05
    if (missing(margin)) {
      theta1 <- -0.2 
    } else {
      theta1 <- margin
    }
    if ((theta0 <= theta1) && (theta1 < 0))  # non-inf. error
      stop("Diff. ", theta0, " must be above margin ", theta1,"!", 
           call. = FALSE)
    if ((theta0 < theta1) && (theta1 > 0)) {
      # non-superiority
      # reduce this case to non-inf
      theta1 <- -theta1
      theta0 <- -theta0
    }
    if ((theta0 >= theta1) && (theta1 > 0)) {  # non-sup. error
      stop("Diff. ", theta0, " must be below margin ", theta1,"!", 
           call. = FALSE)
    }
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
  # Note: "exact conditional power for non-inf." can be obtained via 
  #       .exppower.TOST using ltheta2=Inf and cp_method="nct". 
  #       Using cp_method="nct" results in the same as power.noninf()
  #       but is faster than via cp_method="exact".
  .exppower.TOST(alpha = alpha, ltheta1 = ltheta1, ltheta2 = ltheta2, 
                 ldiff = ldiff, se = se, sefac_n = ds_n$sefac, df_n = ds_n$df, 
                 df_m = df_m, sem_m = sem_m, method = method, 
                 prior.type = prior.type, pts = FALSE, cp_method = "nct")
  
}
