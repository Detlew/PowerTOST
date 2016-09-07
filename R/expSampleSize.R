#------------------------------------------------------------------------------
# Sample size based on 'expected' power
# taking into account the uncertainty of CV and / or uncertain theta0 (GMR)
# 
# Author: dlabes
#------------------------------------------------------------------------------
# sample size start for 'expected' power from Julious approx. for the case
# of uncertain CV
.expsampleN0 <- function(alpha=0.05, targetpower, ltheta1, ltheta2, diffm, 
                         se, dfse, steps=2, bk=2)
{
  z1 <- qnorm(1 - alpha)
  if (abs(diffm) > 0.02) {
    tinv <- qt(targetpower, dfse, z1)
  } else {
    tinv <- qt(1-(1-targetpower)/2, dfse, z1) 
    diffm <- 0
  }
  # factor 2 in Julious = bk
  n01 <- bk*(se*tinv/(ltheta1-diffm))^2
  n02 <- bk*(se*tinv/(ltheta2-diffm))^2
  n0  <- ceiling(max(n01, n02))
  #make an even multiple of step (=2 in case of 2x2 cross-over)
  n0 <- steps*trunc(n0/steps)
  #if (n0<4) n0 <- 4   # minimum sample size will be tested outside
  return(n0)
}
#------------------------------------------------------------------------------
# Sample size for a desired expected power 
# see known.designs() for covered experimental designs
# Only for log-transformed data
# leave upper BE margin (theta2) empty and the function will use 1/lower
# CV and df can be vectors -> pooled CV and df will be calculated
expsampleN.TOST <- function(alpha = 0.05, targetpower = 0.8, logscale = TRUE, 
                            theta0, theta1, theta2, CV, design = "2x2", 
                            robust = FALSE, 
                            dfCV, # to be removed in next versions 
                            prior.type = c("CV", "theta0", "both"), 
                            prior.parm = list(), 
                            method = c("exact", "approx"), 
                            print = TRUE, details) {
  # Error handling
  stopifnot(is.numeric(alpha), alpha >= 0, alpha <= 1, targetpower > 0,
            targetpower < 1, is.logical(logscale), is.numeric(CV), CV > 0, 
            is.character(design), is.logical(robust), is.logical(print))
  if (!missing(details))
    stopifnot(is.logical(details))
  # Data log-normal or normal?
  if (logscale) {
    if (missing(theta0)) theta0 <- 0.95
    if (missing(theta1)) theta1 <- 0.8
    if (missing(theta2)) theta2 <- 1/theta1
    if ((theta0 <= theta1) || (theta0 >= theta2)) {
      stop("Ratio ",theta0," not between margins ", theta1," / ", theta2,"!", 
           call. = FALSE)
    }
    ltheta1 <- log(theta1)
    ltheta2 <- log(theta2)
    diffm   <- log(theta0)
    se      <- CV2se(CV)
  } else {
    if (missing(theta0)) theta0 <- 0.05
    if (missing(theta1)) theta1 <- -0.2
    if (missing(theta1) && missing(theta2)) theta1 <- -0.2
    if (missing(theta0)) theta0 <- 0.05
    if (missing(theta2)) theta2 <- -theta1
    if ((theta0 <= theta1) || (theta0 >= theta2)) {
      stop("Diff. ", theta0, " not between margins ", theta1," / ", theta2, "!", 
           call. = FALSE)
    }
    ltheta1 <- theta1
    ltheta2 <- theta2
    diffm   <- theta0
    se      <- CV
  }
  
  # Check arguments method and uncertainty type
  method <- tolower(method)
  method <- match.arg(method)
  prior.type <- match.arg(prior.type)
  # If details is not given, set default value depending on prior.type
  # TRUE for prior.type != "CV" as sample size search may take several sec.
  if (missing(details)) {
    details <- if (prior.type == "CV") FALSE else TRUE
  }
  # Check if design is implemented and eventually get its properties
  d.no <- .design.no(design)
  if (is.na(d.no)) 
    stop("Design ", design, " not known!", call. = FALSE)
  ades <- .design.props(d.no)
  d.name <- ades$name  # nice name of design
  # Get the df for the design as an unevaluated expression
  dfe <- .design.df(ades, robust = robust)
  steps <- ades$steps  # step size
  bk <- ades$bk  # get design constant
  
  # Check how argument prior.parm was specified
  names(prior.parm) <- tolower(names(prior.parm))  # also allow "sem"
  if (length(prior.parm$sem) > 1)
    warning("SEM must be scalar, only first entry used.", call. = FALSE)
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
    df_m <- dfCV
  }
  if (sum(nms_match[3:4]) == 2) {  # m and design given
    if ((prior.parm$design == "parallel" && design != "parallel") ||
        (prior.parm$design != "parallel" && design == "parallel"))
      warning(paste0("The meaning of CV in an intra-individual design ",
                     "is not the same as in a parallel group design.", 
                     " The result may not be meaningful."), call. = FALSE)
    ds_m <- get_df_sefac(n = prior.parm$m, 
                         design = prior.parm$design, robust = robust)
    df_m <- ds_m$df  # a vector may be allowed (in contrast to exppower.TOST)
    sem_m <- ds_m$sefac * se[[1]]
  }
  if (prior.type == "CV") {
    if (nms_match[1]) {  # df given
      df_m <- prior.parm$df  # possibly overwrite df_m
    }
    if (any(is.na(df_m)))
      stop("df or combination m & design must be supplied to prior.parm!", 
           call. = FALSE)
  }
  if (sum(nms_match[1:2]) == 2) {  # df and SEM given
    df_m <- prior.parm$df  # possibly overwrite df_m
    sem_m <- prior.parm$sem[[1]]  # and sem_m
  }
  # If prior.type="CV", then CV and prior.parm$df are allowed to be vectors
  # in order to be able to calculate pooled data
  lenCV <- length(CV)
  is_pooled <- FALSE
  if (prior.type == "CV" && lenCV > 1 && length(df_m) > 1) {
    is_pooled <- TRUE
    if (length(df_m) != lenCV)
      stop("CV and df must have equal number of entries!", call. = FALSE)
    # Derive pooled SD, see e.g. Section 5.6 in Patterson and Jones
    CV <- if (logscale) CV2mse(CV) else CV^2
    CV <- CV * df_m
    df_m <- sum(df_m)
    CV <- sqrt(sum(CV) / df_m)  # pooled sd 
    se <- CV  # re-calculate standard deviation for later use
    if (logscale) CV <- se2CV(se)  # back to CV
  } else { 
    if (lenCV > 1 || length(df_m) > 1 || length(sem_m) > 1)
      warning("CV, df and SEM must be scalar, only first entry used.", 
              call. = FALSE)
  }
  # If df and CV not already scalar, do it now
  # (sem_m should and is already scalar)
  df_m <- df_m[[1]]
  CV <- CV[[1]]
  if (is.na(df_m) && is.na(sem_m))
    stop("Combination df & SEM or m & design must be supplied to prior.parm!", 
         call. = FALSE)
  # No check for !is.na(df_m) needed as it should be different by now
  if (!is.numeric(df_m) || any(df_m <= 4))
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
    if (sem_m > 0 && any(abs(semc - sem_m) / sem_m > 0.5))
      warning(paste0("Input values df and SEM do not appear to be consistent. ", 
                     "While the formal calculations may be correct, ", 
                     "the resulting power may not be meaningful."), 
              call. = FALSE)
  }
  
  # Derive minimum sample size
  n  <- 0
  df <- 0
  while (df < 1) {
    n <- n + 1
    df <- eval(dfe)
  }
  # Make it a multiple of steps
  nmin <- as.integer(steps * trunc(n / steps)) 
  nmin <- nmin + steps * (nmin < n)
  
  # Check if targetpower can be achieved at all
  # - expected power is bounded above by so-called probability of
  #   technical success (PTS)
  # - PTS calculation does not depend on n, any n does the trick here
  pts <- suppressWarnings(.exppower.TOST(alpha, ltheta1, ltheta2, diffm, se, 
                                         sqrt(bk/nmin), df, df_m, sem_m, 
                                         method, prior.type, pts = TRUE))
  if (pts <= targetpower)
    stop(paste0("Targetpower cannot be achieved because the expected power\n",
                "is bounded above by the probability of technical success (",
                formatC(pts, digits = 6, format = "f"), ").\n"), sep = "", 
         call. = FALSE)
  
  # Start printing the configuration
  if (print) {
    type.txt <- if (prior.type == "both") "CV and theta0" else prior.type
    cat("\n++++++++++++ Equivalence test - TOST ++++++++++++\n")
    blanks <- if (prior.type == "CV") 
      "       " else if (prior.type == "theta0") "     " else "  "
    cat(paste0(blanks, "Sample size est. with uncertain ", type.txt, "\n"))
    cat("-------------------------------------------------\n")
    cat("Study design: ",d.name,"\n")
    if (print) {
      if (logscale)
        cat("log-transformed data (multiplicative model)\n\n")
      else
        cat("untransformed data (additive model)\n\n")
    }
    if (details) { 
      cat("Design characteristics:\n")
      if (robust && (ades$df2 != ades$df)) {
        cat("df = ", ades$df2," (robust)", sep = "") 
      } else cat("df = ", ades$df, sep = "")
      cat(", design const. = ", bk, ", step = ", steps,"\n\n",sep = "")
    }
    cat("alpha = ", alpha,", target power = ", targetpower,"\n", sep = "")
    cat("BE margins =", theta1,"...", theta2, "\n")
    if (prior.type == "CV") {
      if (logscale) cat("Ratio = ", theta0, "\n", sep = "")
      else  cat("Diff.  = ",theta0, "\n", sep = "")
    } else {
      if (logscale) cat("Ratio = ", theta0, " with ", df_m, " df\n", sep = "")
      else  cat("Diff. = ",theta0, " with ", df_m, " df\n", sep = "")
    }
    if (is_pooled) {
      cat("CV(pooled) = ", CV, " with ", df_m," df\n", sep = "")
    } else {
      if (prior.type == "theta0") {
        cat("CV = ", CV, "\n", sep = "")
      } else {
        cat("CV = ", CV, " with ", df_m," df\n", sep = "")
      }
    }
    # Print also PTS if details = TRUE
    if (details)
      cat("\nUpper bound of expected power = ", 
          formatC(pts, digits = 6, format = "f"), "\n", sep = "")
  }
  
  # Sample size search
  # Get starting value
  if (prior.type == "CV") {
    # Use Julious approximation for this case
    n0 <- .expsampleN0(alpha, targetpower, ltheta1, ltheta2, diffm, se, df_m, 
                       steps, bk)
  } else {
    # For uncertainty wrt theta0 use starting value as in sampleN.TOST 
    # but with a shifted theta0
    # For uncertainty wrt both use starting value as in sampleN.TOST
    # but with a shifted theta0 and sigma
    sea <- if (prior.type == "theta0") se else 
      CV2se(CVCL(CV, df = df_m, side = "upper", alpha = 0.3)[[2]])
    diffma <- diffm + ifelse(diffm > 0, 1, -1) * qnorm(1 - 0.3) * sem_m
    # correct the targetpower wrt to PTS
    # doesn't has so much effect, usually 1-2 iterations
    tp <- 1 - (pts - targetpower)
    #tp <- targetpower
    n0 <- .sampleN0_3(alpha, targetpower=tp, ltheta1, ltheta2, diffma, sea, steps, bk)
  }
  if (n0 < nmin) n0 <- nmin
  if (details) {
    cat("\nSample size search (ntotal)\n")
    cat(" n   exp. power\n")
  }
  # Define sample size search function of which a root should be found
  pdiff_n <- function(n) {
    pow <- suppressWarnings(.exppower.TOST(alpha, ltheta1, ltheta2, diffm, se,
                                           sqrt(bk/n), eval(dfe), df_m, sem_m, 
                                           method, prior.type))
    if (details && !is.na(pow))
      cat(n, " ", formatC(pow, digits = 6, format = "f"), "\n")
    pow - targetpower
  }
  
  # Calculate power at n0
  #   if pow > targetpower then use c(nmin, n0) as search interval
  #   if pow <= targetpower then use c(n0, 1e+06) as search interval
  n <- n0
  pow <- .exppower.TOST(alpha, ltheta1, ltheta2, diffm, se, sqrt(bk/n), 
                        eval(dfe), df_m, sem_m, method, prior.type)
  use_nmin <- FALSE
  if (pow <= targetpower) {
    search_int <- c(n0, 1e+07)
    step.up <- TRUE
  } else {
    search_int <- c(nmin, n0)
    step.up <- FALSE
    if (search_int[1] >= search_int[2])
      use_nmin <- TRUE
  }
  # Modify step size for integer search using uniroot.step()
  # If uncertainty is higher n0 not that close, need higher step size
  # Otherwise smaller step size suffices
  step.pwr <- 2
  if (prior.type != "CV") {
    step.pwr <- if (sem_m <= 0.05) 2 else if (sem_m > 0.05 && sem_m <= 0.1) 4 else 7
    if (!step.up) # n0 with exp. power > targetpower, but is usually very close
      step.pwr <- 2
  }

  # Find sample size
  # Use uniroot.step, a variant of ssanv::uniroot.integer which restricts the
  # search to integers which are a multiple of steps
  # This should generally be faster than ssanv::uniroot.integer or stats::uniroot
  iter <- 1
  if (search_int[1]!=search_int[2]){
    n <- tryCatch({ 
      uniroot.step(f = pdiff_n, interval = search_int, step = steps, 
                   step.power = step.pwr, step.up = step.up, pos.side = TRUE)
    }, error = function(e) {
      message(e)
      return(NA)
    })
    if (all(!is.na(n))){
      pow  <- n$f.root + targetpower  # pdiff_n() substracts targetpower
      iter <- n$iter
      n <- n$root
    }
  }
  else {
    # else we have already the result n=nmin
    # print the first step n=nmin
    if (details) cat(n, " ", formatC(pow, digits = 6, format = "f"), "\n")
  }
  # Store result and print details about it
  if (!is.na(n)) {
    if (print && !details) {
      cat("\nSample size (ntotal)\n")
      cat(" n   exp. power\n")
      cat(n, " ", formatC(pow, digits = 6, format = "f"), "\n")
    }
    # print always the last step, except if n=nmin where n, power is already printed
    if(details & n!=nmin){
      cat("\n")
      cat(iter, "iterations\n")
      cat(n, " ", formatC(pow, digits = 6, format = "f"), "\n")
    }
    # give a message? in other sample size functions such a message is not given
    if (n==nmin)
      message(paste0("\nMinimum sample size already gives expected power > ",
                     "targetpower.\nSample size was set to this minimum."))

    # Print calculation method
    if(print && details) {
      if (method == "exact" || (prior.type != "CV"))
        cat("\nExact expected power calculation.\n")
      if (method == "approx" && prior.type == "CV") {
        approx <- "Approximate expected power calculation \nacc. to Julious/Owen."
        cat("\n", approx, "\n", sep = "")
      }
    }
    if(print) cat("\n")
  } else {
    cat("\nSample size search failed!\n")
    pow <- NA  # n is already NA
  }
  
  # Return results as data frame
  res <- data.frame(design = design, alpha = alpha, CV = CV, prior.df = df_m, 
                    prior.sem = sem_m, theta0 = theta0, theta1 = theta1, 
                    theta2 = theta2, n = n, power = pow, 
                    targetpower = targetpower, method = method,
                    prior.type = prior.type)
  names(res) <- c("Design", "alpha", "CV", "prior df", "prior SEM", "theta0", 
                  "theta1", "theta2", "Sample size", "Achieved power", 
                  "Target power", "Power method", "prior type")
  if (print) 
    return(invisible(res)) 
  else 
    return(res)
}
