#-----------------------------------------------------------------------
# Sample size for partial and full replicate design and scaled ABE
# via simulated (empirical) power. Alpha is adjusted to maintain the
# empiric TIE <= nominal alpha.
#
# Author: Helmut Schuetz
#-----------------------------------------------------------------------
sampleN.scABEL.ad <- function(alpha = 0.05, targetpower = 0.8, theta0,
                              theta1, theta2, CV,
                              design = c("2x3x3", "2x2x4", "2x2x3"),
                              regulator, nstart = NA, nsims = 1e6, imax=100,
                              tol, print = TRUE, details = FALSE,
                              alpha.pre = 0.05, setseed = TRUE,
                              sdsims = FALSE, progress)
{
  env <- as.character(Sys.info()[1]) # get info about the OS
  ifelse ((env == "Windows") || (env == "Darwin"), flushable <- TRUE,
    flushable <- FALSE) # suppress flushing on other OS's
  # acceptance range defaults
  if (missing(theta1) && missing(theta2)) theta1 <- 0.8
  if (missing(theta1)) theta1 = 1/theta2
  if (missing(theta2)) theta2 = 1/theta1
  # check theta0
  if (missing(theta0)) theta0 <- 0.9
  if (theta0 < theta1 || theta0 > theta2)
    stop("theta0 must be within [theta1, theta2]")
  # check regulator arg
  if (missing(regulator)) regulator <- "EMA"
  reg <- reg_check(regulator, choices=c("EMA", "HC", "GCC"))
  if (regulator == "HC" & sdsims)
    stop("Subject data simulations are not supported for regulator=\"HC\".")
  if (length(nstart) == 2) nstart <- sum(nstart)
  design <- match.arg(design)
  if (missing(CV)) stop("CV must be given!")
  CVwT <- CV[1]
  if (length(CV) == 2) CVwR <- CV[2] else CVwR <- CVwT
  if(missing(progress)) progress <- FALSE
  if (!is.na(nstart) && nstart < 12)
    warning("Requested sample size below regulatory minimum.")
  if (!is.na(targetpower) && (targetpower < 0 || targetpower >= 1))
    stop("targetpower must be within 0 <= 1.")
  if (alpha.pre > alpha) {
    warning(paste0("alpha.pre > alpha does not make sense.",
                   "\nalpha.pre was set to alpha."))
    alpha.pre <- alpha
  }
  seqs <- as.numeric(substr(design, 3, 3))  # subjects / sequence
  sig  <- binom.test(x = round(alpha*nsims, 0), n = nsims, alternative = "less",
                     conf.level = 1 - alpha)$conf.int[2]
  method <- "ABE"
  if (CVwR > reg$CVswitch) method <- "ABEL"
  # define the data.frame of results
  res <- data.frame(design = design, regulator = reg$name,
                    method = method, eval = reg$est_method,
                    theta0 = theta0, CVwT = CVwT, CVwR = CVwR,
                    alpha = alpha, alpha.pre = alpha.pre,
                    alpha.adj = NA, TIE = NA, n = NA,
                    targetpower = targetpower, power = NA)
  names(res) <- c("Design", "Regulator", "Method", "Eval", "theta0",
                  "CVwT", "CVwR", "alpha", "alpha.pre", "alpha.adj",
                  "TIE", "Sample size", "Target power", "Achieved power")
  if (alpha.pre == alpha) res[["alpha.pre"]] <- NA
  limits <- as.numeric(scABEL(CV = CVwR, regulator = reg))
  U <- limits[2] # Simulate at the upper (expanded) limit. For CVwR
                 # 30% that's 1.25. Due to the symmetry simulations
                 # at the lower limit (0.80) would work as well.
  if (is.na(alpha.pre) || (alpha.pre != alpha)) {
    al <- alpha.pre # If pre-specified, use alpha.pre.
  } else {
    al <- alpha     # If not, use alpha (commonly 0.05).
  }
  designs <- c("2x2x4", "2x2x3", "2x3x3")
  type    <- c("4 period full replicate",
               "3 period full replicate",
               "partial replicate") # clear words
  if (print) { # Show input to keep the spirits of the user high.
    if (sdsims) cat("\nBe patient. Simulating subject data will take a good while!\n")  
    cat("\n+++++++++++ scaled (widened) ABEL ++++++++++++\n")
    cat("            Sample size estimation\n")
    cat("        for iteratively adjusted alpha\n")
    if (regulator %in% c("EMA", "GCC")) cat("   (simulations based on ANOVA evaluation)\n")
    if (regulator == "HC") cat("(simulations based on intra-subject contrasts)\n")
    cat("----------------------------------------------\n")
    cat("Study design: ")
    cat(paste0(design, " (", type[match(design, designs)], ")\n"))
    cat("log-transformed data (multiplicative model)\n")
    cat(formatC(nsims, format = "d", big.mark = ",", decimal.mark = "."),
        "studies in each iteration simulated.\n\n")
    txt <- paste0("Assumed CVwR ", sprintf("%.4g", CVwR))
    txt <- paste0(txt, ", CVwT ", sprintf("%.4g", CVwT))
    cat(txt)
    txt <- paste0("\nNominal alpha      : ", signif(alpha, 5))
    if (!is.na(alpha.pre) && (alpha.pre != alpha)) {
      txt <- paste0(txt, ", pre-specified alpha ", alpha.pre, "\n")
    } else {
      txt <- paste(txt, "\n")
    }
    cat(txt)
    cat("True ratio         :", sprintf("%.4f", theta0), "\n")
    cat("Target power       :", sprintf("%.3g", targetpower), "\n")
    cat(paste0("Regulatory settings: ", reg$name, " (", method, ")\n"))

    # better theta1, theta2 as BE limits, PE constraint?
    if (CVwR <= reg$CVswitch) {
      cat("Switching CVwR     : ", reg$CVswitch, "\n",
          "BE limits          : 0.8000 ... 1.2500\n", sep="")
    } else {
      cat(paste("Switching CVwR     :", reg$CVswitch,
                "\nRegulatory constant:", signif(reg$r_const, 5), "\n"))
      if (!regulator == "GCC") {
        cat(sprintf("%s    : %.4f%s%.4f%s", "Expanded limits",
                    limits[1], " ... ", limits[2], "\n"))
      } else {
        cat(sprintf("%s    : %.4f%s%.4f%s", "Widened limits ",
                    limits[1], " ... ", limits[2], "\n"))
      }
    }
    if (!regulator == "GCC")
      cat("Upper scaling cap  : CVwR >", reg$CVcap, "\n")
    cat("PE constraints     : 0.8000 ... 1.2500\n")
    if (progress) cat("Progress of each iteration:\n")
    if (flushable) flush.console()
  }
  if (details) ptm <- proc.time()
  no <- 0 # Simulation counter.
  if (is.na(nstart)) { # If sample size is not given, estimate one.
    unadj.n  <- sampleN.scABEL(alpha = al, targetpower = targetpower,
                               theta0 = theta0, CV = CV, design = design,
                               regulator = reg, imax=imax,
                               print = FALSE, details = FALSE, nsims = 1e5,
                               setseed = setseed)[["Sample size"]]
    no <- 1e5
  } else {             # Start with the specified sample size.
    unadj.n <- nstart
  }
  # Get results for the sample size.
  x <- scABEL.ad(alpha = alpha, theta0 = theta0, CV = CV,
                 design = design, regulator = reg, n = unadj.n,
                 nsims = nsims, imax=imax, print = FALSE, details = FALSE,
                 alpha.pre = alpha.pre, setseed = setseed,
                 sdsims = sdsims, progress = progress)
  alpha.adj <- x[["alpha.adj"]]
  if (is.na(alpha.adj)) { # No adjustment was necessary:
    if (alpha.pre != alpha) {
      alpha.adj <- alpha.pre # If pre-specified, use alpha.pre.
    } else {
      alpha.adj <- alpha     # If not, use alpha (commonly 0.05).
    }
  }
  TIE.unadj <- x[["TIE.unadj"]]
  pwr.unadj <- x[["pwr.unadj"]]
  TIE.adj   <- x[["TIE.adj"]]
  pwr.adj   <- x[["pwr.adj"]]
  no        <- no + x[["sims"]]
  if (print) {
    if (alpha.pre == alpha)
      al.txt <- "nomin. alpha:" else al.txt <- " spec. alpha:"
    if (details) {
      cat(sprintf("%s %3d, %s %.4f %s %.4f%s %.4f%s", "\nn", unadj.n,
                  al.txt, al, "(power", pwr.unadj, "), TIE:", TIE.unadj,
                  "\n"))
    }
  }
  step.1 <- FALSE # Check conditions for stopping below:
  if (TIE.unadj <= alpha && pwr.unadj > targetpower) step.1 <- TRUE
  if (!is.na(TIE.adj)) {
    if (TIE.adj <= alpha && pwr.adj > targetpower) step.1 <- TRUE
  }
  if (step.1 && is.na(TIE.adj)) { # Stop: Nothing to do.
    if (!details && print) { # only if we don't have this info already
      cat(sprintf("%s %3d, %s %.4f %s %.4f%s %.4f%s", "\nn", unadj.n,
                  al.txt, al, "(power", pwr.unadj, "), TIE:",
                  TIE.unadj, "\n"))
    }
    if (print) {
      cat("No inflation of the TIE expected; ")
      if (alpha.pre != alpha) {
        cat("the chosen pre-specified alpha is justified.\n")
      } else {
        cat("hence, no adjustment of alpha required.\n")
      }
    }
    res[["TIE"]]            <- signif(TIE.unadj, 5)
    res[["Sample size"]]    <- unadj.n
    res[["Achieved power"]] <- signif(pwr.unadj, 5)
    return(invisible(res))
  }
  if (print && details && (alpha.adj != alpha)) { # Some information.
    cat("\nSample size search and iteratively adjusting alpha")
    cat(sprintf("%s %3d, %s ", "\nn", unadj.n, "  adj. alpha:"))
    cat(sprintf("%.5f %s %.4f%s %.2f%%%s", alpha.adj, "(power", pwr.adj,
                "), rel. impact on power:",
                100*(pwr.adj - pwr.unadj)/pwr.unadj, "\n"))
    if (flushable) flush.console()
  }
  # Increase the sample size /and/ adjust alpha until achieved
  # power is at least the target power and the TIE does not
  # exceed (nominal) alpha.
  pwr <- iter <- 0
  while (pwr < targetpower) {
    if (iter == 0) { # Get the sample size for the (first!) adjusted
                     # alpha obtained from scABEL.ad(...) above.
                     # Faster than scABEL.ad() because only 1e5 sim's.
      n.new <- sampleN.scABEL(alpha = alpha.adj, CV = CV, theta0 = theta0,
                              targetpower = targetpower, design = design,
                              regulator = reg, imax=imax, print = FALSE,
                              details = FALSE, setseed = setseed)[["Sample size"]]
      no <- no + 1e5
      step.1 <- TRUE
    } else {         # In later iterations use scABEL.ad().
      if (step.1) {             # Prevents overshooting power in
        n.new <- n.new - seqs   # the 1st step, e.g, aim /lower/
        step.1 <- FALSE         # to be on the safe side!
      } else {
        n.new <- n.new + seqs   # Increase n in further steps.
      }
    }
    if (alpha.adj != alpha.pre) { # Adjust alpha (general case).
      x  <- scABEL.ad(alpha = alpha, regulator = reg, design = design,
                      CV = CV, n = n.new, theta0 = theta0, imax=imax,
                      tol = tol, print = FALSE, details = FALSE,
                      nsims = nsims, setseed = setseed,
                      sdsims = sdsims, progress = progress)
    } else {                      # Do /not/ adjust pre-specified alpha!
      x  <- scABEL.ad(regulator = reg, design = design, CV = CV,
                      n = n.new, theta0 = theta0, imax=imax, tol = tol,
                      print = FALSE, details = FALSE, nsims = nsims,
                      setseed = setseed, alpha.pre = alpha.adj,
                      sdsims = sdsims, progress = progress)
    }
    no <- no + x$sims
    if (is.na(x[["alpha.adj"]])) { # No adjustment was necessary!
      pwr       <- x[["pwr.unadj"]]
      TIE       <- x[["TIE.unadj"]]
    } else {
      pwr       <- x[["pwr.adj"]]
      alpha.adj <- x[["alpha.adj"]]
      TIE       <- x[["TIE.adj"]]
    }
    iter <- iter + 1
    # browser()
    if (pwr < targetpower && iter >= 1) { # Show intermediate steps.
      if (print && details) {
        cat(sprintf("%s %3d, %s %.5f %s %.4f%s", "n", n.new,
                    "  adj. alpha:", alpha.adj, "(power", pwr, ")\n"))        
        if (flushable) flush.console() # advance console output.
      }
    }
  }
  if (details) run.time <- proc.time() - ptm
  if (print) {
    cat(sprintf("%s %3d, %s ", "n", n.new, "  adj. alpha:"))
    cat(sprintf("%.5f %s %.4f%s %.5f%s", alpha.adj, "(power", pwr,
                "), TIE:", TIE, "\n"))
    if (details) {
      cat("Compared to nominal alpha's sample size increase of",
          sprintf("%.1f%%", 100*(n.new - unadj.n)/unadj.n),
          "(~study costs).\n\n")
    } else {
      cat("\n\n")
    }
  }
  if (print && details) {
    cat("Runtime    :", signif(run.time[3], 3), "seconds",
        "\nSimulations:", formatC(no, format = "d", big.mark = ",",
                                  decimal.mark = "."), "\n\n")
  }
  if (TIE > sig) { # Happens very, very rarely.
    warning(paste0("Algorithm failed. ",
                   "Try to restart with at least 'nstart = ",
                   n.new + seqs, "'."))
  }
  if (alpha.adj == alpha.pre) {
    res[["alpha.adj"]] <- NA
  } else {
    res[["alpha.adj"]] <- signif(alpha.adj, 5)
  }
  res[["TIE"]]            <- signif(TIE, 5)
  res[["Sample size"]]    <- n.new
  res[["Achieved power"]] <- signif(pwr, 5)
  if (print || details) {
    return(invisible(res))
  } else {
    return(res)
  }
}
