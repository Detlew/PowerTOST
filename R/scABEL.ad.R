#-----------------------------------------------------------------------
# Iteratively adjust alpha to maintain the TIE <= nominal alpha for
# partial and full replicate design and scaled ABE via simulated
# (empirical) power
#
# Author: Helmut Schuetz
#-----------------------------------------------------------------------
scABEL.ad <-function(alpha = 0.05, theta0, theta1, theta2, CV,
                     design = c("2x3x3", "2x2x4", "2x2x3"), regulator, n,
                     alpha.pre = 0.05, imax = 100, tol, print = TRUE,
                     details = FALSE, setseed = TRUE, nsims = 1e6,
                     sdsims = FALSE, progress)
{
  ## Arguments:
  ##   alpha      Nominal alpha (in BE generally fixed to 0.05).
  ##              Lower value only if needed (e.g. to correct for
  ##              multiplicity).
  ##   theta0     If given, power is estimated for this expected ratio.
  ##              Defaults to 0.90.
  ##   theta1     Lower margin. Defaults to 0.80.
  ##   theta2     Upper margin. Defaults to 1/theta1.
  ##   CV         Intra-subject CV(s) obtained in a replicate design.
  ##              (ratio, /not/ percent).
  ##              If given as a scalar, the CV of R.
  ##              If given as a vector, CV[1] /must/ be the CV of T and
  ##              CV[2] the CV of R. Important!
  ##   design     "2x2x4", "2x2x3", "2x3x3". Defaults to "2x3x3".
  ##   regulator  "EMA", "HC", "FDA".
  ##   n          Total sample size or a vector of subjects/sequences.
  ##   nsims      Simulations for the TIE. Should not be <1E6.
  ##   imax       Max. steps in sample size search.
  ##   tol        Desired accuracy (convergence tolerance of uniroot);
  ##              defaults to 1E-6.
  ##   print      Logical (FALSE returns a list of results).
  ##   details    Logical (runtime, number of simulations).
  ##   alpha.pre  Pre-specified level.
  ##   setseed    Logical (default TRUE uses set.seed(123456)).
  ##   sdsims     Logical (default FALSE). If TRUE subject data are
  ##              simulated. Consider especially for the partial replicate,
  ##              low sample sizes and/or heteroscedasticity. Time consuming!
  ##   progress   Displays a progress bar (only if sdsims=TRUE).
  ## Returns:
  ##   alpha.adj  Iteratively adjusted alpha which does not inflate the
  ##              TIE (for given CVwR and n).
  ##   CI.adj     The adjusted confidence interval in percent, where
  ##              CI.adj = 100(1-2*alpha.adj).
  ##   TIE.unadj  The empiric Type I Error based on nominal alpha.
  ##   TIE.adj    TIE based on adjusted alpha.
  ##   rel.change Relative change in risk (%) compared to nominal alpha.
  ##   If theta0 is given:
  ##   pwr.unadj  Power for alpha (or, if given, alpha.pre).
  ##   pwr.adj    Power for adjusted alpha.
  ##   rel.loss   Relative loss in power if the sample size was planned
  ##              for alpha and will be evaluated with alpha.adj,
  ##              where rel.loss = 100(pwr.adj - pwr.unadj)/pwr.unadj
  ##   If alpha.pre is given:
  ##   Assessment of TIE; alpha.pre is justified if not > alpha.
  ################################################################
  ## Tested on Win 7 Pro SP1 64bit                              ##
  ##   R 3.4.2 64bit (2017-09-28), PowerTOST 1.4-6 (2017-08-17) ##
  ################################################################
  env <- as.character(Sys.info()[1]) # get info about the OS
  ifelse ((env == "Windows") || (env == "Darwin"), flushable <- TRUE,
          flushable <- FALSE) # supress flushing on other OS's
  # acceptance range defaults
  if (missing(theta1) && missing(theta2)) theta1 <- 0.8
  if (missing(theta2)) theta2 = 1/theta1
  # check theta0
  if (missing(theta0)) theta0 <- 0.9
  if (theta0 < theta1 || theta0 > theta2)
    stop("theta0 must be within [theta1, theta2]")
  # check regulator arg
  if(missing(regulator)) regulator <- "EMA"
  reg <- reg_check(regulator, choices=c("EMA", "HC", "FDA"))
  if (regulator %in% c("HC", "FDA") && sdsims)
    stop("Subject data simulations are not supported for regulator=\'HC\' or \'FDA\'.")
  # set iteration tolerance for uniroot().
  if (missing(tol)) tol <- 1e-6
  design <- match.arg(design)
  if (missing(CV)) stop("CV must be given!")
  CVwT <- CV[1]
  if (length(CV) == 2) CVwR <- CV[2] else CVwR <- CVwT
  if(missing(progress)) progress <- FALSE
  no <- 0 # simulation counter
  if (details) ptm <- proc.time()
  if (missing(n) || is.na(n)) {
    if (is.na(alpha.pre) || (alpha.pre != alpha)) {
      al <- alpha.pre # If pre-specified, use alpha.pre
    } else {
      al <- alpha    # If not, use alpha (commonly 0.05)
    }
    # sample size for targetpower 0.8
    if (!reg$name == "FDA") { # EMA or HC
      n <- sampleN.scABEL(alpha = al, CV = CV, theta0 = theta0,
                          theta1 = theta1, theta2 = theta2,
                          design = design, regulator = reg,
                          imax = imax, print = FALSE, details = FALSE,
                          setseed = setseed)[["Sample size"]]
    } else {                 # FDA
      n <- sampleN.RSABE(alpha = al, CV = CV, theta0 = theta0,
                         theta1 = theta1, theta2 = theta2,
                         design = design, regulator = reg,
                         imax = imax, print = FALSE, details = FALSE,
                         setseed = setseed)[["Sample size"]]
    }
    if (is.na(n))
      stop(paste0("Sample size search in sampleN.scABEL() or sampleN.RSABE() failed.\n",
                  "Restart with an explicitly high n (>1000)."))
    no <- 1e5
  }
  if (sum(n) < 6) stop("Sample size too low.")
  if (alpha.pre > alpha) {
    warning(paste0("alpha.pre > alpha does not make sense.\n",
                   "alpha.pre was set to alpha."))
    alpha.pre <- alpha
  }
  seqs <- as.numeric(substr(design, 3, 3)) # subjects / sequence
  if (length(n) == 1) n <- nvec(n, seqs)   # vectorize n
  
  # here we go!
  TIE <- pwr <- rep(NA, 2) # initialize vectors: TIE and pwr
  alpha.adj <- NA          # adjusted alpha
  # Finds adjusted alpha which gives TIE as close as possible to alpha.
  # Simulate underlying statistics (if sdsims=FALSE)
  opt1 <- function(x) {
    if (!reg$name == "FDA") { # EMA or HC
      power.scABEL(alpha = x, CV = CV, theta0 = U, n = n,
                   regulator = reg, design = design,
                   nsims = nsims, setseed = setseed) - alpha
    } else {                  # FDA
      power.RSABE(alpha = x, CV = CV, theta0 = U, n = n,
                  regulator = reg, design = design,
                  nsims = nsims, setseed = setseed) - alpha
    }
  }
  # Simulate subject data (if sdsims=TRUE)
  opt2 <- function(x) power.scABEL.sdsims(alpha = x, CV = CV, theta0 = U,
                                          n = n, regulator = reg,
                                          design = design, nsims = nsims,
                                          setseed = setseed,
                                          progress = progress) - alpha
  sig <- binom.test(x = round(alpha*nsims, 0), n = nsims,
                    alternative = "less",
                    conf.level = 1 - alpha)$conf.int[2]
  method <- "ABE"
  if (CVwR > reg$CVswitch) {
    if (regulator != "FDA") method <- "ABEL"
    if (regulator == "FDA") method <- "RSABE"
  }
  limits <- as.numeric(scABEL(CV = CVwR, regulator = reg))
  U <- limits[2] # Simulate at the upper (expanded) limit. For CVwR
  # 30% that's 1.25. Due to the symmetry simulations
  # at the lower limit (0.80) should work as well.
  if (alpha.pre != alpha) {
    al <- alpha.pre # If pre-specified, use alpha.pre.
  } else {
    al <- alpha     # If not, use alpha (commonly 0.05).
  }
  designs <- c("2x2x4", "2x2x3", "2x3x3")
  type    <- c("TRTR|RTRT", "TRT|RTR", "TRR|RTR|RRT") # clear words
  if (print) { # Show input to keep the spirits of the user high.
    if (sdsims) cat("Be patient. Simulating subject data will take a good while!\n\n")  
    cat("+++++++++++ scaled (widened) ABEL ++++++++++++\n")
    cat("         iteratively adjusted alpha\n")
    if (reg$name == "EMA") cat("   (simulations based on ANOVA evaluation)\n")
    if (reg$name == "HC") cat("(simulations based on intra-subject contrasts)\n")
    cat("----------------------------------------------\n")
    cat("Study design: ")
    cat(paste0(design, " (", type[match(design, designs)], ")\n"))
    cat("log-transformed data (multiplicative model)\n")
    cat(formatC(nsims, format = "d", big.mark = ",", decimal.mark = "."),
        "studies in each iteration simulated.\n\n")
    txt <- paste0("CVwR ", sprintf("%.4g", CVwR))
    if (length(CV) == 2) {
      txt <- paste0(txt, ", CVwT ", sprintf("%.4g", CVwT), ", ")
    } else {
      txt <- paste0(txt, ", ")
    }
    cat(paste0(txt, "n(i) ", paste0(n, collapse = "|"), " (N ", sum(n),
               ")\n"))
    txt <- paste0("Nominal alpha                 : ", signif(alpha, 5))
    if (!is.na(alpha.pre) && (alpha.pre != alpha)) {
      txt <- paste0(txt, ", pre-specified alpha ", alpha.pre, "\n")
    } else {
      txt <- paste(txt, "\n")
    }
    cat(txt)
    cat("True ratio                    :", sprintf("%.4f", theta0), "\n")
    cat(paste0("Regulatory settings           : ", reg$name, " (",
               method, ")\n"))
    
    # better theta1, theta2 as BE limits, PE constraint?
    if (CVwR <= reg$CVswitch) {
      cat("Switching CVwR                : ", reg$CVswitch, "\n",
          "BE limits                     : 0.8000 ... 1.2500\n", sep="")
    } else {
      cat(paste("Switching CVwR                :", reg$CVswitch,
                "\nRegulatory constant           :", signif(reg$r_const, 4), "\n"))
      if (!reg$name == "FDA") { # EMA or HC
        cat(sprintf("%s               : %.4f%s%.4f%s", "Expanded limits",
                    limits[1], " ... ", limits[2], "\n"))
      } else {                  # FDA
        cat(sprintf("%s             : %.4f%s%.4f%s", "Implied BE limits",
                    limits[1], " ... ", limits[2], "\n"))
      }
    }
    if (regulator != "FDA")
      cat("Upper scaling cap             : CVwR >", reg$CVcap, "\n")
    cat("PE constraints                : 0.8000 ... 1.2500\n")
    if (progress) cat("Progress of each iteration:\n")
    if (flushable) flush.console() # advance console output.
  }
  if (!sdsims) { # simulate underlying statistics
    if (!reg$name == "FDA") { # EMA or HC
      TIE[1] <- power.scABEL(alpha = al, CV = CV, theta0 = U, n = n,
                             design = design, regulator = reg,
                             nsims = nsims, setseed = setseed)
    } else {                  # FDA
      TIE[1] <- power.RSABE(alpha = al, CV = CV, theta0 = U, n = n,
                            design = design, regulator = reg,
                            nsims = nsims, setseed = setseed)
    }
    no <- no + nsims
    if (!reg$name == "FDA") { # EMA or HC
      pwr[1] <- power.scABEL(alpha = al, CV = CV, theta0 = theta0,
                             n = n, design = design, regulator = reg,
                             setseed = setseed)
    } else {                  # FDA
      pwr[1] <- power.RSABE(alpha = al, CV = CV, theta0 = theta0,
                            n = n, design = design, regulator = reg,
                            setseed = setseed)
    }
    no <- no + 1e5
  } else { # simulate subject data
    TIE[1] <- power.scABEL.sdsims(alpha = al, CV = CV, theta0 = U, n = n,
                                  design = design, regulator = reg,
                                  nsims = nsims, setseed = setseed,
                                  progress = progress)
    no <- no + nsims
    pwr[1] <- power.scABEL.sdsims(alpha = al, CV = CV, theta0 = theta0,
                                  n = n, design = design, regulator = reg,
                                  setseed = setseed, progress = progress)
    no <- no + 1e5
  }
  if (TIE[1] > alpha) { # adjust only if needed (> nominal alpha)
    if (!sdsims) { # simulate underlying statistics
      x <- uniroot(opt1, interval = c(0, alpha), tol = tol)
    } else { # simulate subject data
      x <- uniroot(opt2, interval = c(0, alpha), tol = tol)
    }
    alpha.adj <- x$root
    if (!sdsims) { # simulate underlying statistics
      if (!reg$name == "FDA") { # EMA or HC
        TIE[2] <- power.scABEL(alpha = alpha.adj, CV = CV, theta0 = U, n = n,
                               design = design, regulator = reg,
                               nsims = nsims, setseed = setseed)
        pwr[2] <- power.scABEL(alpha = alpha.adj, CV = CV, theta0 = theta0,
                               n = n, design = design, regulator = reg,
                               setseed = setseed)
      } else {                  # FDA
        TIE[2] <- power.RSABE(alpha = alpha.adj, CV = CV, theta0 = U, n = n,
                              design = design, regulator = reg,
                              nsims = nsims, setseed = setseed)
        pwr[2] <- power.RSABE(alpha = alpha.adj, CV = CV, theta0 = theta0,
                              n = n, design = design, regulator = reg,
                              setseed = setseed)
      }
    } else { # simulate subject data
      TIE[2] <- power.scABEL.sdsims(alpha = alpha.adj, CV = CV, theta0 = U,
                                    n = n, design = design, regulator = reg,
                                    nsims = nsims, setseed = setseed,
                                    progress = progress)
      pwr[2] <- power.scABEL.sdsims(alpha = alpha.adj, CV = CV,
                                    theta0 = theta0, n = n, design = design,
                                    regulator = reg, setseed = setseed,
                                    progress = progress)
    }
  }
  if (!is.na(alpha.adj)) no <- no + nsims*x$iter
  if (details) run.time <- proc.time() - ptm
  if (print) { # fetch and print results
    txt <- paste0("Empiric TIE for alpha ", sprintf("%.4f", al), "  : ",
                  sprintf("%.5f", TIE[1]))
    if (TIE[1] > alpha || alpha.pre != alpha) {
      rel.change <- 100*(TIE[1] - alpha)/alpha
      if (details) {
        txt <- paste0(txt, " (rel. change of risk: ",
                      sprintf("%+1.3g%%", rel.change), ")")
      }
    }
    if (!is.na(pwr[1])) {
      pwr.unadj <- pwr[1]
      txt <- paste0(txt, "\nPower for theta0 ", sprintf("%.4f", theta0),
                    "       : ", sprintf("%.3f", pwr.unadj))
    }
    if (TIE[1] > alpha) {
      txt <- paste0(txt, "\nIteratively adjusted alpha    : ",
                    sprintf("%.5f", alpha.adj),
                    "\nEmpiric TIE for adjusted alpha: ",
                    sprintf("%.5f", TIE[2]))
      if (!is.na(pwr[2])) {
        pwr.adj <- pwr[2]
        txt <- paste0(txt, "\nPower for theta0 ",sprintf("%.4f", theta0),
                      "       : ", sprintf("%.3f", pwr.adj))
        if (details) {
          txt <- paste0(txt, " (rel. impact: ", sprintf("%+1.3g%%",
                                                        100*(pwr[2] - pwr[1])/pwr[1]), ")")
        }
        txt <- paste(txt, "\n\n")
      } else {
        txt <- paste0(txt, "\n\n")
      }
      if (details) {
        txt <- paste0(txt, "Runtime    : ", signif(run.time[3], 3),
                      " seconds\nSimulations: ",
                      formatC(no, format = "d", big.mark = ",",
                              decimal.mark = "."), " (",
                      (no - no %% nsims - nsims) / nsims,
                      " iterations)\n\n")
      }
    } else {
      txt <- paste0(txt, "\nTIE not > nominal alpha; ")
      ifelse(alpha.pre == alpha,
             txt <- paste0(txt, "no adjustment of alpha is required.\n\n"),
             txt <- paste0(txt, "the chosen pre-specified alpha is ",
                           "justified.\n\n"))
    }
    cat(txt)
  } else { # print=FALSE
    # Prepare and return list of results.
    res <- list(regulator = reg$name, method = method, design = design,
                type = type[match(design, designs)], eval = reg$est_method,
                alpha = alpha,
                alpha.pre = ifelse(alpha.pre == alpha, NA, alpha.pre),
                CVwT = CV[1], CVwR = ifelse(length(CV)==1, CV[1], CV[2]),
                N = sum(n), theta0 = theta0, TIE.unadj = signif(TIE[1], 5),
                rel.change = ifelse(!is.na(alpha.adj),
                                    signif(100*(TIE[1] - alpha)/alpha, 5), NA),
                pwr.unadj = signif(pwr[1], 5), alpha.adj = signif(alpha.adj, 5),
                TIE.adj = signif(TIE[2], 5), pwr.adj = signif(pwr[2], 5),
                rel.loss = ifelse(!is.na(pwr[2]),
                                  signif(100*(pwr[2] - pwr[1])/pwr[1], 5), NA),
                sims = no)
    return(res)
  }
}
