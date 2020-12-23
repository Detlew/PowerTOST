# ----------------------------------------------------------------------------
# Sample size estimation based on replicate design subject data simulations
# and evaluation via EMA ABEL method (ANOVA & average BE with expanding limits)
#
# Author H.Schütz based on code by D.Labes
# ----------------------------------------------------------------------------
sampleN.scABEL.sdsims <- function(alpha=0.05, targetpower=0.8, theta0, theta1,
                                  theta2, CV,
                                  design=c("2x3x3", "2x2x4", "2x2x3"),
                                  regulator, nsims=1e5, nstart, imax=100,
                                  print=TRUE, details=TRUE, setseed=TRUE,
                                  progress) {
  if (missing(theta1) & missing(theta2)) theta1 <- 0.8
  if (missing(theta2)) theta2 <- 1/theta1
  if (missing(theta1)) theta1 <- 1/theta2
  if (missing(theta0)) theta0 <- 0.9
  # theta0 in range
  if ( (theta0 < theta1) | abs(theta0-theta1) <1e-8  
       | (theta0 > theta2) | abs(theta0-theta2) <1e-8 ){
    stop("True ratio ", theta0," not within margins ", theta1," ... ",
         theta2,"!", call.=FALSE)
  }

  if (missing(CV)) stop("CV must be given!", call.=FALSE)
  if (missing(design)) design <- "2x3x3"
  if (missing(regulator)) regulator <- "EMA"
  reg <- reg_check(regulator)
  if (reg$est_method=="ISC")
    stop("ISC evaluation not allowed in this function.", call.=FALSE)
  CVwT <- CV[1]
  if (length(CV) == 2) CVwR <- CV[2] else CVwR <- CVwT
  # start with key stat sims
  n <- sampleN.scABEL(alpha=alpha,targetpower=targetpower, theta0=theta0,
                      theta1=theta1, theta2=theta2, CV=CV, design=design,
                      regulator=regulator, nsims=nsims, nstart=nstart,
                      imax=imax, print=FALSE, details=FALSE,
                      setseed=setseed)[["Sample size"]]
  if (missing(progress)) {
    progress <- FALSE
    if (nsims >= 5e5 | (nsims >= 1e5 & n > 72)) progress <- TRUE
  }
  if (print){
    designs <- c("2x2x4", "2x2x3", "2x3x3")
    type    <- c("4 period full replicate",
                 "3 period full replicate",
                 "partial replicate") # clear words
    cat("Be patient. Simulating subject data may take a while!\n")
    cat("\n+++++++++++ scaled (widened) ABEL +++++++++++\n")
    cat("            Sample size estimation\n")
    cat("   (simulation based on ANOVA evaluation)\n")
    cat("---------------------------------------------\n")
    cat("Study design: ")
    cat(paste0(design, " (", type[match(design, designs)], ")\n"))
    cat("log-transformed data (multiplicative model)\n")
    cat(nsims, "studies for each step simulated.\n\n")
    cat("alpha  = ", alpha,", target power = ", targetpower,"\n", sep="")
    cat("CVw(T) = ", CVwT,"; CVw(R) = ", CVwR,"\n", sep="")
    cat("True ratio = ",theta0,"\n", sep="")
    cat("ABE limits / PE constraint =", theta1, "...", theta2, "\n")
    if (reg$name == "GCC" & CVwR > 0.3)
      cat("Widened limits =", 0.75, "...", 1/0.75, "\n")
    if (details | reg$name == "USER") { 
      print(reg)
      cat("\n")
    } else {
      cat("Regulatory settings:", reg$name, "\n")
    }
  }
  pd <- max(4, round(log10(nsims), 0) - 1) # digits for power
  if (details) {
    cat("Sample size search\n")
    cat(" n    power\n")
  }
  repeat {
    pwr <- power.scABEL.sdsims(alpha=alpha, theta1=theta1, theta2=theta2,
                               theta0=theta0, CV=CV, n=n, design=design,
                               regulator=regulator, nsims=nsims, details=FALSE,
                               setseed=setseed, progress=progress)
    if (details) {
      cat(n, " ", formatC(pwr, digits = pd, format = "f"), "\n")
    }
    if (pwr >= targetpower) {
      break
    } else {
      if (design == "2x3x3") {
        n <- n + 3
      } else {
        n <- n + 2
      }
    }
  }
  if (print && !details) {
    cat("\nSample size\n")
    cat(" n   power\n")
    cat(n," ", formatC(pwr, digits = pd, format = "f"),"\n")
    if (is.na(n)) cat("Sample size search failed!\n")
  }
  if (print) cat("\n")
  res <- data.frame(design = design, alpha = alpha, CVwT = CVwT,
                    CVwR = CVwR, theta0 = theta0, theta1 = theta1,
                    theta2 = theta2, n = n, power = pwr,
                    targetpower = targetpower, nlast = n)
  names(res) <-c("Design", "alpha", "CVwT", "CVwR", "theta0",
                 "theta1", "theta2", "Sample size", "Achieved power",
                 "Target power", "nlast")
  if (print | details) return(invisible(res)) else return(res)
}

# alias 
sampleN.scABEL.sds <- sampleN.scABEL.sdsims