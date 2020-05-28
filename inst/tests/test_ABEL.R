# Author: H. Schuetz based on scripts by D. Labes
# -----------------------------------------------

library(PowerTOST)
# function to create sample size tables
sampsiz <- function(alpha = 0.05, power, CV, GMR, theta1,
                    design = design, sdsims = FALSE) {
# 'cartesian' product
  cells <- length(CV)*length(GMR)*length(power)
  data  <- merge(CV, GMR)
  names(data) <- c("CV", "GMR")
  tbl   <- data.frame()
  pb    <- txtProgressBar(min = 0, max = 1, style = 3)
  cell  <- 0
  for (j in seq_along(power)) {
    data$n <- 1
    data$power <- power[j]
    for (i in seq_along(data$n)) {
      cell <- cell + 1
      setTxtProgressBar(pb, cell/cells)
      if (!sdsims) { # simulations based on key statistics (default)
        n <- sampleN.scABEL(alpha = alpha, CV = data[i, "CV"],
                            theta0 = data$GMR[i], targetpower = power[j],
                            theta1 = theta1, design = design, print = FALSE,
                            details = FALSE)[, "Sample size"]
      } else {       # subject simulations like in the paper
        n <- sampleN.scABEL.sdsims(alpha = alpha, CV = data[i, "CV"],
                            theta0 = data$GMR[i], targetpower = power[j],
                            theta1 = theta1, design = design, print = FALSE,
                            details = FALSE)[, "Sample size"]
      }
      data$n[i] <- n
    }
    data2 <- reshape(data, v.names="n", idvar=c("power", "CV"),
                     timevar="GMR", direction="wide")
    names(data2) <- gsub("n.", "R", names(data2))
    names(data2)[names(data2) == "R1"] <- "R1.0"
    names(data2)[names(data2) == "R0"] <- "R0.0"
    tbl <- rbind(tbl, data2)
  }
  # shift power to first column
  tbl <- tbl[, c(2, 1, 3:ncol(tbl))]
  return(invisible(tbl))
  close(pb)
}

CVs   <- seq(0.3, 0.8, 0.05)
GMRs  <- c(seq(0.85, 1, 0.05), seq(1.05, 1.2, 0.05))
power <- c(0.8, 0.9)
#####################################################
# Tothfalusi L, Endrenyi L.                         #
# Sample Sizes for Designing Bioequivalence Studies #
# for Highly Variable Drugs                         #
# J Pharm Pharmacol Sci. 2012;15(1):73-84.          #
#####################################################
txt0 <- paste('\n\nTothfalusi L, Endrenyi L.',
'\nSample Sizes for Designing Bioequivalence Studies for Highly Variable Drugs',
'\nJ Pharm Pharmacol Sci. 2012;15(1):73-84.\n\n')
note0 <- paste("\nNote: Small discrepancies since in the reference only",
               "10,000 studies\nwere simulated and sampleN.scABEL() simulates",
               "100,000 by default.\nFurthermore, sampleN.scABEL() rounds the",
               "sample size always up to\nobtain balanced sequences")
#########################################################
# Note: Set sdsims = TRUE for subject simulations. Much #
#       slower and practically identical for            #
#       homogenicity (CVwT = CVwR) like in this case.   #
#########################################################
# Appendix Table A1. Sample sizes for the requirements of EMA in 3-period studies
txt1  <- paste(txt0, "APPENDIX Table A1. EMA partial replicate (TRR|RTR|RRT)\n")
note1 <- paste(note0, "(i.e., the sample size is a multiple of 3).\n\n")
tA.1  <- sampsiz(power = power, CV = CVs, GMR = GMRs,
                 theta1 = 0.8, design = "2x3x3",
                 sdsims = FALSE) # Patience please if sdsims=TRUE!
cat(txt1); print(tA.1, row.names=FALSE); cat(note1)
# Appendix Table A2. Sample sizes for the requirements of EMA in 4-period studies
txt2  <- paste(txt0, "APPENDIX Table A2. EMA 4-period full replicate (TRTR|RTRT)\n")
note2 <- paste(note0, "(i.e., the sample size is a multiple of 2).\n\n")
tA.2  <- sampsiz(power = power, CV = CVs, GMR = GMRs,
                 theta1 = 0.8, design = "2x2x4",
                 sdsims = FALSE) # Patience please if sdsims=TRUE!
cat(txt2); print(tA.2, row.names=FALSE); cat(note2)
