# Author: H. Schuetz based on scripts by D. Labes
# -----------------------------------------------

library(PowerTOST)
# function to create sample size tables
sampsiz <- function(alpha=0.05, power, CV, GMR, theta1, design=design)
{
# 'cartesian' product
  cells <- length(CV)*length(GMR)*length(power)
  data  <- merge(CV, GMR)
  names(data) <- c("CV", "GMR")
  tbl   <- data.frame()
  pb    <- txtProgressBar(min=0, max=1, style=3)
  cell  <- 0
  for (j in seq_along(power))
  {
    data$n <- 1
    data$power <- power[j]
    for (i in seq_along(data$n)) {
      cell <- cell + 1
      setTxtProgressBar(pb, cell/cells)
      data$n[i] <- sampleN.RSABE(alpha=alpha, CV=data[i, "CV"], regulator="FDA",
                                 theta0=data$GMR[i], targetpower=power[j],
                                 theta1=theta1, design=design,
                                 print=FALSE, details=FALSE)[, "Sample size"]
    }
    data2 <- reshape(data, v.names="n", idvar=c("power", "CV"), timevar="GMR",
                     direction="wide")
    names(data2) <- gsub("n.", "R", names(data2))
    names(data2)[names(data2)=="R1"] <- "R1.0"
    names(data2)[names(data2)=="R0"] <- "R0.0"
    #cat("Power", power[j], "\n")
    #print(data2[, -2], row.names=FALSE)
    #cat("\n")
    tbl <- rbind(tbl, data2)
  }
  #shift power to first column
  tbl <- tbl[, c(2, 1, 3:ncol(tbl))]
  return(invisible(tbl))
  close(pb)
}

CVs   <- seq(0.3, 0.8, 0.05)
GMRs  <- c(seq(0.85, 1, 0.05), seq(1.05, 1.2, 0.05))
power <- c(0.8, 0.9)
# Tothfalusi L, Endrenyi L.
# "Sample Sizes for Designing Bioequivalence Studies for Highly Variable Drugs"
# J Pharm Pharmacol Sci. 2012;15(1):73-84.
txt0 <- paste('\n\nTothfalusi L, Endrenyi L.',
'\nSample Sizes for Designing Bioequivalence Studies for Highly Variable Drugs',
'\nJ Pharm Pharmacol Sci. 2012;15(1):73-84.\n\n')
note0 <- paste("\nNote: Small discrepancies since in the reference only 10,000 studies",
"\nwere simulated and sampleN.RSABE() simulates 100,000 by default.",
"\nFurthermore, sampleN.RSABE() always rounds the sample size up to",
"\nobtain balanced sequences")
# Appendix Table A3. Sample sizes for 2x3x3 designs, FDA
txt3  <- paste(txt0, "APPENDIX Table A3. FDA partial replicate (TRR|RTR|RRT)\n")
note3 <- paste(note0, "(i.e., sample size is a multiple of 3).\n\n")
tA.3  <- sampsiz(power=power, CV=CVs, GMR=GMRs, theta1=0.8, design="2x3x3")
cat(txt3);print(tA.3, row.names=FALSE);cat(note3)
# Appendix Table A4. Sample sizes for 2x2x4 designs, FDA
txt4  <- paste(txt0, "APPENDIX Table A4. FDA 4-period full replicate (TRTR|RTRT)\n")
note4 <- paste(note0, "(i.e., sample size is a multiple of 2).\n\n")
tA.4  <- sampsiz(power=power, CV=CVs, GMR=GMRs, theta1=0.8, design="2x2x4")
cat(txt4);print(tA.4, row.names=FALSE);cat(note4)
