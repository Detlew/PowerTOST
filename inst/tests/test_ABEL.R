# Author: H. Schuetz based on scripts by D. Labes
# -----------------------------------------------

library(PowerTOST)
# function to create sample size tables
sampsiz <- function(alpha=0.05, power, CV, GMR, theta1, design=design)
{
# 'cartesian' product
  data <- merge(CV, GMR)
  names(data) <- c("CV", "GMR")
  tbl  <- data.frame()
  for (j in seq_along(power))
  {
    data$n <- 1
    data$power <- power[j]
    for (i in seq_along(data$n)) {
      data$n[i] <- sampleN.scABEL(alpha=alpha, CV=data[i, "CV"], regulator="EMA",
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
} 

CVs   <- seq(0.3, 0.8, 0.05)
GMRs  <- c(seq(0.85, 1, 0.05), seq(1.05, 1.2, 0.05))
power <- c(0.8, 0.9)
# Tothfalusi L, Endrenyi L.
# "Sample Sizes for Designing Bioequivalence Studies for Highly Variable Drugs"
# J Pharm Pharmacol Sci. 2012;15(1):73-84.
txt0 <- paste('\nTothfalusi L, Endrenyi L.',
'\nSample Sizes for Designing Bioequivalence Studies for Highly Variable Drugs',
'\nJ Pharm Pharmacol Sci. 2012;15(1):73-84.\n\n')
note0 <- paste("\nNote: Small discrepancies since in the reference only 10,000 studies",
"\nwere simulated and sampleN.scABEL() simulates 100,000 by default.",
"\nFurthermore, sampleN.scABEL() always rounds the sample size up to",
"\nobtain balanced sequences")
# Appendix Table A1. Sample sizes for 2x3x3 designs, EMA
txt1  <- paste(txt0, "APPENDIX Table A1. EMA partial replicate (TRR|RTR|RRT)\n")
note1 <- paste(note0, "(i.e., the sample size is a multiple of 3).\n\n")
tA.1  <- sampsiz(power=power, CV=CVs, GMR=GMRs, theta1=0.8, design="2x3x3")
cat(txt1);print(tA.1, row.names=FALSE);cat(note1)
# Appendix Table A2. Sample sizes for 2x2x4 designs, EMA
txt2  <- paste(txt0, "APPENDIX Table A2. EMA 4-period full replicate (TRTR|RTRT)\n")
note2 <- paste(note0, "(i.e., the sample size is a multiple of 2).\n\n")
tA.2  <- sampsiz(power=power, CV=CVs, GMR=GMRs, theta1=0.8, design="2x2x4")
cat(txt2);print(tA.2, row.names=FALSE);cat(note2)
