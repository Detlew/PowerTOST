## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#"
)


## ----setup--------------------------------------------------------------------
library(PowerTOST) # attach the library


## ---- echo=FALSE--------------------------------------------------------------
L <- c(sprintf("%.4f", scABEL(CV = 0.5, regulator = "EMA")[["lower"]]),
       sprintf("%.4f", scABEL(CV = 0.57382, regulator = "HC")[["lower"]]),
       "none")
U <- c(sprintf("%.4f", scABEL(CV = 0.5, regulator = "EMA")[["upper"]]),
       sprintf("%.4f", scABEL(CV = 0.57382, regulator = "HC")[["upper"]]),
       "none")
res       <- data.frame(regulator = c("EMA", "HC", "FDA"),
                        CVswitch = 0.30,  CVcap = NA, r_const = NA,
                        L = L, U = U, pe_constr = NA, method = NA,
                        stringsAsFactors = FALSE)
x         <- unlist(reg_const(regulator = "EMA"))
res[1, c(2:4, 7:8)]  <- x[c(2, 4, 3, 5:6)]
x         <- unlist(reg_const(regulator = "HC"))
res[2, c(2:4, 7:8)]  <- x[c(2, 4, 3, 6:5)]
x         <- unlist(reg_const(regulator = "FDA"))
res[3, c(2:4, 7:8)]  <- x[c(2, 4, 3, 6:5)]
res[, 3]  <- sprintf("%.5f", as.numeric(res[, 3]))
res[3, 3] <- "none"
res[, 4]  <- signif(as.numeric(res[, 4]), 5)
print(res, row.names = FALSE)


## -----------------------------------------------------------------------------
sampleN.scABEL(CV = 0.55)


## -----------------------------------------------------------------------------
CV <- signif(CVp2CV(CV = 0.55, ratio = 1.5), 4)
sampleN.scABEL(CV = CV, details = FALSE)


## -----------------------------------------------------------------------------
sampleN.scABEL.sdsims(CV = CV, details = FALSE)


## -----------------------------------------------------------------------------
CVp <- seq(0.40, 0.70, 0.05)
CV  <- signif(CVp2CV(CV = CVp, ratio = 2.5), 4)
res <- data.frame(CVp = CVp, CVwT = CV[, 1], CVwR = CV[, 2],
                  f4.key = NA, f4.ss = NA, # 4-period full replicate
                  f3.key = NA, f3.ss = NA, # 3-period full replicate
                  p3.key = NA, p3.ss = NA) # 3-period partial replicate
for (i in seq_along(CVp)) {
  res$f4.key[i] <- sampleN.scABEL(CV = CV[i, ], design = "2x2x4",
                                  print = FALSE,
                                  details = FALSE)[["Sample size"]]
  res$f4.ss[i]  <- sampleN.scABEL.sdsims(CV = CV[i, ], design = "2x2x4",
                                         print = FALSE,
                                         details = FALSE,
                                         progress = FALSE)[["Sample size"]]
  res$f3.key[i] <- sampleN.scABEL(CV = CV[i, ], design = "2x2x3",
                                  print = FALSE,
                                  details = FALSE)[["Sample size"]]
  res$f3.ss[i]  <- sampleN.scABEL.sdsims(CV = CV[i, ], design = "2x2x3",
                                         print = FALSE,
                                         details = FALSE,
                                         progress = FALSE)[["Sample size"]]
  res$p3.key[i] <- sampleN.scABEL(CV = CV[i, ], design = "2x3x3",
                                  print = FALSE,
                                  details = FALSE)[["Sample size"]]
  res$p3.ss[i]  <- sampleN.scABEL.sdsims(CV = CV[i, ], design = "2x3x3",
                                         print = FALSE,
                                         details = FALSE,
                                         progress = FALSE)[["Sample size"]]
}
print(res, row.names = FALSE)


## -----------------------------------------------------------------------------
sampleN.scABEL(CV = 0.55, regulator = "HC")


## -----------------------------------------------------------------------------
sampleN.RSABE(CV = 0.55)


## -----------------------------------------------------------------------------
CV <- signif(CVp2CV(CV = 0.125, ratio = 2.5), 4)
n  <- sampleN.NTIDFDA(CV = CV, details = FALSE)[["Sample size"]]
suppressMessages(power.NTIDFDA(CV = CV, n = n, details = TRUE))


## -----------------------------------------------------------------------------
CV <- 0.125
n  <- sampleN.NTIDFDA(CV = CV, details = FALSE)[["Sample size"]]
suppressMessages(power.NTIDFDA(CV = CV, n = n, details = TRUE))


## -----------------------------------------------------------------------------
sampleN.HVNTID(CV = 0.30, details = FALSE)


## -----------------------------------------------------------------------------
CV <- signif(CVp2CV(CV = 0.125, ratio = 2.5), 4)
sampleN.HVNTID(CV = CV, details = FALSE)


## -----------------------------------------------------------------------------
CV      <- 0.35
res     <- data.frame(n = NA, CV = CV, TIE = NA)
res$n   <- sampleN.scABEL(CV = CV, design = "2x2x4", print = FALSE,
                          details = FALSE)[["Sample size"]]
U       <- scABEL(CV = CV)[["upper"]]
res$TIE <- power.scABEL(CV = CV, n = res$n, theta0 = U, design = "2x2x4")
print(res, row.names = FALSE)


## -----------------------------------------------------------------------------
res <- data.frame(CV = sort(c(seq(0.25, 0.32, 0.01), se2CV(0.25))),
                  impl.L = NA, impl.U = NA, impl.TIE = NA,
                  des.L = NA, des.U = NA, des.TIE = NA)
for (i in 1:nrow(res)) {
  res[i, 2:3] <- scABEL(CV = res$CV[i], regulator = "FDA")
  if (CV2se(res$CV[i]) <= 0.25) {
    res[i, 5:6] <- c(0.80, 1.25)
  } else {
    res[i, 5:6] <- exp(c(-1, +1)*(log(1.25)/0.25)*CV2se(res$CV[i]))
  }
  res[i, 4] <- power.RSABE(CV = res$CV[i], theta0 = res[i, 3],
                           design = "2x2x4", n = 32, nsims = 1e6)
  res[i, 7] <- power.RSABE(CV = res$CV[i], theta0 = res[i, 5],
                           design = "2x2x4", n = 32, nsims = 1e6)
}
print(signif(res, 4), row.names = FALSE)



## -----------------------------------------------------------------------------
CV <- 0.45
n  <- sampleN.scABEL(CV = CV, design = "2x2x4", print = FALSE,
                     details = FALSE)[["Sample size"]]
scABEL.ad(CV = CV, design = "2x2x4", n = n)


## -----------------------------------------------------------------------------
CV <- 0.35
n  <- sampleN.scABEL(CV = CV, design = "2x2x4", print = FALSE,
                     details = FALSE)[["Sample size"]]
scABEL.ad(CV = CV, design = "2x2x4", n = n)


## -----------------------------------------------------------------------------
CV <- 0.35
sampleN.scABEL.ad(CV = CV, design = "2x2x4")


## ---- echo=FALSE--------------------------------------------------------------
CV  <- 0.35
des <- "2x2x4"
n   <- sampleN.scABEL(CV = CV, design = des, print = FALSE,
                      details = FALSE)[["Sample size"]]
res <- data.frame(method = c("EMA (nominal alpha)",
                             "Labes and Schütz",
                             "Molins et al."), adj = "yes",
                  alpha = 0.05, TIE = NA, power = NA,
                  stringsAsFactors = FALSE)
x <- scABEL.ad(CV = CV, n = n, design = "2x2x4", print = FALSE)
res$adj[1]   <- "no"
res$TIE[1]   <- x$TIE.unadj
res$power[1] <- x$pwr.unadj
res$alpha[2] <- x$alpha.adj
res$TIE[2]   <- x$TIE.adj
res$power[2] <- x$pwr.adj
x <- scABEL.ad(CV = 0.30, design = des, n = n, print = FALSE)
res$alpha[3] <- x$alpha.adj
res$TIE[3]   <- x$TIE.adj
res$power[3] <- power.scABEL(alpha = x$alpha.adj, CV = CV,
                             n = n, design = des)
res[, 3]     <- round(res[, 3], 5)
res[, 4]     <- round(res[, 4], 4)
res[, 5]     <- round(res[, 5], 3)
cat(paste0("CV = ", CV, ", n = ", n, ", design = \"", des, "\"\n")); print(res, row.names = FALSE)


## ---- echo=FALSE--------------------------------------------------------------
CV  <- 0.80
des <- "2x2x4"
n   <- sampleN.scABEL(CV = CV, design = des, print = FALSE,
                      details = FALSE)[["Sample size"]]
res <- data.frame(method = c("Labes and Schütz", "Molins et al."),
                  adj = "no", alpha = 0.05, TIE = NA, power = NA,
                  stringsAsFactors = FALSE)
x   <- scABEL.ad(CV = CV, n = n, design = des, print = FALSE)
res$TIE[1]   <- x$TIE.unadj
res$power[1] <- x$pwr.unadj
x   <- scABEL.ad(CV = 0.30, n = n, design = des, print = FALSE)
res$adj[2]   <- "yes"
res$alpha[2] <- x$alpha.adj
res$TIE[2]   <- x$TIE.adj
res$power[2] <- power.scABEL(alpha = x$alpha.adj, CV = CV,
                             n = n, design = "2x2x4")
res[, 3]     <- round(res[, 3], 5)
res$TIE      <- round(res$TIE, 4)
res$power    <- round(res$power, 3)
cat(paste0("CV = ", CV, ", n = ", n, ", design = \"", des, "\"\n")); print(res, row.names = FALSE)


## -----------------------------------------------------------------------------
CV <- 0.35
n  <- sampleN.scABEL(CV = CV, design = "2x2x4", print = FALSE,
                     details = FALSE)[["Sample size"]]
U  <- scABEL(CV = CV)[["upper"]]
# subject simulations and therefore, relatively slow
power.RSABE2L.sds(CV = CV, design = "2x2x4", theta0 = U,
                  n = n, SABE_test = "exact", nsims = 1e6,
                  progress = FALSE)


## -----------------------------------------------------------------------------
CV  <- c(0.30, 0.40898, 0.50, 0.57382)
res <- data.frame(CV = CV, EMA.L = NA, EMA.U = NA, EMA.cap = "",
                  HC.L = NA, HC.U = NA, HC.cap = "",
                  stringsAsFactors = FALSE)
for (i in seq_along(CV)) {
  res[i, 2:3] <- sprintf("%.4f", scABEL(CV[i], regulator = "EMA"))
  res[i, 5:6] <- sprintf("%.3f", scABEL(CV[i], regulator = "HC"))
}
res$EMA.cap[res$CV <= 0.30]   <- res$HC.cap[res$CV <= 0.30] <- "lower"
res$EMA.cap[res$CV >= 0.50]   <- "upper"
res$HC.cap[res$CV >= 0.57382] <- "upper"
print(res, row.names = FALSE)


## -----------------------------------------------------------------------------
res <- data.frame(CV = c(0.25, se2CV(0.25), 0.275, 0.3, 0.5, 1.0),
                  impl.L = NA, impl.U = NA, cap = "",
                  stringsAsFactors = FALSE)
for (i in 1:nrow(res)) {
  res[i, 2:3] <- sprintf("%.4f", scABEL(CV = res$CV[i],
                                        regulator = "FDA"))
}
res$cap[res$CV <= 0.30] <- "lower"
res$CV <- sprintf("%.3f", res$CV)
print(res, row.names = FALSE)


## -----------------------------------------------------------------------------
res <- data.frame(CV = c(0.25, se2CV(0.25), 0.275, 0.3, 0.5, 1.0),
                  des.L = NA, des.U = NA, cap = "",
                  stringsAsFactors = FALSE)
for (i in 1:nrow(res)) {
  if (CV2se(res$CV[i]) <= 0.25) {
    res[i, 2:3] <- sprintf("%.4f", c(0.80, 1.25))
  } else {
    res[i, 2:3] <- sprintf("%.4f",
                     exp(c(-1, +1)*(log(1.25)/0.25)*CV2se(res$CV[i])))
  }
}
res$cap[res$CV <= 0.30] <- "lower"
res$CV <- sprintf("%.3f", res$CV)
print(res, row.names = FALSE)


## -----------------------------------------------------------------------------
reg <- c("EMA", "HC", "FDA")
for (i in 1:3) {
  print(reg_const(regulator = reg[i]))
  cat("\n")
}

