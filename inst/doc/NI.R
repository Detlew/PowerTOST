## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#"
)


## ----setup--------------------------------------------------------------------
library(PowerTOST) # attach the library


## -----------------------------------------------------------------------------
sampleN.noninf(CV = 0.25)


## -----------------------------------------------------------------------------
sampleN.noninf(CV = 0.25, details = FALSE, print = FALSE)[["Sample size"]]


## -----------------------------------------------------------------------------
power.noninf(CV = 0.25, n = 35)


## -----------------------------------------------------------------------------
sampleN.noninf(CV = 0.25, margin = 1.25, theta0 = 1/0.95)


## -----------------------------------------------------------------------------
res <- data.frame(design = "2x2x4", metric = c("Cmin", "Cmax"),
                  margin = c(0.80, 1.25), CV = c(0.35, 0.20),
                  theta0 = c(0.95, 1.05), n = NA, power = NA)
for (i in 1:2) {
  res[i, 6:7] <- sampleN.noninf(design = res$design[i],
                                margin = res$margin[i],
                                theta0 = res$theta0[i],
                                CV = res$CV[i],
                                details = FALSE,
                                print = FALSE)[6:7]
}
print(res, row.names = FALSE)


## -----------------------------------------------------------------------------
power.noninf(design = "2x2x4", margin = 1.25, CV = 0.20,
             theta0 = 1.05, n = 32)


## -----------------------------------------------------------------------------
power.noninf(design = "2x2x4", margin = 1.25, CV = 0.25,
             theta0 = 1.10, n = 32) # higher CV, worse theta0


## -----------------------------------------------------------------------------
res <- data.frame(design = "2x2x4", indended = c("ABEL", "ABE"),
                  metric = c("Cmin", "Cmax"), CV = c(0.35, 0.20),
                  theta0 = c(0.90, 1.05), n = NA, power = NA,
                  stringsAsFactors = FALSE)
res[1, 6:7] <- sampleN.scABEL(CV = res$CV[1], theta0 = res$theta0[1],
                              design = res$design[1], print = FALSE,
                              details = FALSE)[8:9]
res[2, 6:7] <- sampleN.TOST(CV = res$CV[2], theta0 = res$theta0[2],
                            design = res$design[2], print = FALSE,
                            details = FALSE)[7:8]
print(res, row.names = FALSE)


## -----------------------------------------------------------------------------
n <- sampleN.scABEL(CV = 0.35, theta0 = 0.90, design = "2x2x4",
                    print = FALSE, details = FALSE)[["Sample size"]]
# CV and theta0 of both metrics worse than assumed
res <- data.frame(design = "2x2x4", indended = c("ABEL", "ABE"),
                  metric = c("Cmin", "Cmax"), CV = c(0.50, 0.25),
                  theta0 = c(0.88, 1.12), n = n, power = NA,
                  stringsAsFactors = FALSE)
res[1, 7] <- power.scABEL(CV = res$CV[1], theta0 = res$theta0[1],
                            design = res$design[1], n = n)
res[2, 7] <- power.TOST(CV = res$CV[2], theta0 = res$theta0[2],
                          design = res$design[2], n = n)
print(res, row.names = FALSE)
