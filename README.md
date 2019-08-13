README
================

<!-- README.md is generated from README.Rmd. Please edit that file -->

## PowerTOST

The package contains functions to calculate power and estimate sample
size for various study designs used in (not only bio-) equivalence
studies.

### Supported designs

    #>    design                        name    df
    #>  parallel           2 parallel groups   n-2
    #>       2x2               2x2 crossover   n-2
    #>     2x2x2             2x2x2 crossover   n-2
    #>       3x3               3x3 crossover 2*n-4
    #>     3x6x3             3x6x3 crossover 2*n-4
    #>       4x4               4x4 crossover 3*n-6
    #>     2x2x3   2x2x3 replicate crossover 2*n-3
    #>     2x2x4   2x2x4 replicate crossover 3*n-4
    #>     2x4x4   2x4x4 replicate crossover 3*n-4
    #>     2x3x3   partial replicate (2x3x3) 2*n-3
    #>     2x4x2            Balaam's (2x4x2)   n-2
    #>    2x2x2r Liu's 2x2x2 repeated x-over 3*n-2
    #>    paired                paired means   n-1

### Purpose

For various methods power can be *calculated* based on

  - nominal *α*, coefficient of variation (*CV*), deviation of test from
    reference (*θ*<sub>0</sub>), acceptance limits {*θ*<sub>1</sub>,
    *θ*<sub>2</sub>}, sample size (*n*), and design.

For all methods the sample size can be *estimated* based on

  - nominal *α*, coefficient of variation (*CV*), deviation of test from
    reference (*θ*<sub>0</sub>), acceptance limits {*θ*<sub>1</sub>,
    *θ*<sub>2</sub>}, target power, and design.

### Supported

#### Power and sample size

  - Average bioequivalence (with arbitrary *fixed* limits).
  - Two simultaneous TOST procedures.
  - Non-inferiority *t*-test.
  - Ratio of two means with normally distributed data on the original
    scale based on Fieller’s (‘fiducial’) confidence interval.
  - ‘Expected’ power in case of uncertain (estimated) variability and/or
    uncertain *θ*<sub>0</sub>.
  - Reference-scaled bioequivalence based on simulations.  
    EMA: Average Bioequivalence with Expanding Limits (ABEL).  
    FDA: Reference-scaled Average Bioequivalence (RSABE) for Highly
    Variable Drugs / Drug Products and Narrow Therapeutic Index Drugs
    (NTIDs).
  - Iteratively adjust *α* to control the type I error in ABEL and
    RSABE.
  - Dose-proportionality using the power model.

#### Methods

  - Exact
      - Owen’s Q.
      - Direct integration of the bivariate non-central
        *t*-distribution.
  - Approximations
      - Non-central *t*-distribution.
      - ‘Shifted’ central *t*-distribution.

#### Helpers

  - Calculate *CV* from *MSE* or *SE* (and vice versa).
  - Calculate *CV* from given confidence interval.
  - Calculate *CV<sub>wR</sub>* from the upper expanded limit of an ABEL
    study.
  - Confidence interval of *CV*.
  - Pool *CV* from several studies.
  - Confidence interval for given *α*, *CV*, point estimate, sample
    size, and design.
  - Calculate *CV<sub>wT</sub>* and *CV<sub>wR</sub>* from a (pooled)
    *CV<sub>w</sub>* assuming a ratio of intra-subject variances.
  - *p*-values of the TOST procedure.
  - Analysis tool for exploration/visualization of the impact of
    expected values (*CV*, *θ*<sub>0</sub>, reduced sample size due to
    drop-outs) on power of BE decision.
  - Construct design matrices of incomplete block designs.

### Examples

1.  Power of a paralled design (group sizes 75 and 70), total *CV* 0.36,
    *θ*<sub>0</sub> 0.925. The defaults of *α* (0.05) and
    {*θ*<sub>1</sub>, *θ*<sub>2</sub>} (0.80, 1.25), and the exact
    method are automatically employed.

<!-- end list -->

``` r
PowerTOST::power.TOST(CV = 0.36, theta0 = 0.925,
                      n = c(75, 70), design = "parallel")
#> [1] 0.8009425
```

2.  Sample size for a 2×2×2 crossover study, assumed intra-subject *CV*
    0.20, and desired power 0.80. The defaults of *α* (0.05),
    *θ*<sub>0</sub> (0.95), {*θ*<sub>1</sub>, *θ*<sub>2</sub>} (0.80,
    1.25), targetpower (0.80), design (“2x2”), and the exact method are
    automatically employed.

<!-- end list -->

``` r
PowerTOST::sampleN.TOST(CV = 0.20)
#> 
#> +++++++++++ Equivalence test - TOST +++++++++++
#>             Sample size estimation
#> -----------------------------------------------
#> Study design:  2x2 crossover 
#> log-transformed data (multiplicative model)
#> 
#> alpha = 0.05, target power = 0.8
#> BE margins = 0.8 ... 1.25 
#> True ratio = 0.95,  CV = 0.2
#> 
#> Sample size (total)
#>  n     power
#> 20   0.834680
```

3.  As above but intra-subject *CV* 0.17. Explore sample sizes and
    achieved power for the supported methods (the 1<sup>st</sup> one is
    the default).

<!-- end list -->

``` r
expl <- data.frame(method = c("owenq", "mvt", "noncentral", "shifted"),
                   n = NA, power = NA, seconds = NA)
for (i in 1:nrow(expl)) {
  start <- proc.time()[[3]]
  expl[i, 2:3] <- PowerTOST::sampleN.TOST(CV = 0.17,
                                          method = expl$method[i],
                                          print = FALSE)[7:8]
  expl[i, 4] <- proc.time()[[3]] - start
}
print(expl, digits = 6, row.names = FALSE)
#>      method  n    power seconds
#>       owenq 14 0.805683    0.00
#>         mvt 14 0.805690    0.13
#>  noncentral 14 0.805683    0.00
#>     shifted 16 0.852301    0.00
```

  - The 2<sup>nd</sup> exact method is substantially slower than the
    1<sup>st</sup>. The approximation based on the noncentral
    *t*-distribution is – sometimes – slightly slower but matches the
    1<sup>st</sup> exact method closely. However, it is up to 100times
    faster in replicate designs and therefore, the default in methods
    for reference-scaling. The approximation based on the shifted
    central *t*-distribution might estimate a sample size higher than
    necessary. Hence, it should be used only for comparative purposes.

<!-- end list -->

4.  Sample size for a study intended for the EMA’s ABEL assuming
    homoscedasticity (*CV<sub>wT</sub>* = *CV<sub>wR</sub>* = 0.45). The
    defaults of *α* (0.05), *θ*<sub>0</sub> (0.90 for HVD(P)s),
    targetpower (0.80), the partial replicate design (TRR|RTR|RRT), and
    the approximation based on the non-central *t*-distribution method
    are automatically employed. The expanded limits \[*L*, *U*\] are
    based on *CV<sub>wR</sub>*. 100,000 studies are simulated by
    default; details of the regulatory settings and the sample size
    search shown.

<!-- end list -->

``` r
PowerTOST::sampleN.scABEL(CV = 0.45, details = TRUE)
#> 
#> +++++++++++ scaled (widened) ABEL +++++++++++
#>             Sample size estimation
#>    (simulation based on ANOVA evaluation)
#> ---------------------------------------------
#> Study design:  2x3x3 (partial replicate) 
#> log-transformed data (multiplicative model)
#> 1e+05 studies for each step simulated.
#> 
#> alpha  = 0.05, target power = 0.8
#> CVw(T) = 0.45; CVw(R) = 0.45
#> True ratio = 0.9
#> ABE limits / PE constraint = 0.8 ... 1.25 
#> EMA regulatory settings
#> - CVswitch            = 0.3 
#> - cap on scABEL if CVw(R) > 0.5
#> - regulatory constant = 0.76 
#> - pe constraint applied
#> 
#> 
#> Sample size search
#>  n     power
#> 36   0.7755 
#> 39   0.8059
```

5.  Sample size for a four-period full replicate study (TRTR|RTRT)
    intended for the FDA’s RSABE assuming heteroscedasticity
    (*CV<sub>wT</sub>* 0.40, *CV<sub>wR</sub>* 0.50). The defaults of
    *α* (0.05), *θ*<sub>0</sub> (0.90), targetpower (0.80), and the
    approximation based on the non-central *t*-distribution method are
    automatically employed. Scaling is based on *CV<sub>wR</sub>*.
    100,000 studies are simulated by default; details of the sample size
    search suppressed.

<!-- end list -->

``` r
PowerTOST::sampleN.RSABE(CV = c(0.40, 0.50), design = "2x2x4",
                         details = FALSE)
#> 
#> ++++++++ Reference scaled ABE crit. +++++++++
#>            Sample size estimation
#> ---------------------------------------------
#> Study design:  2x2x4 (full replicate) 
#> log-transformed data (multiplicative model)
#> 1e+05 studies for each step simulated.
#> 
#> alpha  = 0.05, target power = 0.8
#> CVw(T) = 0.4; CVw(R) = 0.5
#> True ratio = 0.9
#> ABE limits / PE constraints = 0.8 ... 1.25 
#> Regulatory settings: FDA 
#> 
#> Sample size
#>  n    power
#> 20   0.81509
```

6.  Sample size for a study intended for the FDA’s RSABE for NTIDs
    assuming heteroscedasticity (*CV<sub>wT</sub>* 0.15,
    *CV<sub>wR</sub>* 0.10). The defaults of *α* (0.05), *θ*<sub>0</sub>
    (0.975 for NTIDs), targetpower (0.80), TRTR|RTRT design, and the
    approximation based on the non-central *t*-distribution method are
    automatically employed. Scaling is based on *CV<sub>wR</sub>*.
    100,000 studies are simulated by default. Assess which of the three
    conditions (scaled, ABE, *s<sub>wT</sub>*/*s<sub>wR</sub>* ratio)
    drives the sample size.

<!-- end list -->

``` r
CV <- c(0.15, 0.10)
n  <- PowerTOST::sampleN.NTIDFDA(CV = CV, details = FALSE)[["Sample size"]]
#> 
#> +++++++++++ FDA method for NTIDs ++++++++++++
#>            Sample size estimation
#> ---------------------------------------------
#> Study design:  2x2x4 
#> log-transformed data (multiplicative model)
#> 1e+05 studies for each step simulated.
#> 
#> alpha  = 0.05, target power = 0.8
#> CVw(T) = 0.15, CVw(R) = 0.1
#> True ratio     = 0.975 
#> ABE limits     = 0.8 ... 1.25 
#> Regulatory settings: FDA 
#> 
#> Sample size
#>  n     power
#> 32   0.817560
suppressMessages(PowerTOST::power.NTIDFDA(CV = CV, n = n, details = TRUE))
#>        p(BE)  p(BE-sABEc)    p(BE-ABE) p(BE-sratio) 
#>      0.81756      0.92270      1.00000      0.87137
```

  - Here the *s<sub>wT</sub>*/*s<sub>wR</sub>* component has the lowest
    power and hence, drives the sample size.  
    Compare that with homoscedasticity (*CV<sub>wT</sub>* =
    *CV<sub>wR</sub>* = 0.10):

<!-- end list -->

``` r
CV <- 0.10
n  <- PowerTOST::sampleN.NTIDFDA(CV = CV, details = FALSE)[["Sample size"]]
#> 
#> +++++++++++ FDA method for NTIDs ++++++++++++
#>            Sample size estimation
#> ---------------------------------------------
#> Study design:  2x2x4 
#> log-transformed data (multiplicative model)
#> 1e+05 studies for each step simulated.
#> 
#> alpha  = 0.05, target power = 0.8
#> CVw(T) = 0.1, CVw(R) = 0.1
#> True ratio     = 0.975 
#> ABE limits     = 0.8 ... 1.25 
#> Regulatory settings: FDA 
#> 
#> Sample size
#>  n     power
#> 18   0.841790
suppressMessages(PowerTOST::power.NTIDFDA(CV = CV, n = n, details = TRUE))
#>        p(BE)  p(BE-sABEc)    p(BE-ABE) p(BE-sratio) 
#>      0.84179      0.85628      1.00000      0.97210
```

  - Now the scaled ABE component has the lowest power and drives the –
    much lower – sample size.

<!-- end list -->

7.  Sample size for a dose-proportionality study assessed by the power
    model, doses 1, 2, 8 units, assumed *CV* 0.20, slope 1, target power
    0.90. The defaults of *α* (0.05), {*θ*<sub>1</sub>, *θ*<sub>2</sub>}
    (0.80, 1.25), and a crossover design are automatically employed.

<!-- end list -->

``` r
PowerTOST::sampleN.dp(CV = 0.20, doses = c(1, 2, 8),
                      beta0 = 1, targetpower = 0.90)
#> 
#> ++++ Dose proportionality study, power model ++++
#>             Sample size estimation
#> -------------------------------------------------
#> Study design: crossover (3x3 Latin square) 
#> alpha = 0.05, target power = 0.9
#> Equivalence margins of R(dnm) = 0.8 ... 1.25 
#> Doses = 1 2 8 
#> True slope = 1, CV = 0.2
#> Slope acceptance range = 0.89269 ... 1.1073 
#> 
#> Sample size (total)
#>  n     power
#> 18   0.915574
```

### Installation

You can install the released version of PowerTOST from
[CRAN](https://CRAN.R-project.org) with:

``` r
package <- "PowerTOST"
inst    <- package %in% installed.packages()
if (length(package[!inst]) > 0) install.packages(package[!inst])
```

And the development version from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("Detlew/PowerTOST")
```
