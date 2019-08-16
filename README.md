README
================
Helmut Schütz

  - [PowerTOST](#powertost)
      - [Supported Designs](#supported-designs)
      - [Purpose](#purpose)
      - [Supported](#supported)
          - [Power and Sample Size](#power-and-sample-size)
          - [Methods](#methods)
          - [Helpers](#helpers)
      - [Defaults](#defaults)
          - [Average Bioequivalence](#average-bioequivalence)
          - [Reference-Scaled Average
            Bioequivalence](#reference-scaled-average-bioequivalence)
          - [Dose-Proportionality](#dose-proportionality)
          - [Power Analysis](#power-analysis)
      - [Examples](#examples)
          - [Parallel Design](#parallel-design)
          - [Crossover Design](#crossover-design)
          - [Replicate Designs](#replicate-designs)
          - [Dose-Proportionality](#dose-proportionality-1)
          - [Power Analysis](#power-analysis-1)
          - [Speed Comparisons](#speed-comparisons)
      - [Installation](#installation)

<!-- README.md is generated from README.Rmd. Please edit that file -->

# PowerTOST

The package contains functions to calculate power and estimate sample
size for various study designs used in (not only bio-) equivalence
studies.  
Built 2019-08-16 with R 3.6.1.

## Supported Designs

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

Although some replicate designs are more ‘popular’ than others, sample
size estimations are valid for *all* of the following designs:

| design  |  type   | sequences       |
| :-----: | :-----: | --------------- |
| `2x2x4` |  full   | TRTR / RTRT     |
| `2x2x4` |  full   | TRRT / RTTR     |
| `2x2x4` |  full   | TTRR / RRTT     |
| `2x2x3` |  full   | TRT / RTR       |
| `2x2x3` |  full   | TRT / RTR       |
| `2x3x3` | partial | TRR / RTR / RRT |

## Purpose

For various methods power can be *calculated* based on

  - nominal *α*, coefficient of variation (*CV*), deviation of test from
    reference (*θ*<sub>0</sub>), acceptance limits {*θ*<sub>1</sub>,
    *θ*<sub>2</sub>}, sample size (*n*), and design.

For all methods the sample size can be *estimated* based on

  - nominal *α*, coefficient of variation (*CV*), deviation of test from
    reference (*θ*<sub>0</sub>), acceptance limits {*θ*<sub>1</sub>,
    *θ*<sub>2</sub>}, target power, and design.

## Supported

### Power and Sample Size

Power covers balanced as well as unbalanced sequences in crossover or
replicate designs and equal/unequal group sizes in two-group parallel
designs. Sample sizes are always rounded up to achieve balanced
sequences or equal group sizes.

  - Average Bioequivalence (with arbitrary *fixed* limits).
  - Two simultaneous TOST procedures.
  - Non-inferiority *t*-test.
  - Ratio of two means with normally distributed data on the original
    scale based on Fieller’s (‘fiducial’) confidence interval.
  - ‘Expected’ power in case of uncertain (estimated) variability and/or
    uncertain *θ*<sub>0</sub>.
  - Reference-scaled bioequivalence based on simulations.
      - EMA: Average Bioequivalence with Expanding Limits (ABEL).  
      - FDA: Reference-scaled Average Bioequivalence (RSABE) for Highly
        Variable Drugs / Drug Products and Narrow Therapeutic Index
        Drugs (NTIDs).  
  - Iteratively adjust *α* to control the type I error in ABEL and
    RSABE.
  - Dose-Proportionality using the power model.

### Methods

  - Exact
      - Owen’s Q.
      - Direct integration of the bivariate non-central
        *t*-distribution.
  - Approximations
      - Non-central *t*-distribution.
      - ‘Shifted’ central *t*-distribution.

### Helpers

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
    dropouts) on power of BE decision.
  - Construct design matrices of incomplete block designs.

## Defaults

  - *α* 0.05, {*θ*<sub>1</sub>, *θ*<sub>2</sub>} (0.80, 1.25). Details
    of the sample size search (and the regulatory settings in
    reference-scaled average bioequivalence) are printed.
  - Note: In all functions values have to be given as ratios, not in
    percent.

### Average Bioequivalence

*θ*<sub>0</sub> 0.95, target power 0.80, design "2x2" (TR|RT), exact
method (Owen’s Q).

### Reference-Scaled Average Bioequivalence

*α* 0.05, point estimate constraint (0.80, 1.25), homoscedasticity
(*CV<sub>wT</sup>* = *CV<sub>wR</sup>*), scaling is based on
*CV<sub>wR</sub>*, target power 0.80, design "2x3x3" (TRR|RTR|RRT),
approximation by the non-central *t*-distribution, 100,000 simulations.

  - EMA, WHO, Health Canada, and many others: Average bioequivalence
    with expanding limits (ABEL).
  - FDA: RSABE.

#### Highly Variable Drugs / Drug Products

*θ*<sub>0</sub> 0.90.

###### EMA

Regulatory constant `0.76`, upper cap of scaling at *CV<sub>wR</sup>*
50%, evaluation by ANOVA.

###### Health Canada

Regulatory constant `0.76`, upper cap of scaling at *CV<sub>wR</sup>*
\~57.4%, evaluation by intra-subject contrasts.

###### FDA

Regulatory constant `log(1.25)/0.25`, linearized scaled ABE (Howe’s
approximation).

#### Narrow Therapeutic Index Drugs (FDA)

*θ*<sub>0</sub> 0.975, regulatory constant `log(1.11111)/0.1`, upper cap
of scaling at *CV<sub>wR</sup>* \~21.4%, design "2x2x4" (TRTR|RTRT),
linearized scaled ABE (Howe’s approximation), upper limit of the
confidence interval of *s<sub>wT</sup>*/*s<sub>wR</sup>* ≤2.5.

### Dose-Proportionality

*β*<sub>0</sub> (slope) `1+log(0.95)/log(rd)` where `rd` is the ratio of
the highest and lowest dose, target power 0.80, crossover design,
details of the sample size search suppressed.

### Power Analysis

Minimum acceptable power 0.70. *θ*<sub>0</sub>, design, conditions, and
sample size method depend on defaults of the respective approaches (ABE,
ABEL, RSABE, NTID).

## Examples

Before running the examples attach the library.

``` r
library(PowerTOST)
```

If not noted otherwise, defaults are employed.

### Parallel Design

Power for total *CV* 0.35, *θ*<sub>0</sub> 0.95, group sizes 52 and 49,
design "parallel".

``` r
power.TOST(CV = 0.35, theta0 = 0.95, n = c(52, 49), design = "parallel")
#> [1] 0.8011186
```

### Crossover Design

Sample size for assumed intra-subject *CV* 0.20.

``` r
sampleN.TOST(CV = 0.20)
#> 
#> +++++++++++ Equivalence test - TOST +++++++++++
#>             Sample size estimation
#> -----------------------------------------------
#> Study design: 2x2 crossover 
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

Sample size for equivalence of the ratio of two means with normality on
original scale based on Fieller’s (‘fiducial’) confidence interval.
*CV<sub>w</sub>* 0.20, *CV<sub>b</sub>* 0.40.  
Note the default *α* 0.025 (95% CI) of this function because it is
intended for studies with clinical endpoints.

``` r
sampleN.RatioF(CV = 0.20, CVb = 0.40)
#> 
#> +++++++++++ Equivalence test - TOST +++++++++++
#>     based on Fieller's confidence interval
#>             Sample size estimation
#> -----------------------------------------------
#> Study design: 2x2 crossover
#> Ratio of means with normality on original scale
#> alpha = 0.025, target power = 0.8
#> BE margins = 0.8 ... 1.25 
#> True ratio = 0.95,  CVw = 0.2,  CVb = 0.4
#> 
#> Sample size
#>  n     power
#> 28   0.807774
```

### Replicate Designs

#### ABE

Sample size for assumed intra-subject *CV* 0.45, *θ*<sub>0</sub> 0.90,
3-period full replicate design "2x2x3" (TRT|RTR).

``` r
sampleN.TOST(CV = 0.45, theta0 = 0.90, design = "2x2x3")
#> 
#> +++++++++++ Equivalence test - TOST +++++++++++
#>             Sample size estimation
#> -----------------------------------------------
#> Study design: 2x2x3 replicate crossover 
#> log-transformed data (multiplicative model)
#> 
#> alpha = 0.05, target power = 0.8
#> BE margins = 0.8 ... 1.25 
#> True ratio = 0.9,  CV = 0.45
#> 
#> Sample size (total)
#>  n     power
#> 124   0.800125
```

Note that the conventional model assumes homoscedasticity. For
heteroscedasticity we can ‘switch off’ all conditions of one of the
methods for reference-scaled ABE. We assume a σ<sup>2</sup> ratio of ⅔
(*i.e.*, T has a lower variability than R). Only relevant columns of the
data.frame shown.

``` r
reg <- reg_const("USER", r_const = NA, CVswitch = Inf,
                 CVcap = Inf, pe_constr = FALSE)
CV  <- round(CVp2CV(CV = 0.45, ratio = 2/3), 4)
res <- sampleN.scABEL(CV=CV, design = "2x2x3", regulator = reg,
                      details = FALSE, print = FALSE)
print(res[c(3:4, 8:9)], row.names = FALSE)
#>    CVwT   CVwR Sample size Achieved power
#>  0.3987 0.4977         126        0.80515
```

Similar sample size because the pooled *CV* is still 0.45.

#### ABEL

Sample size assuming homoscedasticity (*CV<sub>w</sub>* = 0.45).

``` r
sampleN.scABEL(CV = 0.45, details = TRUE)
#> 
#> +++++++++++ scaled (widened) ABEL +++++++++++
#>             Sample size estimation
#>    (simulation based on ANOVA evaluation)
#> ---------------------------------------------
#> Study design: 2x3x3 (TRT|RTR) 
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

#### RSABE

#### HVD(P)s

Sample size for a four-period full replicate study (TRTR|RTRT) assuming
heteroscedasticity (*CV<sub>wT</sub>* 0.40, *CV<sub>wR</sub>* 0.50).
Details of the sample size search suppressed.

``` r
sampleN.RSABE(CV = c(0.40, 0.50), design = "2x2x4", details = FALSE)
#> 
#> ++++++++ Reference scaled ABE crit. +++++++++
#>            Sample size estimation
#> ---------------------------------------------
#> Study design: 2x2x4 (TRTR|RTRT) 
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

#### NTIDs

Sample size assuming heteroscedasticity (*CV<sub>w</sub>* 0.125,
σ<sup>2</sup> ratio 2.5, *i.e.*, T has a substantially higher
variability than R). Assess additionally which one of the three
components (scaled, ABE, *s<sub>wT</sub>*/*s<sub>wR</sub>* ratio) drives
the sample size.

``` r
CV <- signif(CVp2CV(CV = 0.125, ratio = 2.5), 4)
n  <- sampleN.NTIDFDA(CV = CV)[["Sample size"]]
#> 
#> +++++++++++ FDA method for NTIDs ++++++++++++
#>            Sample size estimation
#> ---------------------------------------------
#> Study design:  2x2x4 (TRTR|RTRT) 
#> log-transformed data (multiplicative model)
#> 1e+05 studies for each step simulated.
#> 
#> alpha  = 0.05, target power = 0.8
#> CVw(T) = 0.1497, CVw(R) = 0.09433
#> True ratio     = 0.975 
#> ABE limits     = 0.8 ... 1.25 
#> Implied scABEL = 0.9056 ... 1.1043 
#> Regulatory settings: FDA 
#> - Regulatory const. = 1.053605 
#> - 'CVcap'           = 0.2142 
#> 
#> Sample size search
#>  n     power
#> 28   0.665530 
#> 30   0.701440 
#> 32   0.734240 
#> 34   0.764500 
#> 36   0.792880 
#> 38   0.816080
suppressMessages(power.NTIDFDA(CV = CV, n = n, details = TRUE))
#>        p(BE)  p(BE-sABEc)    p(BE-ABE) p(BE-sratio) 
#>      0.81608      0.93848      1.00000      0.85794
```

The *s<sub>wT</sub>*/*s<sub>wR</sub>* component shows the lowest power
and hence, drives the sample size.  
Compare that with homoscedasticity (*CV<sub>wT</sub>* =
*CV<sub>wR</sub>* = 0.125):

``` r
CV <- 0.125
n  <- sampleN.NTIDFDA(CV = CV, details = FALSE)[["Sample size"]]
#> 
#> +++++++++++ FDA method for NTIDs ++++++++++++
#>            Sample size estimation
#> ---------------------------------------------
#> Study design:  2x2x4 (TRTR|RTRT) 
#> log-transformed data (multiplicative model)
#> 1e+05 studies for each step simulated.
#> 
#> alpha  = 0.05, target power = 0.8
#> CVw(T) = 0.125, CVw(R) = 0.125
#> True ratio     = 0.975 
#> ABE limits     = 0.8 ... 1.25 
#> Regulatory settings: FDA 
#> 
#> Sample size
#>  n     power
#> 16   0.822780
suppressMessages(power.NTIDFDA(CV = CV, n = n, details = TRUE))
#>        p(BE)  p(BE-sABEc)    p(BE-ABE) p(BE-sratio) 
#>      0.82278      0.84869      1.00000      0.95128
```

Here the scaled ABE component shows the lowest power and drives the
sample size, which is much lower than in the previous example.

### Dose-Proportionality

*CV* 0.20, Doses 1, 2, and 8 units, *β*<sub>0</sub> 1, target power
0.90.

``` r
sampleN.dp(CV = 0.20, doses = c(1, 2, 8), beta0 = 1, targetpower = 0.90)
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

Note that the acceptance range of the slope depends on the ratio of the
highest and lowest doses (*i.e.*, it gets tighter for wider dose ranges
and therefore, higher sample sizes will be required).  
In an exploratory setting wider equivalence margins {*θ*<sub>1</sub>,
*θ*<sub>2</sub>} (0.50, 2.00) are recommended, which would translate in
this example to an acceptance range of `0.66667 ... 1.3333` and a sample
size of only six subjects.

### Power Analysis

Explore impact of deviations from assumptions (higher *CV*, higher
deviation of *θ*<sub>0</sub> from 1, dropouts) on power. Assumed
intra-subject *CV* 0.20, target power 0.90.

``` r
res <- pa.ABE(CV = 0.20, targetpower = 0.90)
print(res)
#> Sample size plan ABE
#>  Design alpha  CV theta0 theta1 theta2 Sample size Achieved power
#>     2x2  0.05 0.2   0.95    0.8   1.25          26      0.9176333
#> 
#> Power analysis
#> CV, theta0 and number of subjects which lead to min. acceptable power of at least 0.7:
#>  CV= 0.2729, theta0= 0.9044
#>  N = 16 (power= 0.7354)
```

If the study starts with 26 subjects (power \~0.92), the *CV* can
increase to \~0.27 **or** *θ*<sub>0</sub> decrease to \~0.90 **or** the
sample size decrease to 10 whilst power will still be ≥0.70.  
However, this is **not** a substitute for the “Sensitivity Analysis”
recommended in
[ICH-E9](https://www.ich.org/fileadmin/Public_Web_Site/ICH_Products/Guidelines/Efficacy/E9/Step4/E9_Guideline.pdf),
since in a real study a combination of all effects occurs
simultaneously. It is up to *you* to decide on reasonable combinations
and analyze their respective power.

### Speed Comparisons

#### ABE

"2x2" crossover design, intra-subject *CV* 0.17. Explore sample sizes
and achieved power for the supported methods (the 1<sup>st</sup> one is
the default).

``` r
CV   <- 0.17
expl <- data.frame(method = c("owenq", "mvt", "noncentral", "shifted"),
                   n = NA, power = NA, seconds = NA)
runs <- 20
for (i in 1:nrow(expl)) {
  start <- proc.time()[[3]]
  for (j in 1:runs) { # repeat to get better estimate of run times
    expl[i, 2:3] <- sampleN.TOST(CV = CV, method = expl$method[i],
                                 print = FALSE)[7:8]
  }
  expl[i, 4] <- (proc.time()[[3]] - start) / runs
}
print(expl, digits = 6, row.names = FALSE)
#>      method  n    power seconds
#>       owenq 14 0.805683  0.0015
#>         mvt 14 0.805690  0.1215
#>  noncentral 14 0.805683  0.0010
#>     shifted 16 0.852301  0.0005
```

The 2<sup>nd</sup> exact method is substantially slower than the
1<sup>st</sup>. The approximation based on the noncentral
*t*-distribution is slightly faster but matches the 1<sup>st</sup> exact
method closely. The approximation based on the shifted central
*t*-distribution is the fastest but might estimate a sample size higher
than necessary. Hence, it should be used only for comparative purposes.

#### ABEL

"2x2x4" full replicate design (TRTR|RTRT), homogenicity
(*CV<sub>wT</sub>* = *CV<sub>wR</sub>* 0.45). Explore sample sizes and
achieved power for the supported methods (‘key’ statistics or subject
simulations).

``` r
CV           <- c(0.45, 0.45)
design       <- "2x2x4"
expl         <- data.frame(method = c("key statistics", "subject simulations"),
                           n = NA, power = NA, seconds = NA)
start        <- proc.time()[[3]]
expl[1, 2:3] <- sampleN.scABEL(CV = CV, design = design,
                               print = FALSE, details = FALSE)[8:9]
expl[1, 4]   <- proc.time()[[3]] - start
start        <- proc.time()[[3]]
expl[2, 2:3] <- sampleN.scABEL.sdsims(CV = CV, design = design,
                                      print = FALSE, details = FALSE)[8:9]
expl[2, 4]   <- proc.time()[[3]] - start
print(expl, row.names = FALSE)
#>               method  n   power seconds
#>       key statistics 28 0.81116    0.15
#>  subject simulations 28 0.81196    2.46
```

Simulating via the ‘key’ statistics is the method of choice for speed
reasons. However, subject simulations are recommended if

  - the partial replicate design (TRR|RTR|RRT) is planned *and*
  - the special case of heterogenicity *CV<sub>wT</sub>* \>
    *CV<sub>wR</sub>* is expected.

## Installation

You can install the released version of PowerTOST from
[CRAN](https://CRAN.R-project.org) with:

``` r
package <- "PowerTOST"
inst    <- package %in% installed.packages()
if (length(package[!inst]) > 0) install.packages(package[!inst])
```

And the development version from [GitHub](https://github.com/) with:

    # install.packages("devtools")
    devtools::install_github("Detlew/PowerTOST")
