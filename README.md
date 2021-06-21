PowerTOST
================

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
  - [Session Information](#session-information)

<!-- README.md is generated from README.Rmd. Please edit that file -->

[![License: GPL
v2](https://img.shields.io/badge/License-GPLv2-blue.svg)](https://www.gnu.org/licenses/gpl-2.0)
[![License: GPL
v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
![repo](https://img.shields.io/badge/repo%20since-Jun%202016-brightgreen)
![active](https://www.repostatus.org/badges/latest/active.svg) ![repo
size](https://img.shields.io/github/repo-size/Detlew/PowerTOST?color=yellow)
![code
size](https://img.shields.io/github/languages/code-size/Detlew/PowerTOST?color=green)
![first](https://img.shields.io/badge/CRAN%20since-May%202010-brightgreen)
![on CRAN](https://www.r-pkg.org/badges/version-ago/PowerTOST) [![cran
checks](https://cranchecks.info/badges/worst/PowerTOST)](https://cran.r-project.org/web/checks/check_results_PowerTOST.html)
[![CRAN RStudio mirror
downloads](https://cranlogs.r-pkg.org/badges/grand-total/PowerTOST?color=blue)](https://r-pkg.org/pkg/PowerTOST)
[![CRAN RStudio mirror
downloads](https://cranlogs.r-pkg.org/badges/last-month/PowerTOST?color=green)](https://r-pkg.org/pkg/PowerTOST)

The package contains functions to calculate power and estimate sample
size for various study designs used in (not only bio-) equivalence
studies.  
Version 1.5.3.9000 built 2021-05-26 with R 4.1.0 (development version
not on CRAN).

## Supported Designs

    #R>    design                        name    df
    #R>  parallel           2 parallel groups   n-2
    #R>       2x2               2x2 crossover   n-2
    #R>     2x2x2             2x2x2 crossover   n-2
    #R>       3x3               3x3 crossover 2*n-4
    #R>     3x6x3             3x6x3 crossover 2*n-4
    #R>       4x4               4x4 crossover 3*n-6
    #R>     2x2x3   2x2x3 replicate crossover 2*n-3
    #R>     2x2x4   2x2x4 replicate crossover 3*n-4
    #R>     2x4x4   2x4x4 replicate crossover 3*n-4
    #R>     2x3x3   partial replicate (2x3x3) 2*n-3
    #R>     2x4x2            Balaam's (2x4x2)   n-2
    #R>    2x2x2r Liu's 2x2x2 repeated x-over 3*n-2
    #R>    paired                paired means   n-1

Codes of designs follow this pattern: `treatments x sequences x
periods`.

Although some replicate designs are more ‘popular’ than others, sample
size estimations are valid for *all* of the following designs:

| <small>design</small> | <small>type</small> | <small>sequences</small> |                             | <small>periods</small> |
| :-------------------: | :-----------------: | :----------------------: | --------------------------- | :--------------------: |
|        `2x2x4`        |     <small>full     |     <small>2</small>     | `TRTR\\|RTRT`               |    <small>4</small>    |
|        `2x2x4`        |     <small>full     |     <small>2</small>     | `TRRT\\|RTTR`               |    <small>4</small>    |
|        `2x2x4`        |     <small>full     |     <small>2</small>     | `TTRR\\|RRTT`               |    <small>4</small>    |
|        `2x4x4`        |     <small>full     |     <small>4</small>     | `TRTR\\|RTRT\\|TRRT\\|RTTR` |    <small>4</small>    |
|        `2x4x4`        |     <small>full     |     <small>4</small>     | `TRRT\\|RTTR\\|TTRR\\|RRTT` |    <small>4</small>    |
|        `2x2x3`        |     <small>full     |     <small>2</small>     | `TRT\\|RTR`                 |    <small>3</small>    |
|        `2x2x3`        |     <small>full     |     <small>2</small>     | `TRR\\|RTT`                 |    <small>3</small>    |
|        `2x4x2`        |     <small>full     |     <small>4</small>     | `TR\\|RT\\|TT\\|RR`         |    <small>2</small>    |
|        `2x3x3`        |   <small>partial    |     <small>3</small>     | `TRR\\|RTR\\|RRT`           |    <small>3</small>    |
|        `2x2x3`        |   <small>partial    |     <small>2</small>     | `TRR\\|RTR`                 |    <small>3</small>    |

Balaam’s design `TR|RT|TT|RR` should be avoided due to its poor power
characteristics. The three period partial replicate design with two
sequences `TRR|RTR` (<span title="also known as">a.k.a.</span>
extra-reference design) should be avoided because it is biased in the
presence of period effects.

<small>[TOC ↩](#powertost)</small>

## Purpose

For various methods power can be *calculated* based on

  - nominal *α*, coefficient of variation (*CV*), deviation of test from
    reference (*θ*<sub>0</sub>), acceptance limits {*θ*<sub>1</sub>,
    *θ*<sub>2</sub>}, sample size (*n*), and design.

For all methods the sample size can be *estimated* based on

  - nominal *α*, coefficient of variation (*CV*), deviation of test from
    reference (*θ*<sub>0</sub>), acceptance limits {*θ*<sub>1</sub>,
    *θ*<sub>2</sub>}, target (*i.e.*, desired) power, and design.

<small>[TOC ↩](#powertost)</small>

## Supported

### Power and Sample Size

Power covers balanced as well as unbalanced sequences in crossover or
replicate designs and equal/unequal group sizes in two-group parallel
designs. Sample sizes are always rounded up to achieve balanced
sequences or equal group sizes.

  - Average Bioequivalence (with arbitrary *fixed* limits).
  - Two simultaneous <span title="Two One-Sided Tests">TOST</span>
    procedures.
  - Non-inferiority *t*-test.
  - Ratio of two means with normally distributed data on the original
    scale based on Fieller’s (‘fiducial’) confidence interval.
  - ‘Expected’ power in case of uncertain (estimated) variability and/or
    uncertain *θ*<sub>0</sub>.
  - Reference-scaled bioequivalence based on simulations.
      - <span title="European Medicines Agency">EMA</span>,
        <span title="World Health Organization">WHO</span> and many
        others: Average Bioequivalence with Expanding Limits (ABEL).
      - <span title="Gulf Co-operation Council">GCC</span>: Average
        Bioequivalence with *fixed* widened limits of 75.00–133.33% if
        *CV*<sub>wR</sub> \>30%.  
      - U.S. <span title="Food and Drug Administration">FDA</span>,
        China <span title="Centre for Drug Evaluation">CDE</span>:
        Reference-scaled Average Bioequivalence (RSABE) for Highly
        Variable Drugs / Drug Products and Narrow Therapeutic Index
        Drugs (NTIDs). In China the former is required and the latter
        currently under discussion.
  - Iteratively adjust *α* to control the type I error in ABEL and
    RSABE.
  - U.S. <span title="Food and Drug Administration">FDA</span>:
    <span title="Average Bioequivalence">ABE</span> for highly variable
    <span title="Narrow Therapeutic Index Drugs">NTIDs</span> by
    simulations.
  - Dose-Proportionality using the power model.

<small>[TOC ↩](#powertost)</small>

### Methods

  - Exact
      - Owen’s Q.
      - Direct integration of the bivariate non-central
        *t*-distribution.
  - Approximations
      - Non-central *t*-distribution.
      - ‘Shifted’ central *t*-distribution.

<small>[TOC ↩](#powertost)</small>

### Helpers

  - Calculate *CV* from *MSE* or *SE* (and vice versa).
  - Calculate *CV* from given confidence interval.
  - Calculate *CV*<sub>wR</sub> from the upper expanded limit of an ABEL
    study.
  - Confidence interval of *CV*.
  - Pool *CV* from several studies.
  - Confidence interval for given *α*, *CV*, point estimate, sample
    size, and design.
  - Calculate *CV*<sub>wT</sub> and *CV*<sub>wR</sub> from a (pooled)
    *CV*<sub>w</sub> assuming a ratio of intra-subject variances.
  - *p*-values of the <span title="Two One-Sided Tests">TOST</span>
    procedure.
  - Analysis tool for exploration/visualization of the impact of
    expected values (*CV*, *θ*<sub>0</sub>, reduced sample size due to
    dropouts) on power of BE decision.
  - Construct design matrices of incomplete block designs.

<small>[TOC ↩](#powertost)</small>

## Defaults

  - *α* 0.05, {*θ*<sub>1</sub>, *θ*<sub>2</sub>} (0.80, 1.25). Details
    of the sample size search (and the regulatory settings in
    reference-scaled average bioequivalence) are shown in the console.
  - Note: In all functions values have to be given as ratios, not in
    percent.

### Average Bioequivalence

*θ*<sub>0</sub> 0.95, target power 0.80, design "2x2" (TR|RT), exact
method (Owen’s Q).

### Reference-Scaled Average Bioequivalence

*α* 0.05, point estimate constraint (0.80, 1.25), homoscedasticity
(*CV*<sub>wT</sub> = *CV*<sub>wR</sub>), scaling is based on
*CV*<sub>wR</sub>, target power 0.80, design "2x3x3" (TRR|RTR|RRT),
approximation by the non-central *t*-distribution, 100,000 simulations.

  - <span title="European Medicines Agency">EMA</span>,
    <span title="World Health Organization">WHO</span>, Health Canada,
    and many other jurisdictions: Average Bioequivalence with Expanding
    Limits (ABEL).
  - U.S. <span title="Food and Drug Administration">FDA</span>, China
    <span title="Centre for Drug Evaluation">CDE</span>: RSABE.

#### Highly Variable Drugs / Drug Products

*θ*<sub>0</sub> 0.90 as recommended by [Tóthfalusi and
Endrényi](https://ejournals.library.ualberta.ca/index.php/JPPS/article/download/11612/9489)
(2011).

###### EMA

Regulatory constant `0.76`, upper cap of scaling at *CV*<sub>wR</sub>
50%, evaluation by <span title="Analysis of Variance">ANOVA</span>.

###### Health Canada

Regulatory constant `0.76`, upper cap of scaling at *CV*<sub>wR</sub>
\~57.4%, evaluation by intra-subject contrasts.

###### Gulf Cooperation Council

Regulatory constant `log(1/0.75)/sqrt(log(0.3^2+1))`, widened limits
75.00–133.33% if *CV*<sub>wR</sub> \>30%, no upper cap of scaling,
evaluation by <span title="Analysis of Variance">ANOVA</span>.

###### FDA

Regulatory constant `log(1.25)/0.25`, linearized scaled
<span title="Average Bioequivalence">ABE</span> (Howe’s approximation).

#### Narrow Therapeutic Index Drugs (FDA)

*θ*<sub>0</sub> 0.975, regulatory constant `log(1.11111)/0.1`, upper cap
of scaling at *CV*<sub>wR</sub> \~21.4%, design "2x2x4" (TRTR|RTRT),
linearized scaled <span title="Average Bioequivalence">ABE</span>
(Howe’s approximation), upper limit of the confidence interval of
*s*<sub>wT</sub>/*s*<sub>wR</sub> ≤2.5.

### Dose-Proportionality

*β*<sub>0</sub> (slope) `1+log(0.95)/log(rd)` where `rd` is the ratio of
the highest and lowest dose, target power 0.80, crossover design,
details of the sample size search suppressed.

### Power Analysis

Minimum acceptable power 0.70. *θ*<sub>0</sub>, design, conditions, and
sample size method depend on defaults of the respective approaches (ABE,
ABEL, RSABE, NTID, HVNTID).

<small>[TOC ↩](#powertost)</small>

## Examples

Before running the examples attach the library.

``` r
library(PowerTOST)
```

If not noted otherwise, the functions’ [defaults](#defaults) are
employed.

### Parallel Design

Power for total *CV* 0.35 (35%), group sizes 52 and 49, design
"parallel".

``` r
power.TOST(CV = 0.35, n = c(52, 49), design = "parallel")
#R> [1] 0.8011186
```

### Crossover Design

Sample size for assumed within- (intra-) subject *CV* 0.20 (20%).

``` r
sampleN.TOST(CV = 0.20)
#R> 
#R> +++++++++++ Equivalence test - TOST +++++++++++
#R>             Sample size estimation
#R> -----------------------------------------------
#R> Study design: 2x2 crossover 
#R> log-transformed data (multiplicative model)
#R> 
#R> alpha = 0.05, target power = 0.8
#R> BE margins = 0.8 ... 1.25 
#R> True ratio = 0.95,  CV = 0.2
#R> 
#R> Sample size (total)
#R>  n     power
#R> 20   0.834680
```

Sample size for assumed within- (intra-) subject *CV* 0.40 (40%),
*θ*<sub>0</sub> 0.90, four period full replicate “2x2x4” study. Wider
acceptance range for *C*<sub>max</sub> (South Africa).

``` r
sampleN.TOST(CV = 0.40, theta0 = 0.90, theta1 = 0.75, design = "2x2x4")
#R> 
#R> +++++++++++ Equivalence test - TOST +++++++++++
#R>             Sample size estimation
#R> -----------------------------------------------
#R> Study design: 2x2x4 (4 period full replicate) 
#R> log-transformed data (multiplicative model)
#R> 
#R> alpha = 0.05, target power = 0.8
#R> BE margins = 0.75 ... 1.333333 
#R> True ratio = 0.9,  CV = 0.4
#R> 
#R> Sample size (total)
#R>  n     power
#R> 30   0.822929
```

<small>[TOC ↩](#powertost)</small>

Sample size for assumed within- (intra-) subject *CV* 0.125 (12.5%),
*θ*<sub>0</sub> 0.975. Acceptance range for
<span title="Narrow Therapeutic Index Drugs">NTIDs</span> (most
jurisdictions).

``` r
sampleN.TOST(CV = 0.125, theta0 = 0.975, theta1 = 0.90)
#R> 
#R> +++++++++++ Equivalence test - TOST +++++++++++
#R>             Sample size estimation
#R> -----------------------------------------------
#R> Study design: 2x2 crossover 
#R> log-transformed data (multiplicative model)
#R> 
#R> alpha = 0.05, target power = 0.8
#R> BE margins = 0.9 ... 1.111111 
#R> True ratio = 0.975,  CV = 0.125
#R> 
#R> Sample size (total)
#R>  n     power
#R> 32   0.800218
```

<small>[TOC ↩](#powertost)</small>

Sample size for equivalence of the ratio of two means with normality on
the original scale based on [Fieller’s (‘fiducial’) confidence
interval](https://doi.org/10.1111/j.2517-6161.1954.tb00159.x). Within-
(intra-) subject *CV*<sub>w</sub> 0.20 (20%), between- (inter-) subject
*CV*<sub>b</sub> 0.40 (40%).  
Note the default *α* 0.025 (95% CI) of this function because it is
intended for studies with clinical endpoints.

``` r
sampleN.RatioF(CV = 0.20, CVb = 0.40)
#R> 
#R> +++++++++++ Equivalence test - TOST +++++++++++
#R>     based on Fieller's confidence interval
#R>             Sample size estimation
#R> -----------------------------------------------
#R> Study design: 2x2 crossover
#R> Ratio of means with normality on original scale
#R> alpha = 0.025, target power = 0.8
#R> BE margins = 0.8 ... 1.25 
#R> True ratio = 0.95,  CVw = 0.2,  CVb = 0.4
#R> 
#R> Sample size
#R>  n     power
#R> 28   0.807774
```

<small>[TOC ↩](#powertost)</small>

### Replicate Designs

#### ABE

Sample size for assumed within- (intra-) subject *CV* 0.45 (45%),
*θ*<sub>0</sub> 0.90, three period full replicate study "2x2x3"
(TRT|RTR *or* TRR|RTT).

``` r
sampleN.TOST(CV = 0.45, theta0 = 0.90, design = "2x2x3")
#R> 
#R> +++++++++++ Equivalence test - TOST +++++++++++
#R>             Sample size estimation
#R> -----------------------------------------------
#R> Study design: 2x2x3 (3 period full replicate) 
#R> log-transformed data (multiplicative model)
#R> 
#R> alpha = 0.05, target power = 0.8
#R> BE margins = 0.8 ... 1.25 
#R> True ratio = 0.9,  CV = 0.45
#R> 
#R> Sample size (total)
#R>  n     power
#R> 124   0.800125
```

Note that the conventional model assumes homoscedasticity (equal
variances of treatments). For heteroscedasticity we can ‘switch off’ all
conditions of one of the methods for reference-scaled
<span title="Average Bioequivalence">ABE</span>. We assume a
σ<sup>2</sup> ratio of ⅔ (*i.e.*, the test has a lower variability than
the reference). Only relevant columns of the data.frame shown.

``` r
reg <- reg_const("USER", r_const = NA, CVswitch = Inf,
                 CVcap = Inf, pe_constr = FALSE)
CV  <- CVp2CV(CV = 0.45, ratio = 2/3)
res <- sampleN.scABEL(CV=CV, design = "2x2x3", regulator = reg,
                      details = FALSE, print = FALSE)
print(res[c(3:4, 8:9)], digits = 5, row.names = FALSE)
#R>    CVwT    CVwR Sample size Achieved power
#R>  0.3987 0.49767         126         0.8052
```

Similar sample size because the pooled *CV* is still 0.45.

<small>[TOC ↩](#powertost)</small>

#### ABEL

Sample size assuming homoscedasticity (*CV*<sub>wT</sub> =
*CV*<sub>wR</sub> = 0.45).

``` r
sampleN.scABEL(CV = 0.45)
#R> 
#R> +++++++++++ scaled (widened) ABEL +++++++++++
#R>             Sample size estimation
#R>    (simulation based on ANOVA evaluation)
#R> ---------------------------------------------
#R> Study design: 2x3x3 (partial replicate)
#R> log-transformed data (multiplicative model)
#R> 1e+05 studies for each step simulated.
#R> 
#R> alpha  = 0.05, target power = 0.8
#R> CVw(T) = 0.45; CVw(R) = 0.45
#R> True ratio = 0.9
#R> ABE limits / PE constraint = 0.8 ... 1.25 
#R> EMA regulatory settings
#R> - CVswitch            = 0.3 
#R> - cap on scABEL if CVw(R) > 0.5
#R> - regulatory constant = 0.76 
#R> - pe constraint applied
#R> 
#R> 
#R> Sample size search
#R>  n     power
#R> 36   0.7755 
#R> 39   0.8059
```

<small>[TOC ↩](#powertost)</small>

Iteratively adjust *α* to control the Type I Error ([Labes,
Schütz](https://doi.org/10.1007/s11095-016-2006-1)). Slight
heteroscedasticity (*CV*<sub>wT</sub> 0.30, *CV*<sub>wR</sub> 0.35),
four period full replicate "2x2x4" study, 30 subjects, balanced
sequences.

``` r
scABEL.ad(CV = c(0.30, 0.35), design = "2x2x4", n = 30)
#R> +++++++++++ scaled (widened) ABEL ++++++++++++
#R>          iteratively adjusted alpha
#R>    (simulations based on ANOVA evaluation)
#R> ----------------------------------------------
#R> Study design: 2x2x4 (4 period full replicate)
#R> log-transformed data (multiplicative model)
#R> 1,000,000 studies in each iteration simulated.
#R> 
#R> CVwR 0.35, CVwT 0.3, n(i) 15|15 (N 30)
#R> Nominal alpha                 : 0.05 
#R> True ratio                    : 0.9000 
#R> Regulatory settings           : EMA (ABEL)
#R> Switching CVwR                : 0.3 
#R> Regulatory constant           : 0.76 
#R> Expanded limits               : 0.7723 ... 1.2948
#R> Upper scaling cap             : CVwR > 0.5 
#R> PE constraints                : 0.8000 ... 1.2500
#R> Empiric TIE for alpha 0.0500  : 0.06651
#R> Power for theta0 0.9000       : 0.814
#R> Iteratively adjusted alpha    : 0.03540
#R> Empiric TIE for adjusted alpha: 0.05000
#R> Power for theta0 0.9000       : 0.771
```

With the nominal *α* 0.05 the Type I Error will be inflated (0.0665).
With the adjusted *α* 0.0354 (*i.e.*, a 92.92%
<span title="Confidence Interval">CI</span>) the
<span title="Type I Error">TIE</span> will be controlled, although with
a slight loss in power (decreases from 0.814 to 0.771).  
Consider `sampleN.scABEL.ad(CV = c(0.30, 0.35), design = "2x2x4")` to
estimate the sample size which both controls the
<span title="Type I Error">TIE</span> and maintains the target power. In
this example 34 subjects would be required.

<small>[TOC ↩](#powertost)</small>

<span title="Average Bioequivalence with Expanded Limits">ABEL</span>
cannot be applied for *AUC* (except for the
<span title="World Health Organization">WHO<span>). Hence, in many cases
<span title="Average Bioequivalence">ABE</span> drives the sample size.
Three period full replicate "2x2x3" study (TRT|RTR *or* TRR|RTT).

``` r
PK  <- c("Cmax", "AUC")
CV  <- c(0.45, 0.30)
# extract sample sizes and power
r1  <- sampleN.scABEL(CV = CV[1], theta0 = 0.90, design = "2x2x3",
                      print = FALSE, details = FALSE)[8:9]
r2  <- sampleN.TOST(CV = CV[2], theta0 = 0.90, design = "2x2x3",
                    print = FALSE, details = FALSE)[7:8]
n   <- as.numeric(c(r1[1], r2[1]))
pwr <- signif(as.numeric(c(r1[2], r2[2])), 5)
# compile results
res <- data.frame(PK = PK, method = c("ABEL", "ABE"), n = n, power = pwr)
print(res, row.names = FALSE)
#R>    PK method  n   power
#R>  Cmax   ABEL 42 0.80017
#R>   AUC    ABE 60 0.81002
```

<small>[TOC ↩](#powertost)</small>

Sample size assuming homoscedasticity (*CV*<sub>wT</sub> =
*CV*<sub>wR</sub> = 0.45) for the widened limits of the Gulf Cooperation
Council.

``` r
sampleN.scABEL(CV = 0.45, regulator = "GCC", details = FALSE)
#R> 
#R> +++++++++++ scaled (widened) ABEL +++++++++++
#R>             Sample size estimation
#R>    (simulation based on ANOVA evaluation)
#R> ---------------------------------------------
#R> Study design: 2x3x3 (partial replicate)
#R> log-transformed data (multiplicative model)
#R> 1e+05 studies for each step simulated.
#R> 
#R> alpha  = 0.05, target power = 0.8
#R> CVw(T) = 0.45; CVw(R) = 0.45
#R> True ratio = 0.9
#R> ABE limits / PE constraint = 0.8 ... 1.25 
#R> Widened limits = 0.75 ... 1.333333 
#R> Regulatory settings: GCC 
#R> 
#R> Sample size
#R>  n     power
#R> 54   0.8123
```

<small>[TOC ↩](#powertost)</small>

#### RSABE

#### HVD(P)s

Sample size for a four period full replicate "2x2x4" study (any of
TRTR|RTRT, TRRT|RTTR, TTRR|RRTT) assuming heteroscedasticity
(*CV*<sub>wT</sub> 0.40, *CV*<sub>wR</sub> 0.50). Details of the sample
size search suppressed.

``` r
sampleN.RSABE(CV = c(0.40, 0.50), design = "2x2x4", details = FALSE)
#R> 
#R> ++++++++ Reference scaled ABE crit. +++++++++
#R>            Sample size estimation
#R> ---------------------------------------------
#R> Study design: 2x2x4 (4 period full replicate)
#R> log-transformed data (multiplicative model)
#R> 1e+05 studies for each step simulated.
#R> 
#R> alpha  = 0.05, target power = 0.8
#R> CVw(T) = 0.4; CVw(R) = 0.5
#R> True ratio = 0.9
#R> ABE limits / PE constraints = 0.8 ... 1.25 
#R> Regulatory settings: FDA 
#R> 
#R> Sample size
#R>  n    power
#R> 20   0.81509
```

<small>[TOC ↩](#powertost)</small>

#### NTIDs

Sample size assuming heteroscedasticity (*CV*<sub>w</sub> 0.10,
σ<sup>2</sup> ratio 2.5, *i.e.*, the test treatment has a substantially
higher variability than the reference). TRTR|RTRT according to the
[FDA’s
guidance](https://www.accessdata.fda.gov/drugsatfda_docs/psg/Warfarin_Sodium_tab_09218_RC12-12.pdf).
Assess additionally which one of the three components (scaled
<span title="Average Bioequivalence">ABE</span>, conventional
<span title="Average Bioequivalence">ABE</span>,
*s*<sub>wT</sub>/*s*<sub>wR</sub> ratio) drives the sample size.

``` r
CV <- signif(CVp2CV(CV = 0.10, ratio = 2.5), 4)
n  <- sampleN.NTIDFDA(CV = CV)[["Sample size"]]
#R> 
#R> +++++++++++ FDA method for NTIDs ++++++++++++
#R>            Sample size estimation
#R> ---------------------------------------------
#R> Study design:  2x2x4 (TRTR|RTRT) 
#R> log-transformed data (multiplicative model)
#R> 1e+05 studies for each step simulated.
#R> 
#R> alpha  = 0.05, target power = 0.8
#R> CVw(T) = 0.1197, CVw(R) = 0.07551
#R> True ratio     = 0.975 
#R> ABE limits     = 0.8 ... 1.25 
#R> Implied scABEL = 0.9236 ... 1.0827 
#R> Regulatory settings: FDA 
#R> - Regulatory const. = 1.053605 
#R> - 'CVcap'           = 0.2142 
#R> 
#R> Sample size search
#R>  n     power
#R> 32   0.699120 
#R> 34   0.730910 
#R> 36   0.761440 
#R> 38   0.785910 
#R> 40   0.809580
suppressMessages(power.NTIDFDA(CV = CV, n = n, details = TRUE))
#R>        p(BE)  p(BE-sABEc)    p(BE-ABE) p(BE-sratio) 
#R>      0.80958      0.90966      1.00000      0.87447
```

The *s*<sub>wT</sub>/*s*<sub>wR</sub> component shows the lowest
probability to demonstrate <span title="Bioequivalence">BE</span> and
hence, drives the sample size.

<small>[TOC ↩](#powertost)</small>

Compare that with homoscedasticity (*CV*<sub>wT</sub> =
*CV*<sub>wR</sub> = 0.10):

``` r
CV <- 0.10
n  <- sampleN.NTIDFDA(CV = CV, details = FALSE)[["Sample size"]]
#R> 
#R> +++++++++++ FDA method for NTIDs ++++++++++++
#R>            Sample size estimation
#R> ---------------------------------------------
#R> Study design:  2x2x4 (TRTR|RTRT) 
#R> log-transformed data (multiplicative model)
#R> 1e+05 studies for each step simulated.
#R> 
#R> alpha  = 0.05, target power = 0.8
#R> CVw(T) = 0.1, CVw(R) = 0.1
#R> True ratio     = 0.975 
#R> ABE limits     = 0.8 ... 1.25 
#R> Regulatory settings: FDA 
#R> 
#R> Sample size
#R>  n     power
#R> 18   0.841790
suppressMessages(power.NTIDFDA(CV = CV, n = n, details = TRUE))
#R>        p(BE)  p(BE-sABEc)    p(BE-ABE) p(BE-sratio) 
#R>      0.84179      0.85628      1.00000      0.97210
```

Here the scaled <span title="Average Bioequivalence">ABE</span>
component shows the lowest probability to demonstrate
<span title="Bioequivalence">BE</span> and drives the sample size –
which is much lower than in the previous example.

<small>[TOC ↩](#powertost)</small>

Comparison with *fixed* narrower limits applicable in other
jurisdictions. Note that a replicate design is not required, reducing
the chance of dropouts.

``` r
CV  <- 0.10
# extract sample sizes and power
r1  <- sampleN.NTIDFDA(CV = CV, print = FALSE, details = FALSE)[8:9]
r2  <- sampleN.TOST(CV = CV, theta0 = 0.975, theta1 = 0.90,
                    design = "2x2x4", print = FALSE, details = FALSE)[7:8]
r3  <- sampleN.TOST(CV = CV, theta0 = 0.975, theta1 = 0.90,
                    design = "2x2x3", print = FALSE, details = FALSE)[7:8]
r4  <- sampleN.TOST(CV = CV, theta0 = 0.975, theta1 = 0.90,
                    print = FALSE, details = FALSE)[7:8]
n   <- as.numeric(c(r1[1], r2[1], r3[1], r4[1]))
pwr <- signif(as.numeric(c(r1[2], r2[2], r3[2], r4[2])), 5)
# compile results
res <- data.frame(method = c("FDA scaled", rep ("fixed narrow", 3)),
                  design = c(rep("2x2x4", 2), "2x2x3", "2x2x2"),
                  n = n, power = pwr, a = n * c(4, 4, 3, 2))
names(res)[5] <- "adm. #"
print(res, row.names = FALSE)
#R>        method design  n   power adm. #
#R>    FDA scaled  2x2x4 18 0.84179     72
#R>  fixed narrow  2x2x4 12 0.85628     48
#R>  fixed narrow  2x2x3 16 0.81393     48
#R>  fixed narrow  2x2x2 22 0.81702     44
```

<small>[TOC ↩](#powertost)</small>

### Dose-Proportionality

*CV* 0.20 (20%), doses 1, 2, and 8 units, assumed slope *β*<sub>0</sub>
1, target power 0.90.

``` r
sampleN.dp(CV = 0.20, doses = c(1, 2, 8), beta0 = 1, targetpower = 0.90)
#R> 
#R> ++++ Dose proportionality study, power model ++++
#R>             Sample size estimation
#R> -------------------------------------------------
#R> Study design: crossover (3x3 Latin square) 
#R> alpha = 0.05, target power = 0.9
#R> Equivalence margins of R(dnm) = 0.8 ... 1.25 
#R> Doses = 1 2 8 
#R> True slope = 1, CV = 0.2
#R> Slope acceptance range = 0.89269 ... 1.1073 
#R> 
#R> Sample size (total)
#R>  n     power
#R> 18   0.915574
```

Note that the acceptance range of the slope depends on the ratio of the
highest and lowest doses (*i.e.*, it gets tighter for wider dose ranges
and therefore, higher sample sizes will be required).  
In an exploratory setting wider equivalence margins {*θ*<sub>1</sub>,
*θ*<sub>2</sub>} (0.50, 2.00) were
[proposed](https://doi.org/10.1002/pst.326), translating in this example
to an acceptance range of `0.66667 ... 1.3333` and a sample size of only
six subjects.

<small>[TOC ↩](#powertost)</small>

### Power Analysis

Explore impact of deviations from assumptions (higher *CV*, higher
deviation of *θ*<sub>0</sub> from 1, dropouts) on power. Assumed
within-subject *CV* 0.20 (20%), target power 0.90. Plot suppressed.

``` r
res <- pa.ABE(CV = 0.20, targetpower = 0.90)
print(res, plotit = FALSE)
#R> Sample size plan ABE
#R>  Design alpha  CV theta0 theta1 theta2 Sample size Achieved power
#R>     2x2  0.05 0.2   0.95    0.8   1.25          26      0.9176333
#R> 
#R> Power analysis
#R> CV, theta0 and number of subjects leading to min. acceptable power of =0.7:
#R>  CV= 0.2729, theta0= 0.9044
#R>  n = 16 (power= 0.7354)
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

<small>[TOC ↩](#powertost)</small>

### Speed Comparisons

Performed on a Xeon E3-1245v3 3.4 GHz, 8 MB cache, 16 GB RAM, R 4.1.0
64 bit on Windows 7.

#### ABE

"2x2" crossover design, *CV* 0.17. Sample sizes and achieved power for
the supported methods (the 1<sup>st</sup> one is the default).

``` 
    method  n   power time (s)
     owenq 14 0.80568  0.00128
       mvt 14 0.80569  0.11778
noncentral 14 0.80568  0.00100
   shifted 16 0.85230  0.00096
```

The 2<sup>nd</sup> exact method is substantially slower than the
1<sup>st</sup>. The approximation based on the noncentral
*t*-distribution is slightly faster but matches the 1<sup>st</sup> exact
method closely. Though the approximation based on the shifted central
*t*-distribution is the fastest, it *might* estimate a larger than
necessary sample size. Hence, it should be used only for comparative
purposes.

#### ABEL

Four period full replicate study, homogenicity (*CV*<sub>wT</sub> =
*CV*<sub>wR</sub> 0.45). Sample sizes and achieved power for the
supported methods.

``` 
              function              method  n   power time (s)
        sampleN.scABEL    ‘key’ statistics 28 0.81116   0.1348
 sampleN.scABEL.sdsims subject simulations 28 0.81196   2.5377
```

Simulating via the ‘key’ statistics is the method of choice for speed
reasons.  
However, subject simulations are recommended *if*

  - the partial replicate design (TRR|RTR|RRT) is planned **and**
  - the special case of heterogenicity *CV*<sub>wT</sub> \>
    *CV*<sub>wR</sub> is expected.

<small>[TOC ↩](#powertost)</small>

## Installation

You can install the released version of PowerTOST from
[CRAN](https://CRAN.R-project.org) with

``` r
package <- "PowerTOST"
inst    <- package %in% installed.packages()
if (length(package[!inst]) > 0) install.packages(package[!inst])
```

… and the development version from [GitHub](https://github.com/) with

    # install.packages("remotes")
    remotes::install_github("Detlew/PowerTOST")

Skips installation from a github remote if the
[SHA-1](https://en.wikipedia.org/wiki/SHA-1) has not changed since last
install. Use `force = TRUE` to force installation.

<small>[TOC ↩](#powertost)</small>

## Session Information

Inspect this information for reproducibility. Of particular importance
are the versions of R and the packages used to create this workflow. It
is considered good practice to record this information with every
analysis.  
Version 1.5.3.9000 built 2021-05-26 with R 4.1.0.

``` r
options(width = 60)
devtools::session_info()
#R> - Session info -------------------------------------------
#R>  setting  value                       
#R>  version  R version 4.1.0 (2021-05-18)
#R>  os       Windows 7 x64 SP 1          
#R>  system   x86_64, mingw32             
#R>  ui       RTerm                       
#R>  language EN                          
#R>  collate  German_Germany.1252         
#R>  ctype    German_Germany.1252         
#R>  tz       Europe/Vienna               
#R>  date     2021-05-26                  
#R> 
#R> - Packages -----------------------------------------------
#R>  package       * version    date       lib source        
#R>  cachem          1.0.5      2021-05-15 [1] CRAN (R 4.1.0)
#R>  callr           3.7.0      2021-04-20 [1] CRAN (R 4.1.0)
#R>  cli             2.5.0      2021-04-26 [1] CRAN (R 4.1.0)
#R>  crayon          1.4.1      2021-02-08 [1] CRAN (R 4.1.0)
#R>  cubature        2.0.4.2    2021-05-13 [1] CRAN (R 4.1.0)
#R>  desc            1.3.0      2021-03-05 [1] CRAN (R 4.1.0)
#R>  devtools        2.4.1      2021-05-05 [1] CRAN (R 4.1.0)
#R>  digest          0.6.27     2020-10-24 [1] CRAN (R 4.1.0)
#R>  ellipsis        0.3.2      2021-04-29 [1] CRAN (R 4.1.0)
#R>  evaluate        0.14       2019-05-28 [1] CRAN (R 4.1.0)
#R>  fastmap         1.1.0      2021-01-25 [1] CRAN (R 4.1.0)
#R>  fs              1.5.0      2020-07-31 [1] CRAN (R 4.1.0)
#R>  glue            1.4.2      2020-08-27 [1] CRAN (R 4.1.0)
#R>  htmltools       0.5.1.1    2021-01-22 [1] CRAN (R 4.1.0)
#R>  knitr           1.33       2021-04-24 [1] CRAN (R 4.1.0)
#R>  lifecycle       1.0.0      2021-02-15 [1] CRAN (R 4.1.0)
#R>  magrittr        2.0.1      2020-11-17 [1] CRAN (R 4.1.0)
#R>  memoise         2.0.0      2021-01-26 [1] CRAN (R 4.1.0)
#R>  mvtnorm         1.1-1      2020-06-09 [1] CRAN (R 4.1.0)
#R>  pkgbuild        1.2.0      2020-12-15 [1] CRAN (R 4.1.0)
#R>  pkgload         1.2.1      2021-04-06 [1] CRAN (R 4.1.0)
#R>  PowerTOST     * 1.5.3.9000 2021-05-26 [1] local         
#R>  prettyunits     1.1.1      2020-01-24 [1] CRAN (R 4.1.0)
#R>  processx        3.5.2      2021-04-30 [1] CRAN (R 4.1.0)
#R>  ps              1.6.0      2021-02-28 [1] CRAN (R 4.1.0)
#R>  purrr           0.3.4      2020-04-17 [1] CRAN (R 4.1.0)
#R>  R6              2.5.0      2020-10-28 [1] CRAN (R 4.1.0)
#R>  Rcpp            1.0.6      2021-01-15 [1] CRAN (R 4.1.0)
#R>  remotes         2.3.0      2021-04-01 [1] CRAN (R 4.1.0)
#R>  rlang           0.4.11     2021-04-30 [1] CRAN (R 4.1.0)
#R>  rmarkdown       2.8        2021-05-07 [1] CRAN (R 4.1.0)
#R>  rprojroot       2.0.2      2020-11-15 [1] CRAN (R 4.1.0)
#R>  sessioninfo     1.1.1      2018-11-05 [1] CRAN (R 4.1.0)
#R>  stringi         1.6.1      2021-05-10 [1] CRAN (R 4.1.0)
#R>  stringr         1.4.0      2019-02-10 [1] CRAN (R 4.1.0)
#R>  TeachingDemos   2.12       2020-04-07 [1] CRAN (R 4.1.0)
#R>  testthat        3.0.2      2021-02-14 [1] CRAN (R 4.1.0)
#R>  usethis         2.0.1      2021-02-10 [1] CRAN (R 4.1.0)
#R>  withr           2.4.2      2021-04-18 [1] CRAN (R 4.1.0)
#R>  xfun            0.23       2021-05-15 [1] CRAN (R 4.1.0)
#R>  yaml            2.2.1      2020-02-01 [1] CRAN (R 4.1.0)
#R> 
#R> [1] D:/Program Files/R/R-4.1.0/library
```

<small>[TOC ↩](#powertost)</small>
