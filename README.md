PowerTOST
================

- [Introduction](#introduction)
- [Supported Designs](#designs)
- [Purpose](#purpose)
- [Supported](#supported)
  - [Power and Sample Size](#power_sample_size)
  - [Methods](#methods)
  - [Helpers](#helpers)
- [Defaults](#defaults)
  - [Average Bioequivalence](#ABE.defaults)
  - [Reference-Scaled Average Bioequivalence](#RSABE.defaults)
  - [Dose-Proportionality](#DP)
  - [Power Analysis](#PA)
- [Examples](#examples)
  - [Parallel Design](#parallel)
  - [Crossover Design](#crossover)
  - [Replicate Designs](#replicate)
  - [Dose-Proportionality](#dose-proportionality)
  - [Power Analysis](#power-analysis)
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
checks](https://badges.cranchecks.info/worst/PowerTOST.svg)](https://cran.r-project.org/web/checks/check_results_PowerTOST.html)
![commits](https://img.shields.io/github/commits-since/Detlew/PowerTOST/latest)
![last
commit](https://img.shields.io/github/last-commit/Detlew/PowerTOST)
[![CRAN RStudio mirror
downloads](https://cranlogs.r-pkg.org/badges/grand-total/PowerTOST?color=blue)](https://r-pkg.org/pkg/PowerTOST)
[![CRAN RStudio mirror
downloads](https://cranlogs.r-pkg.org/badges/last-month/PowerTOST?color=green)](https://r-pkg.org/pkg/PowerTOST)

Version 1.5.6 built 2024-03-18 with R 4.3.3

## Introduction

The package contains functions to calculate power and estimate sample
size for various study designs used in (not only bio-) equivalence
studies.

## Supported Designs

    #    design                        name    df
    #  parallel           2 parallel groups   n-2
    #       2x2               2x2 crossover   n-2
    #     2x2x2             2x2x2 crossover   n-2
    #       3x3               3x3 crossover 2*n-4
    #     3x6x3             3x6x3 crossover 2*n-4
    #       4x4               4x4 crossover 3*n-6
    #     2x2x3   2x2x3 replicate crossover 2*n-3
    #     2x2x4   2x2x4 replicate crossover 3*n-4
    #     2x4x4   2x4x4 replicate crossover 3*n-4
    #     2x3x3   partial replicate (2x3x3) 2*n-3
    #     2x4x2            Balaam's (2x4x2)   n-2
    #    2x2x2r Liu's 2x2x2 repeated x-over 3*n-2
    #    paired                paired means   n-1

Codes of designs follow this pattern:
`treatments x sequences x periods`.

Although some replicate designs are more ‘popular’ than others, sample
size estimations are valid for *all* of the following designs:

| <small>design</small> | <small>type</small> | <small>sequences</small> |                                       | <small>periods</small> |
|:---------------------:|:-------------------:|:------------------------:|---------------------------------------|:----------------------:|
|        `2x2x4`        |     <small>full     |     <small>2</small>     | <small>TRTR\|RTRT</small>             |    <small>4</small>    |
|        `2x2x4`        |     <small>full     |     <small>2</small>     | <small>TRRT\|RTTR</small>             |    <small>4</small>    |
|        `2x2x4`        |     <small>full     |     <small>2</small>     | <small>TTRR\|RRTT</small>             |    <small>4</small>    |
|        `2x4x4`        |     <small>full     |     <small>4</small>     | <small>TRTR\|RTRT\|TRRT\|RTTR</small> |    <small>4</small>    |
|        `2x4x4`        |     <small>full     |     <small>4</small>     | <small>TRRT\|RTTR\|TTRR\|RRTT</small> |    <small>4</small>    |
|        `2x2x3`        |     <small>full     |     <small>2</small>     | <small>TRT\|RTR</small>               |    <small>3</small>    |
|        `2x2x3`        |     <small>full     |     <small>2</small>     | <small>TRR\|RTT</small>               |    <small>3</small>    |
|        `2x4x2`        |     <small>full     |     <small>4</small>     | <small>TR\|RT\|TT\|RR</small>         |    <small>2</small>    |
|        `2x3x3`        |   <small>partial    |     <small>3</small>     | <small>TRR\|RTR\|RRT</small>          |    <small>3</small>    |
|        `2x2x3`        |   <small>partial    |     <small>2</small>     | <small>TRR\|RTR</small>               |    <small>3</small>    |

Balaam’s design TR\|RT\|TT\|RR should be avoided due to its poor power
characteristics. The three period partial replicate design with two
sequences TRR\|RTR (<span title="also known as">a.k.a.</span>
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
- <span title="Average Bioequivalence">ABE</span> for Highly Variable
  Narrow Therapeutic Index Drugs by simulations: U.S.
  <span title="Food and Drug Administration">FDA</span>, China
  <span title="Center for Drug Evaluation">CDE</span>.
- Scaled Average Bioequivalence based on simulations.
  - Average Bioequivalence with Expanding Limits (ABEL) for Highly
    Variable Drugs / Drug Products:
    <span title="European Medicines Agency">EMA</span>,
    <span title="World Health Organization">WHO</span> and many others.
  - Average Bioequivalence with *fixed* widened limits of 75.00–133.33%
    if *CV*<sub>wR</sub> \>30%: Gulf Cooperation Council.  
  - Reference-scaled Average Bioequivalence (RSABE) for
    <span title="Highly Variable Drugs / Drug Products">HVDP(s)</span>:
    U.S. <span title="Food and Drug Administration">FDA</span>, China
    <span title="Center for Drug Evaluation">CDE</span>.
  - Iteratively adjust *α* to control the type I error in
    <span title="Average Bioequivalence with Expanding Limits">ABEL</span>
    and
    <span title="Reference-scaled Average Bioequivalence">RSABE</span>
    for
    <span title="Highly Variable Drugs / Drug Products">HVDP(s)</span>.
  - <span title="Reference-scaled Average Bioequivalence">RSABE</span>
    for <span title="Narrow Therapeutic Index Drugs">NTIDs</span>: U.S.
    <span title="Food and Drug Administration">FDA</span>, China
    <span title="Center for Drug Evaluation">CDE</span>.
- Two simultaneous <span title="Two One-Sided Tests">TOST</span>
  procedures.
- Non-inferiority *t*-test.
- Ratio of two means with normally distributed data on the original
  scale based on Fieller’s (‘fiducial’) confidence interval.
- ‘Expected’ power in case of uncertain (estimated) variability and/or
  uncertain *θ*<sub>0</sub>.
- Dose-Proportionality using the power model.

<small>[TOC ↩](#powertost)</small>

### Methods

- Exact
  - Owen’s Q.
  - Direct integration of the bivariate non-central *t*-distribution.
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
- Confidence interval for given *α*, *CV*, point estimate, sample size,
  and design.
- Calculate *CV*<sub>wT</sub> and *CV*<sub>wR</sub> from a (pooled)
  *CV*<sub>w</sub> assuming a ratio of intra-subject variances.
- *p*-values of the <span title="Two One-Sided Tests">TOST</span>
  procedure.
- Analysis tool for exploration/visualization of the impact of expected
  values (*CV*, *θ*<sub>0</sub>, reduced sample size due to dropouts) on
  power of BE decision.
- Construct design matrices of incomplete block designs.

<small>[TOC ↩](#powertost)</small>

## Defaults

- *α* 0.05, {*θ*<sub>1</sub>, *θ*<sub>2</sub>} (0.80, 1.25), target
  power 0.80. Details of the sample size search (and the regulatory
  settings in reference-scaled average bioequivalence) are shown in the
  console.
- Note: In all functions values have to be given as ratios, not in
  percent.

### Average Bioequivalence

#### Conventional (unscaled)

Design `"2x2"` (TR\|RT), exact method (Owen’s Q).

#### Highly Variable NTIDs (FDA, CDE)

Design `"2x2x4"` (TRTR\|RTRT), upper limit of the confidence interval of
*σ*<sub>wT</sub>/*σ*<sub>wR</sub> ≤2.5, approximation by the non-central
*t*-distribution, 100,000 simulations.

### Reference-Scaled Average Bioequivalence

Point estimate constraints (0.80, 1.25), homoscedasticity
(*CV*<sub>wT</sub> = *CV*<sub>wR</sub>), scaling is based on
*CV*<sub>wR</sub>, design `"2x3x3"` (TRR\|RTR\|RRT), approximation by
the non-central *t*-distribution, 100,000 simulations.

- <span title="European Medicines Agency">EMA</span>,
  <span title="World Health Organization">WHO</span>, Health Canada, and
  many other jurisdictions: Average Bioequivalence with Expanding Limits
  (ABEL).
- U.S. <span title="Food and Drug Administration">FDA</span>, China
  <span title="Center for Drug Evaluation">CDE</span>: RSABE.

#### Highly Variable Drugs / Drug Products

*θ*<sub>0</sub> 0.90.<sup id="a1">[1](#f1)</sup>

###### EMA and many others

Regulatory constant `0.760`, upper cap of scaling at *CV*<sub>wR</sub>
50%, evaluation by <span title="Analysis of Variance">ANOVA</span>.

###### Health Canada

Regulatory constant `0.760`, upper cap of scaling at *CV*<sub>wR</sub>
~57.4%, evaluation by intra-subject contrasts.

###### Gulf Cooperation Council

Regulatory constant `log(1/0.75)/sqrt(log(0.3^2+1))`, widened limits
75.00–133.33% if *CV*<sub>wR</sub> \>30%, no upper cap of scaling,
evaluation by <span title="Analysis of Variance">ANOVA</span>.

###### FDA, CDE

Regulatory constant `log(1.25)/0.25`, no upper cap of scaling,
evaluation by linearized scaled
<span title="Average Bioequivalence">ABE</span> (Howe’s approximation).

#### Narrow Therapeutic Index Drugs (FDA, CDE)

*θ*<sub>0</sub> 0.975, regulatory constant `log(1.11111)/0.1`, implicit
upper cap of scaling at *CV*<sub>wR</sub> ~21.4%, design `"2x2x4"`
(TRTR\|RTRT), evaluation by linearized scaled
<span title="Average Bioequivalence">ABE</span> (Howe’s approximation),
upper limit of the confidence interval of
*σ*<sub>wT</sub>/*σ*<sub>wR</sub> ≤2.5.

### Dose-Proportionality

*β*<sub>0</sub> (slope) `1+log(0.95)/log(rd)` where `rd` is the ratio of
the highest and lowest dose, target power 0.80, crossover design,
details of the sample size search suppressed.

### Power Analysis

Minimum acceptable power 0.70. *θ*<sub>0</sub>; design, conditions, and
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

Power for total *CV* 0.35 (35%), group sizes 52 and 49.

``` r
power.TOST(CV = 0.35, n = c(52, 49), design = "parallel")
# [1] 0.8011186
```

### Crossover Design

Sample size for assumed within- (intra-) subject *CV* 0.20 (20%).

``` r
sampleN.TOST(CV = 0.20)
# 
# +++++++++++ Equivalence test - TOST +++++++++++
#             Sample size estimation
# -----------------------------------------------
# Study design: 2x2 crossover 
# log-transformed data (multiplicative model)
# 
# alpha = 0.05, target power = 0.8
# BE margins = 0.8 ... 1.25 
# True ratio = 0.95,  CV = 0.2
# 
# Sample size (total)
#  n     power
# 20   0.834680
```

<small>[TOC ↩](#powertost)</small>

Sample size for assumed within- (intra-) subject *CV* 0.40 (40%),
*θ*<sub>0</sub> 0.90, four period full replicate study (any of
TRTR\|RTRT, TRRT\|RTTR, TTRR\|RRTT). Wider acceptance range for
*C*<sub>max</sub> (South Africa).

``` r
sampleN.TOST(CV = 0.40, theta0 = 0.90, theta1 = 0.75, design = "2x2x4")
# 
# +++++++++++ Equivalence test - TOST +++++++++++
#             Sample size estimation
# -----------------------------------------------
# Study design: 2x2x4 (4 period full replicate) 
# log-transformed data (multiplicative model)
# 
# alpha = 0.05, target power = 0.8
# BE margins = 0.75 ... 1.333333 
# True ratio = 0.9,  CV = 0.4
# 
# Sample size (total)
#  n     power
# 30   0.822929
```

<small>[TOC ↩](#powertost)</small>

Sample size for assumed within- (intra-) subject *CV* 0.125 (12.5%),
*θ*<sub>0</sub> 0.975. Narrower acceptance range for
<span title="Narrow Therapeutic Index Drugs">NTIDs</span> (most
jurisdictions).

``` r
sampleN.TOST(CV = 0.125, theta0 = 0.975, theta1 = 0.90)
# 
# +++++++++++ Equivalence test - TOST +++++++++++
#             Sample size estimation
# -----------------------------------------------
# Study design: 2x2 crossover 
# log-transformed data (multiplicative model)
# 
# alpha = 0.05, target power = 0.8
# BE margins = 0.9 ... 1.111111 
# True ratio = 0.975,  CV = 0.125
# 
# Sample size (total)
#  n     power
# 32   0.800218
```

<small>[TOC ↩](#powertost)</small>

Sample size for equivalence of the ratio of two means with normality on
the original scale based on Fieller’s (‘fiducial’) confidence
interval.<sup id="a2">[2](#f2)</sup> Within- (intra-) subject
*CV*<sub>w</sub> 0.20 (20%), between- (inter-) subject *CV*<sub>b</sub>
0.40 (40%).  
Note the default *α* 0.025 (95% CI) of this function because it is
intended for studies with clinical endpoints.

``` r
sampleN.RatioF(CV = 0.20, CVb = 0.40)
# 
# +++++++++++ Equivalence test - TOST +++++++++++
#     based on Fieller's confidence interval
#             Sample size estimation
# -----------------------------------------------
# Study design: 2x2 crossover
# Ratio of means with normality on original scale
# alpha = 0.025, target power = 0.8
# BE margins = 0.8 ... 1.25 
# True ratio = 0.95,  CVw = 0.2,  CVb = 0.4
# 
# Sample size
#  n     power
# 28   0.807774
```

<small>[TOC ↩](#powertost)</small>

### Replicate Designs

#### ABE

##### Conventional (unscaled)

Sample size for assumed within- (intra-) subject *CV* 0.45 (45%),
*θ*<sub>0</sub> 0.90, three period full replicate study (TRT\|RTR *or*
TRR\|RTT).

``` r
sampleN.TOST(CV = 0.45, theta0 = 0.90, design = "2x2x3")
# 
# +++++++++++ Equivalence test - TOST +++++++++++
#             Sample size estimation
# -----------------------------------------------
# Study design: 2x2x3 (3 period full replicate) 
# log-transformed data (multiplicative model)
# 
# alpha = 0.05, target power = 0.8
# BE margins = 0.8 ... 1.25 
# True ratio = 0.9,  CV = 0.45
# 
# Sample size (total)
#  n     power
# 124   0.800125
```

Note that the conventional model assumes homoscedasticity (equal
variances of treatments). For heteroscedasticity we can ‘switch off’ all
conditions of one of the methods for reference-scaled
<span title="Average Bioequivalence">ABE</span>. We assume a
*σ*<sup>2</sup>-ratio of ⅔ (*i.e.*, the test has a lower variability
than the reference). Only relevant columns of the data frame shown.

``` r
reg <- reg_const("USER", r_const = NA, CVswitch = Inf,
                 CVcap = Inf, pe_constr = FALSE)
CV  <- CVp2CV(CV = 0.45, ratio = 2/3)
res <- sampleN.scABEL(CV=CV, design = "2x2x3", regulator = reg,
                      details = FALSE, print = FALSE)
print(res[c(3:4, 8:9)], digits = 5, row.names = FALSE)
#    CVwT    CVwR Sample size Achieved power
#  0.3987 0.49767         126         0.8052
```

Similar sample size because the pooled *CV*<sub>w</sub> is still 0.45.

<small>[TOC ↩](#powertost)</small>

##### Highly Variable Narrow Therapeutic Index Drug

Sample size assuming heteroscedasticity (*CV*<sub>w</sub> 0.45,
variance-ratio 2.5, *i.e.*, the test treatment has a substantially
higher variability than the reference). TRTR\|RTRT according to the
FDA’s
guidances.<sup id="a3">[3,](#f3)</sup><sup id="a4">[4,](#f4)</sup><sup id="a5">[5](#f5)</sup>
Assess additionally which one of the components
(<span title="Average Bioequivalence">ABE</span>,
*s*<sub>wT</sub>/*s*<sub>wR</sub>-ratio) drives the sample size.

``` r
CV <- signif(CVp2CV(CV = 0.45, ratio = 2.5), 4)
n  <- sampleN.HVNTID(CV = CV, details = FALSE)[["Sample size"]]
# 
# +++++++++ FDA method for HV NTIDs ++++++++++++
#            Sample size estimation
# ----------------------------------------------
# Study design: 2x2x4 (TRTR|RTRT) 
# log-transformed data (multiplicative model)
# 1e+05 studies for each step simulated.
# 
# alpha  = 0.05, target power = 0.8
# CVw(T) = 0.549, CVw(R) = 0.3334
# True ratio = 0.95 
# ABE limits = 0.8 ... 1.25 
# 
# Sample size
#  n     power
# 50   0.812820
suppressMessages(power.HVNTID(CV = CV, n = n, details = TRUE))
#        p(BE)    p(BE-ABE) p(BE-sratio) 
#      0.81282      0.87052      0.93379
```

The ABE component shows a lower probability to demonstrate BE than the
*s*<sub>wT</sub>/*s*<sub>wR</sub> component and hence, drives the sample
size.

<small>[TOC ↩](#powertost)</small>

#### ABEL

Sample size assuming homoscedasticity (*CV*<sub>wT</sub> =
*CV*<sub>wR</sub> = 0.45).

``` r
sampleN.scABEL(CV = 0.45)
# 
# +++++++++++ scaled (widened) ABEL +++++++++++
#             Sample size estimation
#    (simulation based on ANOVA evaluation)
# ---------------------------------------------
# Study design: 2x3x3 (partial replicate)
# log-transformed data (multiplicative model)
# 1e+05 studies for each step simulated.
# 
# alpha  = 0.05, target power = 0.8
# CVw(T) = 0.45; CVw(R) = 0.45
# True ratio = 0.9
# ABE limits / PE constraint = 0.8 ... 1.25 
# EMA regulatory settings
# - CVswitch            = 0.3 
# - cap on scABEL if CVw(R) > 0.5
# - regulatory constant = 0.76 
# - pe constraint applied
# 
# 
# Sample size search
#  n     power
# 36   0.7755 
# 39   0.8059
```

<small>[TOC ↩](#powertost)</small>

Iteratively adjust *α* to control the Type I
Error.<sup id="a6">[6](#f6)</sup> Heteroscedasticity
(*CV*<sub>wT</sub> 0.30, *CV*<sub>wR</sub> 0.40, *i.e.*, variance-ratio
~0.58), four period full replicate study (any of TRTR\|RTRT, TRRT\|RTTR,
TTRR\|RRTT), 24 subjects, balanced sequences.

``` r
scABEL.ad(CV = c(0.30, 0.40), design = "2x2x4", n = 24)
# +++++++++++ scaled (widened) ABEL ++++++++++++
#          iteratively adjusted alpha
#    (simulations based on ANOVA evaluation)
# ----------------------------------------------
# Study design: 2x2x4 (4 period full replicate)
# log-transformed data (multiplicative model)
# 1,000,000 studies in each iteration simulated.
# 
# CVwR 0.4, CVwT 0.3, n(i) 12|12 (N 24)
# Nominal alpha                 : 0.05 
# True ratio                    : 0.9000 
# Regulatory settings           : EMA (ABEL)
# Switching CVwR                : 0.3 
# Regulatory constant           : 0.76 
# Expanded limits               : 0.7462 ... 1.3402
# Upper scaling cap             : CVwR > 0.5 
# PE constraints                : 0.8000 ... 1.2500
# Empiric TIE for alpha 0.0500  : 0.05953
# Power for theta0 0.9000       : 0.805
# Iteratively adjusted alpha    : 0.03997
# Empiric TIE for adjusted alpha: 0.05000
# Power for theta0 0.9000       : 0.778
```

With the nominal *α* 0.05 the Type I Error will be inflated (0.05953).
With the adjusted *α* 0.03997 (*i.e.*, a ~92%
<span title="Confidence Interval">CI</span>) the
<span title="Type I Error">TIE</span> will be controlled, although with
a slight loss in power (decreases from 0.805 to 0.778).  
Consider `sampleN.scABEL.ad(CV = c(0.30, 0.35), design = "2x2x4")` to
estimate the sample size preserving both the
<span title="Type I Error">TIE</span> and target power. In this example
26 subjects would be required.

<small>[TOC ↩](#powertost)</small>

<span title="Average Bioequivalence with Expanded Limits">ABEL</span>
cannot be applied for *AUC* (except for the
<span title="World Health Organization">WHO<span>). Hence, in many cases
<span title="Average Bioequivalence">ABE</span> drives the sample size.
Four period full replicate study (any of TRTR\|RTRT, TRRT\|RTTR,
TTRR\|RRTT).

``` r
PK  <- c("Cmax", "AUC")
CV  <- c(0.45, 0.30)
# extract sample sizes and power
r1  <- sampleN.scABEL(CV = CV[1], design = "2x2x4",
                      print = FALSE, details = FALSE)[8:9]
r2  <- sampleN.TOST(CV = CV[2], theta0 = 0.90, design = "2x2x4",
                    print = FALSE, details = FALSE)[7:8]
n   <- as.numeric(c(r1[1], r2[1]))
pwr <- signif(as.numeric(c(r1[2], r2[2])), 5)
# compile results
res <- data.frame(PK = PK, method = c("ABEL", "ABE"),
                  n = n, power = pwr)
print(res, row.names = FALSE)
#    PK method  n   power
#  Cmax   ABEL 28 0.81116
#   AUC    ABE 40 0.80999
```

*AUC* drives the sample size.

For Health Canada it is the opposite (ABE for *C*<sub>max</sub> and ABEL
for *AUC*).

``` r
PK  <- c("Cmax", "AUC")
CV  <- c(0.45, 0.30)
# extract sample sizes and power
r1  <- sampleN.TOST(CV = CV[1], theta0 = 0.90, design = "2x2x4",
                    print = FALSE, details = FALSE)[7:8]
r2  <- sampleN.scABEL(CV = CV[2], design = "2x2x4",
                      print = FALSE, details = FALSE)[8:9]
n   <- as.numeric(c(r1[1], r2[1]))
pwr <- signif(as.numeric(c(r1[2], r2[2])), 5)
# compile results
res <- data.frame(PK = PK, method = c("ABE", "ABEL"),
                  n = n, power = pwr)
print(res, row.names = FALSE)
#    PK method  n   power
#  Cmax    ABE 84 0.80569
#   AUC   ABEL 34 0.80281
```

Here *C*<sub>max</sub> drives the sample size.

<small>[TOC ↩](#powertost)</small>

Sample size assuming homoscedasticity (*CV*<sub>wT</sub> =
*CV*<sub>wR</sub> = 0.45) for the widened limits of the Gulf Cooperation
Council.

``` r
sampleN.scABEL(CV = 0.45, regulator = "GCC", details = FALSE)
# 
# +++++++++++ scaled (widened) ABEL +++++++++++
#             Sample size estimation
#    (simulation based on ANOVA evaluation)
# ---------------------------------------------
# Study design: 2x3x3 (partial replicate)
# log-transformed data (multiplicative model)
# 1e+05 studies for each step simulated.
# 
# alpha  = 0.05, target power = 0.8
# CVw(T) = 0.45; CVw(R) = 0.45
# True ratio = 0.9
# ABE limits / PE constraint = 0.8 ... 1.25 
# Widened limits = 0.75 ... 1.333333 
# Regulatory settings: GCC 
# 
# Sample size
#  n     power
# 54   0.8123
```

<small>[TOC ↩](#powertost)</small>

#### RSABE

#### HVD(P)s

Sample size for a four period full replicate study (any of TRTR\|RTRT,
TRRT\|RTTR, TTRR\|RRTT) assuming heteroscedasticity (*CV*<sub>wT</sub>
0.40, *CV*<sub>wR</sub> 0.50, *i.e.*, variance-ratio ~0.67). Details of
the sample size search suppressed.

``` r
sampleN.RSABE(CV = c(0.40, 0.50), design = "2x2x4", details = FALSE)
# 
# ++++++++ Reference scaled ABE crit. +++++++++
#            Sample size estimation
# ---------------------------------------------
# Study design: 2x2x4 (4 period full replicate)
# log-transformed data (multiplicative model)
# 1e+05 studies for each step simulated.
# 
# alpha  = 0.05, target power = 0.8
# CVw(T) = 0.4; CVw(R) = 0.5
# True ratio = 0.9
# ABE limits / PE constraints = 0.8 ... 1.25 
# Regulatory settings: FDA 
# 
# Sample size
#  n    power
# 20   0.81509
```

<small>[TOC ↩](#powertost)</small>

#### NTIDs (FDA, CDE)

Sample size assuming heteroscedasticity (*CV*<sub>w</sub> 0.10,
variance-ratio 2.5, *i.e.*, the test treatment has a substantially
higher variability than the reference). TRTR\|RTRT according to the
FDA’s guidance.<sup id="a7">[7](#f7)</sup> Assess additionally which one
of the three components (scaled
<span title="Average Bioequivalence">ABE</span>, conventional
<span title="Average Bioequivalence">ABE</span>,
*s*<sub>wT</sub>/*s*<sub>wR</sub>-ratio) drives the sample size.

``` r
CV <- signif(CVp2CV(CV = 0.10, ratio = 2.5), 4)
n  <- sampleN.NTID(CV = CV)[["Sample size"]]
# 
# +++++++++++ FDA method for NTIDs ++++++++++++
#            Sample size estimation
# ---------------------------------------------
# Study design:  2x2x4 (TRTR|RTRT) 
# log-transformed data (multiplicative model)
# 1e+05 studies for each step simulated.
# 
# alpha  = 0.05, target power = 0.8
# CVw(T) = 0.1197, CVw(R) = 0.07551
# True ratio     = 0.975 
# ABE limits     = 0.8 ... 1.25 
# Implied scABEL = 0.9236 ... 1.0827 
# Regulatory settings: FDA 
# - Regulatory const. = 1.053605 
# - 'CVcap'           = 0.2142 
# 
# Sample size search
#  n     power
# 32   0.699120 
# 34   0.730910 
# 36   0.761440 
# 38   0.785910 
# 40   0.809580
suppressMessages(power.NTID(CV = CV, n = n, details = TRUE))
#        p(BE)  p(BE-sABEc)    p(BE-ABE) p(BE-sratio) 
#      0.80958      0.90966      1.00000      0.87447
```

The *s*<sub>wT</sub>/*s*<sub>wR</sub> component shows the lowest
probability to demonstrate <span title="Bioequivalence">BE</span> and
hence, drives the sample size.

<small>[TOC ↩](#powertost)</small>

Compare that with homoscedasticity (*CV*<sub>wT</sub> =
*CV*<sub>wR</sub> = 0.10):

``` r
CV <- 0.10
n  <- sampleN.NTID(CV = CV, details = FALSE)[["Sample size"]]
# 
# +++++++++++ FDA method for NTIDs ++++++++++++
#            Sample size estimation
# ---------------------------------------------
# Study design:  2x2x4 (TRTR|RTRT) 
# log-transformed data (multiplicative model)
# 1e+05 studies for each step simulated.
# 
# alpha  = 0.05, target power = 0.8
# CVw(T) = 0.1, CVw(R) = 0.1
# True ratio     = 0.975 
# ABE limits     = 0.8 ... 1.25 
# Regulatory settings: FDA 
# 
# Sample size
#  n     power
# 18   0.841790
suppressMessages(power.NTID(CV = CV, n = n, details = TRUE))
#        p(BE)  p(BE-sABEc)    p(BE-ABE) p(BE-sratio) 
#      0.84179      0.85628      1.00000      0.97210
```

Here the scaled <span title="Average Bioequivalence">ABE</span>
component shows the lowest probability to demonstrate
<span title="Bioequivalence">BE</span> and drives the sample size –
which is much lower than in the previous example.

<small>[TOC ↩](#powertost)</small>

Comparison with *fixed* narrower limits applicable in other
jurisdictions. Note that a replicate design is not mandatory – reducing
the chance of dropouts and requiring less administrations

``` r
CV  <- 0.10
# extract sample sizes and power
r1  <- sampleN.NTID(CV = CV, print = FALSE, details = FALSE)[8:9]
r2  <- sampleN.TOST(CV = CV, theta0 = 0.975, theta1 = 0.90,
                    design = "2x2x4", print = FALSE, details = FALSE)[7:8]
r3  <- sampleN.TOST(CV = CV, theta0 = 0.975, theta1 = 0.90,
                    design = "2x2x3", print = FALSE, details = FALSE)[7:8]
r4  <- sampleN.TOST(CV = CV, theta0 = 0.975, theta1 = 0.90,
                    print = FALSE, details = FALSE)[7:8]
n   <- as.numeric(c(r1[1], r2[1], r3[1], r4[1]))
pwr <- signif(as.numeric(c(r1[2], r2[2], r3[2], r4[2])), 5)
# compile results
res <- data.frame(method = c("FDA/CDE", rep ("fixed narrow", 3)),
                  design = c(rep("2x2x4", 2), "2x2x3", "2x2x2"),
                  n = n, power = pwr, a = n * c(4, 4, 3, 2))
names(res)[5] <- "adm. #" # number of administrations
print(res, row.names = FALSE)
#        method design  n   power adm. #
#       FDA/CDE  2x2x4 18 0.84179     72
#  fixed narrow  2x2x4 12 0.85628     48
#  fixed narrow  2x2x3 16 0.81393     48
#  fixed narrow  2x2x2 22 0.81702     44
```

<small>[TOC ↩](#powertost)</small>

### Dose-Proportionality

*CV* 0.20 (20%), doses 1, 2, and 8 units, assumed slope *β*<sub>0</sub>
1, target power 0.90.

``` r
sampleN.dp(CV = 0.20, doses = c(1, 2, 8), beta0 = 1, targetpower = 0.90)
# 
# ++++ Dose proportionality study, power model ++++
#             Sample size estimation
# -------------------------------------------------
# Study design: crossover (3x3 Latin square) 
# alpha = 0.05, target power = 0.9
# Equivalence margins of R(dnm) = 0.8 ... 1.25 
# Doses = 1 2 8 
# True slope = 1, CV = 0.2
# Slope acceptance range = 0.89269 ... 1.1073 
# 
# Sample size (total)
#  n     power
# 18   0.915574
```

Note that the acceptance range of the slope depends on the ratio of the
highest and lowest doses (*i.e.*, it gets tighter for wider dose ranges
and therefore, higher sample sizes will be required).  
In an exploratory setting wider equivalence margins {*θ*<sub>1</sub>,
*θ*<sub>2</sub>} (0.50, 2.00) were proposed,<sup id="a8">[8](#f8)</sup>
translating in this example to an acceptance range of
`0.66667 ... 1.3333` and a sample size of only six subjects.

<small>[TOC ↩](#powertost)</small>

### Power Analysis

Explore impact of deviations from assumptions (higher *CV*, higher
deviation of *θ*<sub>0</sub> from 1, dropouts) on power. Assumed
within-subject *CV* 0.20 (20%), target power 0.90. Plot suppressed.

``` r
res <- pa.ABE(CV = 0.20, targetpower = 0.90)
print(res, plotit = FALSE)
# Sample size plan ABE
#  Design alpha  CV theta0 theta1 theta2 Sample size Achieved power
#     2x2  0.05 0.2   0.95    0.8   1.25          26      0.9176333
# 
# Power analysis
# CV, theta0 and number of subjects leading to min. acceptable power of ~0.7:
#  CV= 0.2729, theta0= 0.9044
#  n = 16 (power= 0.7354)
```

If the study starts with 26 subjects (power ~0.92), the *CV* can
increase to ~0.27 **or** *θ*<sub>0</sub> decrease to ~0.90 **or** the
sample size decrease to 10 whilst power will still be ≥0.70.  
However, this is **not** a substitute for the ‘Sensitivity Analysis’
recommended in ICH-E9,<sup id="a9">[9](#f9)</sup> since in a real study
a combination of all effects occurs simultaneously. It is up to *you* to
decide on reasonable combinations and analyze their respective power.

<small>[TOC ↩](#powertost)</small>

### Speed Comparisons

Performed on a Xeon E3-1245v3 3.4 GHz, 8 MB cache, 16 GB RAM, R 4.3.3
64 bit on Windows 7.

#### ABE

2×2 crossover design, *CV* 0.17. Sample sizes and achieved power for the
supported methods (the 1<sup>st</sup> one is the default).

        method  n   power time (s)
         owenq 14 0.80568  0.00128
           mvt 14 0.80569  0.11778
    noncentral 14 0.80568  0.00100
       shifted 16 0.85230  0.00096

The 2<sup>nd</sup> exact method is substantially slower than the
1<sup>st</sup>. The approximation based on the noncentral
*t*-distribution is slightly faster but matches the 1<sup>st</sup> exact
method closely. Though the approximation based on the shifted central
*t*-distribution is the fastest, it *might* estimate a larger than
necessary sample size. Hence, it should be used only for comparative
purposes.

#### ABEL

Four period full replicate study (any of TRTR\|RTRT, TRRT\|RTTR,
TTRR\|RRTT), homogenicity (*CV*<sub>wT</sub> = *CV*<sub>wR</sub> 0.45).
Sample sizes and achieved power for the supported methods.

                  function              method  n   power time (s)
            sampleN.scABEL    ‘key’ statistics 28 0.81116   0.1348
     sampleN.scABEL.sdsims subject simulations 28 0.81196   2.5377

Simulating via the ‘key’ statistics is the method of choice for speed
reasons.  
However, subject simulations are recommended **if**

- the partial replicate design (TRR\|RTR\|RRT) is planned **and**
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
Version 1.5.6 built 2024-03-18 with R 4.3.3.

``` r
options(width = 66)
sessionInfo()
# R version 4.3.3 (2024-02-29 ucrt)
# Platform: x86_64-w64-mingw32/x64 (64-bit)
# Running under: Windows 10 x64 (build 19045)
# 
# Matrix products: default
# 
# 
# locale:
# [1] LC_COLLATE=German_Germany.utf8 
# [2] LC_CTYPE=German_Germany.utf8   
# [3] LC_MONETARY=German_Germany.utf8
# [4] LC_NUMERIC=C                   
# [5] LC_TIME=German_Germany.utf8    
# 
# time zone: Europe/Berlin
# tzcode source: internal
# 
# attached base packages:
# [1] stats     graphics  grDevices utils     datasets  methods  
# [7] base     
# 
# other attached packages:
# [1] PowerTOST_1.5-6
# 
# loaded via a namespace (and not attached):
#  [1] cubature_2.1.0    compiler_4.3.3    fastmap_1.1.1    
#  [4] cli_3.6.2         tools_4.3.3       htmltools_0.5.7  
#  [7] rstudioapi_0.15.0 yaml_2.3.8        Rcpp_1.0.12      
# [10] mvtnorm_1.2-4     rmarkdown_2.26    knitr_1.45       
# [13] xfun_0.42         digest_0.6.35     rlang_1.1.3      
# [16] evaluate_0.23
```

<small>[TOC ↩](#powertost)</small>

------------------------------------------------------------------------

<span id="f1" style="font-size:smaller">1. Tóthfalusi L, Endrényi L.
*Sample Sizes for Designing Bioequivalence Studies for Highly Variable
Drugs.* J Pharm Pharmacol Sci. 2012; 15(1): 73–84.
[doi:10.18433/j3z88f](https://doi.org/10.18433/j3z88f). [Open
access](https://ejournals.library.ualberta.ca/index.php/JPPS/article/download/11612/9489).</span>
[↩](#a1)  
<span id="f2" style="font-size:smaller">2. Fieller EC. *Some Problems In
Interval Estimation.* J Royal Stat Soc B. 1954; 16(2): 175–85.
[JSTOR:2984043](https://www.jstor.org/stable/2984043).</span> [↩](#a2)  
<span id="f3" style="font-size:smaller">3. U.S. Food and Drug
Administration, Office of Generic Drugs. *Draft Guidance on Dabigatran
Etexilate Mesylate.* Recommended Jun 2012; Revised Sep 2015, Jul 2017.
[Online](https://www.accessdata.fda.gov/drugsatfda_docs/psg/Dabigatran%20etexilate%20mesylate_oral%20capsule_NDA%20022512_RV05-17.pdf).</span>
[↩](#a3)  
<span id="f4" style="font-size:smaller">4. U.S. Food and Drug
Administration, Office of Generic Drugs. *Draft Guidance on
Rivaroxaban.* Recommended Sep 2015.
[Online](https://www.accessdata.fda.gov/drugsatfda_docs/psg/Rivaroxaban_oral%20tablet_22406_RC09-15.pdf).</span>
[↩](#a4)  
<span id="f5" style="font-size:smaller">5. U.S. Food and Drug
Administration, Office of Generic Drugs. *Draft Guidance on Edoxaban
Tosylate.* Recommended May 2017; Revised Mar 2020.
[Online](https://www.accessdata.fda.gov/drugsatfda_docs/psg/PSG_206316.pdf).</span>
[↩](#a5)  
<span id="f6" style="font-size:smaller">6. Labes D, Schütz H. *Inflation
of Type I Error in the Evaluation of Scaled Average Bioequivalence, and
a Method for its Control.* Pharm Res. 2016; 33(11): 2805–14.
[doi:10.1007/s11095-016-2006-1](https://doi.org/10.1007/s11095-016-2006-1).</span>
[↩](#a6)  
<span id="f7" style="font-size:smaller">7. U.S. Food and Drug
Administration, Center for Drug Evaluation and Research. *Draft Guidance
for Industry. Bioequivalence Studies with Pharmacokinetic Endpoints for
Drugs Submitted Under an ANDA.* August 2021.
[Online](https://www.fda.gov/media/87219/download).</span> [↩](#a7)  
<span id="f8" style="font-size:smaller">8. Hummel J, McKendrick S,
Brindley C, French R. *Exploratory assessment of dose proportionality:
review of current approaches and proposal for a practical criterion.*
Pharm. Stat. 2009; 8(1): 38–49.
[doi:10.1002/pst.326](https://doi.org/10.1002/pst.326).</span>
[↩](#a8)  
<span id="f9" style="font-size:smaller">9. International Conference on
Harmonisation of Technical Requirements for Registration of
Pharmaceuticals for Human Use. *ICH Harmonised Tripartite Guideline. E9.
Statistical Principles for Clinical Trials.* 5 February 1998.
[Online](https://database.ich.org/sites/default/files/E9_Guideline.pdf).</span>
[↩](#a9)
