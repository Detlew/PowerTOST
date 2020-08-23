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
![active](https://www.repostatus.org/badges/latest/active.svg) ![on
CRAN](https://www.r-pkg.org/badges/version-ago/PowerTOST) [![cran
checks](https://cranchecks.info/badges/worst/PowerTOST)](https://cran.r-project.org/web/checks/check_results_PowerTOST.html)
[![CRAN RStudio mirror
downloads](https://cranlogs.r-pkg.org/badges/grand-total/PowerTOST?color=blue)](https://r-pkg.org/pkg/PowerTOST)
[![CRAN RStudio mirror
downloads](https://cranlogs.r-pkg.org/badges/last-month/PowerTOST?color=green)](https://r-pkg.org/pkg/PowerTOST)

The package contains functions to calculate power and estimate sample
size for various study designs used in (not only bio-) equivalence
studies.

    Version 1.5.0.9999 built 2020-08-23 with R 4.0.2 
    (development version not on CRAN).

## Supported Designs

``` 
   design                        name    df
 parallel           2 parallel groups   n-2
      2x2               2x2 crossover   n-2
    2x2x2             2x2x2 crossover   n-2
      3x3               3x3 crossover 2*n-4
    3x6x3             3x6x3 crossover 2*n-4
      4x4               4x4 crossover 3*n-6
    2x2x3   2x2x3 replicate crossover 2*n-3
    2x2x4   2x2x4 replicate crossover 3*n-4
    2x4x4   2x4x4 replicate crossover 3*n-4
    2x3x3   partial replicate (2x3x3) 2*n-3
    2x4x2            Balaam's (2x4x2)   n-2
   2x2x2r Liu's 2x2x2 repeated x-over 3*n-2
   paired                paired means   n-1
```

Codes of designs follow this pattern: `treatments x sequences x
periods`.

Although some replicate designs are more ‘popular’ than others, sample
size estimations are valid for *all* of the following designs:

| design  |  type   | sequences         | periods |
| :-----: | :-----: | ----------------- | :-----: |
| `2x2x4` |  full   | 2 `TRTR\|RTRT`    |    4    |
| `2x2x4` |  full   | 2 `TRRT\|RTTR`    |    4    |
| `2x2x4` |  full   | 2 `TTRR\|RRTT`    |    4    |
| `2x2x3` |  full   | 2 `TRT\|RTR`      |    3    |
| `2x2x3` |  full   | 2 `TRR\|RTT`      |    3    |
| `2x3x3` | partial | 3 `TRR\|RTR\|RRT` |    3    |

Whilst "2x4x4" four period full replicate designs with *four* sequences
(TRTR|RTRT|TRRT|RTTR *or* TRRT|RTTR|TTRR|RRTT) are supported, they
should be avoided due to confounded effects.

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
      - <span title="European Medicines Agency">EMA</span>: Average
        Bioequivalence with Expanding Limits (ABEL).  
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
[1] 0.8011186
```

### Crossover Design

Sample size for assumed within- (intra-) subject *CV* 0.20 (20%).

``` r
sampleN.TOST(CV = 0.20)

+++++++++++ Equivalence test - TOST +++++++++++
            Sample size estimation
-----------------------------------------------
Study design: 2x2 crossover 
log-transformed data (multiplicative model)

alpha = 0.05, target power = 0.8
BE margins = 0.8 ... 1.25 
True ratio = 0.95,  CV = 0.2

Sample size (total)
 n     power
20   0.834680 
```

Sample size for assumed within- (intra-) subject *CV* 0.40 (40%),
*θ*<sub>0</sub> 0.90, four period full replicate “2x2x4” study. Wider
acceptance range for *C*<sub>max</sub> (Gulf Cooperation Council, South
Africa).

``` r
sampleN.TOST(CV = 0.40, theta0 = 0.90, theta1 = 0.75, design = "2x2x4")

+++++++++++ Equivalence test - TOST +++++++++++
            Sample size estimation
-----------------------------------------------
Study design: 2x2x4 (4 period full replicate) 
log-transformed data (multiplicative model)

alpha = 0.05, target power = 0.8
BE margins = 0.75 ... 1.333333 
True ratio = 0.9,  CV = 0.4

Sample size (total)
 n     power
30   0.822929 
```

<small>[TOC ↩](#powertost)</small>

Sample size for assumed within- (intra-) subject *CV* 0.125 (12.5%),
*θ*<sub>0</sub> 0.975. Acceptance range for
<span title="Narrow Therapeutic Index Drugs">NTIDs</span> (most
jurisdictions).

``` r
sampleN.TOST(CV = 0.125, theta0 = 0.975, theta1 = 0.90)

+++++++++++ Equivalence test - TOST +++++++++++
            Sample size estimation
-----------------------------------------------
Study design: 2x2 crossover 
log-transformed data (multiplicative model)

alpha = 0.05, target power = 0.8
BE margins = 0.9 ... 1.111111 
True ratio = 0.975,  CV = 0.125

Sample size (total)
 n     power
32   0.800218 
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

+++++++++++ Equivalence test - TOST +++++++++++
    based on Fieller's confidence interval
            Sample size estimation
-----------------------------------------------
Study design: 2x2 crossover
Ratio of means with normality on original scale
alpha = 0.025, target power = 0.8
BE margins = 0.8 ... 1.25 
True ratio = 0.95,  CVw = 0.2,  CVb = 0.4

Sample size
 n     power
28   0.807774 
```

<small>[TOC ↩](#powertost)</small>

### Replicate Designs

#### ABE

Sample size for assumed within- (intra-) subject *CV* 0.45 (45%),
*θ*<sub>0</sub> 0.90, three period full replicate study "2x2x3"
(TRT|RTR *or* TRR|RTT).

``` r
sampleN.TOST(CV = 0.45, theta0 = 0.90, design = "2x2x3")

+++++++++++ Equivalence test - TOST +++++++++++
            Sample size estimation
-----------------------------------------------
Study design: 2x2x3 (3 period full replicate) 
log-transformed data (multiplicative model)

alpha = 0.05, target power = 0.8
BE margins = 0.8 ... 1.25 
True ratio = 0.9,  CV = 0.45

Sample size (total)
 n     power
124   0.800125 
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
   CVwT    CVwR Sample size Achieved power
 0.3987 0.49767         126         0.8052
```

Similar sample size because the pooled *CV* is still 0.45.

<small>[TOC ↩](#powertost)</small>

#### ABEL

Sample size assuming homoscedasticity (*CV*<sub>wT</sub> =
*CV*<sub>wR</sub> = 0.45).

``` r
sampleN.scABEL(CV = 0.45)

+++++++++++ scaled (widened) ABEL +++++++++++
            Sample size estimation
   (simulation based on ANOVA evaluation)
---------------------------------------------
Study design: 2x3x3 (partial replicate)
log-transformed data (multiplicative model)
1e+05 studies for each step simulated.

alpha  = 0.05, target power = 0.8
CVw(T) = 0.45; CVw(R) = 0.45
True ratio = 0.9
ABE limits / PE constraint = 0.8 ... 1.25 
EMA regulatory settings
- CVswitch            = 0.3 
- cap on scABEL if CVw(R) > 0.5
- regulatory constant = 0.76 
- pe constraint applied


Sample size search
 n     power
36   0.7755 
39   0.8059 
```

<small>[TOC ↩](#powertost)</small>

Iteratively adjust *α* to control the Type I Error ([Labes,
Schütz](https://doi.org/10.1007/s11095-016-2006-1)). Slight
heteroscedasticity (*CV*<sub>wT</sub> 0.30, *CV*<sub>wR</sub> 0.35),
four period full replicate "2x2x4" study, 30 subjects, balanced
sequences.

``` r
scABEL.ad(CV = c(0.30, 0.35), design = "2x2x4", n = 30)
+++++++++++ scaled (widened) ABEL ++++++++++++
         iteratively adjusted alpha
   (simulations based on ANOVA evaluation)
----------------------------------------------
Study design: 2x2x4 (4 period full replicate)
log-transformed data (multiplicative model)
1,000,000 studies in each iteration simulated.

CVwR 0.35, CVwT 0.3, n(i) 15|15 (N 30)
Nominal alpha                 : 0.05 
True ratio                    : 0.9000 
Regulatory settings           : EMA (ABEL)
Switching CVwR                : 0.3 
Regulatory constant           : 0.76 
Expanded limits               : 0.7723 ... 1.2948
Upper scaling cap             : CVwR > 0.5 
PE constraints                : 0.8000 ... 1.2500
Empiric TIE for alpha 0.0500  : 0.06651
Power for theta0 0.9000       : 0.814
Iteratively adjusted alpha    : 0.03540
Empiric TIE for adjusted alpha: 0.05000
Power for theta0 0.9000       : 0.771 
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
   PK method  n   power
 Cmax   ABEL 42 0.80017
  AUC    ABE 60 0.81002
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

++++++++ Reference scaled ABE crit. +++++++++
           Sample size estimation
---------------------------------------------
Study design: 2x2x4 (4 period full replicate)
log-transformed data (multiplicative model)
1e+05 studies for each step simulated.

alpha  = 0.05, target power = 0.8
CVw(T) = 0.4; CVw(R) = 0.5
True ratio = 0.9
ABE limits / PE constraints = 0.8 ... 1.25 
Regulatory settings: FDA 

Sample size
 n    power
20   0.81509 
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

+++++++++++ FDA method for NTIDs ++++++++++++
           Sample size estimation
---------------------------------------------
Study design:  2x2x4 (TRTR|RTRT) 
log-transformed data (multiplicative model)
1e+05 studies for each step simulated.

alpha  = 0.05, target power = 0.8
CVw(T) = 0.1197, CVw(R) = 0.07551
True ratio     = 0.975 
ABE limits     = 0.8 ... 1.25 
Implied scABEL = 0.9236 ... 1.0827 
Regulatory settings: FDA 
- Regulatory const. = 1.053605 
- 'CVcap'           = 0.2142 

Sample size search
 n     power
32   0.699120 
34   0.730910 
36   0.761440 
38   0.785910 
40   0.809580 
suppressMessages(power.NTIDFDA(CV = CV, n = n, details = TRUE))
       p(BE)  p(BE-sABEc)    p(BE-ABE) p(BE-sratio) 
     0.80958      0.90966      1.00000      0.87447 
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

+++++++++++ FDA method for NTIDs ++++++++++++
           Sample size estimation
---------------------------------------------
Study design:  2x2x4 (TRTR|RTRT) 
log-transformed data (multiplicative model)
1e+05 studies for each step simulated.

alpha  = 0.05, target power = 0.8
CVw(T) = 0.1, CVw(R) = 0.1
True ratio     = 0.975 
ABE limits     = 0.8 ... 1.25 
Regulatory settings: FDA 

Sample size
 n     power
18   0.841790 
suppressMessages(power.NTIDFDA(CV = CV, n = n, details = TRUE))
       p(BE)  p(BE-sABEc)    p(BE-ABE) p(BE-sratio) 
     0.84179      0.85628      1.00000      0.97210 
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
       method design  n   power adm. #
   FDA scaled  2x2x4 18 0.84179     72
 fixed narrow  2x2x4 12 0.85628     48
 fixed narrow  2x2x3 16 0.81393     48
 fixed narrow  2x2x2 22 0.81702     44
```

<small>[TOC ↩](#powertost)</small>

### Dose-Proportionality

*CV* 0.20 (20%), doses 1, 2, and 8 units, assumed slope *β*<sub>0</sub>
1, target power 0.90.

``` r
sampleN.dp(CV = 0.20, doses = c(1, 2, 8), beta0 = 1, targetpower = 0.90)

++++ Dose proportionality study, power model ++++
            Sample size estimation
-------------------------------------------------
Study design: crossover (3x3 Latin square) 
alpha = 0.05, target power = 0.9
Equivalence margins of R(dnm) = 0.8 ... 1.25 
Doses = 1 2 8 
True slope = 1, CV = 0.2
Slope acceptance range = 0.89269 ... 1.1073 

Sample size (total)
 n     power
18   0.915574 
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
Sample size plan ABE
 Design alpha  CV theta0 theta1 theta2 Sample size Achieved power
    2x2  0.05 0.2   0.95    0.8   1.25          26      0.9176333

Power analysis
CV, theta0 and number of subjects which lead to min. acceptable power of at least 0.7:
 CV= 0.2729, theta0= 0.9044
 n = 16 (power= 0.7354)
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

Performed on a Xeon E3-1245v3 3.4 GHz, 8 MB cache, 16 GB RAM, R 4.0.2
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
Version 1.5.0.9999 built 2020-08-23 with R 4.0.2.

``` r
options(width = 80)
devtools::session_info()
- Session info ---------------------------------------------------------------
 setting  value                       
 version  R version 4.0.2 (2020-06-22)
 os       Windows 10 x64              
 system   x86_64, mingw32             
 ui       RTerm                       
 language EN                          
 collate  German_Germany.1252         
 ctype    German_Germany.1252         
 tz       Europe/Berlin               
 date     2020-08-23                  

- Packages -------------------------------------------------------------------
 package       * version    date       lib source        
 assertthat      0.2.1      2019-03-21 [1] CRAN (R 4.0.0)
 backports       1.1.7      2020-05-13 [1] CRAN (R 4.0.0)
 callr           3.4.3      2020-03-28 [1] CRAN (R 4.0.0)
 cli             2.0.2      2020-02-28 [1] CRAN (R 4.0.0)
 crayon          1.3.4      2017-09-16 [1] CRAN (R 4.0.0)
 cubature        2.0.4.1    2020-07-06 [1] CRAN (R 4.0.2)
 desc            1.2.0      2018-05-01 [1] CRAN (R 4.0.0)
 devtools        2.3.1      2020-07-21 [1] CRAN (R 4.0.2)
 digest          0.6.25     2020-02-23 [1] CRAN (R 4.0.0)
 ellipsis        0.3.1      2020-05-15 [1] CRAN (R 4.0.0)
 evaluate        0.14       2019-05-28 [1] CRAN (R 4.0.0)
 fansi           0.4.1      2020-01-08 [1] CRAN (R 4.0.0)
 fs              1.5.0      2020-07-31 [1] CRAN (R 4.0.2)
 glue            1.4.1      2020-05-13 [1] CRAN (R 4.0.0)
 htmltools       0.5.0      2020-06-16 [1] CRAN (R 4.0.0)
 knitr           1.29       2020-06-23 [1] CRAN (R 4.0.2)
 magrittr        1.5        2014-11-22 [1] CRAN (R 4.0.0)
 memoise         1.1.0      2017-04-21 [1] CRAN (R 4.0.0)
 mvtnorm         1.1-1      2020-06-09 [1] CRAN (R 4.0.0)
 pkgbuild        1.1.0      2020-07-13 [1] CRAN (R 4.0.2)
 pkgload         1.1.0      2020-05-29 [1] CRAN (R 4.0.0)
 PowerTOST     * 1.5.0.9999 2020-08-23 [1] local         
 prettyunits     1.1.1      2020-01-24 [1] CRAN (R 4.0.0)
 processx        3.4.3      2020-07-05 [1] CRAN (R 4.0.2)
 ps              1.3.4      2020-08-11 [1] CRAN (R 4.0.2)
 R6              2.4.1      2019-11-12 [1] CRAN (R 4.0.0)
 Rcpp            1.0.5      2020-07-06 [1] CRAN (R 4.0.2)
 remotes         2.2.0      2020-07-21 [1] CRAN (R 4.0.2)
 rlang           0.4.7      2020-07-09 [1] CRAN (R 4.0.2)
 rmarkdown       2.3        2020-06-18 [1] CRAN (R 4.0.0)
 rprojroot       1.3-2      2018-01-03 [1] CRAN (R 4.0.0)
 sessioninfo     1.1.1      2018-11-05 [1] CRAN (R 4.0.0)
 stringi         1.4.6      2020-02-17 [1] CRAN (R 4.0.0)
 stringr         1.4.0      2019-02-10 [1] CRAN (R 4.0.0)
 TeachingDemos   2.12       2020-04-07 [1] CRAN (R 4.0.0)
 testthat        2.3.2      2020-03-02 [1] CRAN (R 4.0.0)
 usethis         1.6.1      2020-04-29 [1] CRAN (R 4.0.0)
 withr           2.2.0      2020-04-20 [1] CRAN (R 4.0.0)
 xfun            0.16       2020-07-24 [1] CRAN (R 4.0.2)
 yaml            2.2.1      2020-02-01 [1] CRAN (R 4.0.0)

[1] C:/Program Files/R/library
[2] C:/Program Files/R/R-4.0.2/library
```

<small>[TOC ↩](#powertost)</small>
