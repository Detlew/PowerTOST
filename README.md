README
================

<!-- README.md is generated from README.Rmd. Please edit that file -->

## PowerTOST

The package contains functions to calculate power and sample size for
various study designs used in bioequivalence studies.

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
    reference (*θ*<sub>0</sub>), accteptance limits {*θ*<sub>1</sub>,
    *θ*<sub>2</sub>}, sample size (*n*), and design.

For all methods the sample size can be *estimated* based on

  - nominal *α*, coefficient of variation (*CV*), deviation of test from
    reference(*θ*<sub>0</sub>), accteptance limits {*θ*<sub>1</sub>,
    *θ*<sub>2</sub>}, target power, and design.

### Supported methods

#### Power and sample size

  - Average bioequivalence (with arbitrary *fixed* limits).
  - Reference-scaled bioequivalence (EMA: ABEL, FDA: RSABE for HVD(P)s
    and NTIDs).
  - Iteratively adjust *α* to control the type I error in ABEL and
    RSABE.
  - Two simulataneous TOST procedures.
  - Dose proportionality using the power model.
  - Non-inferiority *t*-test.
  - Ratio of two means with normally distributed data on the original
    scale based on Fieller’s (‘fiducial’) confidence interval.
  - ‘Expected’ power in case of uncertain (estimated) variability and/or
    uncertain *θ*<sub>0</sub>.

#### Helpers

  - Calculate *CV* from *MSE* or *SE* (and vice versa).
  - Calculate *CV* from given confidence interval.
  - Confidence interval of *CV*.
  - Pool *CV* from several studies.
  - Confidence interval for given *α*, *CV*, point estimate, sample
    size, and design.
  - *p*-values of the TOST procedure.
  - Analysis tool for exploration/visualization of the impact of
    expected values (*CV*, *θ*<sub>0</sub>, reduced sample size due to
    drop-outs) on power of BE decision.

### Examples

1.  Power of a paralled design (group sizes 75 and 70), total *CV* 0.36,
    *θ*<sub>0</sub> 0.925. The defaults for *α* (0.05) and
    {*θ*<sub>1</sub>, *θ*<sub>2</sub>} (0.80, 1.25) are automatically
    employed.

<!-- end list -->

``` r
PowerTOST:::power.TOST(CV = 0.36, theta0 = 0.93,
                       n = c(75, 70), design = "parallel")
#> [1] 0.8256362
```

2.  Sample size for a 2×2×2 crossover study, assumed intra-subject *CV*
    0.20, and desired power 0.80. The defaults for *α* (0.05),
    *θ*<sub>0</sub> (0.95), {*θ*<sub>1</sub>, *θ*<sub>2</sub>} (0.80,
    1.25), targetpower (0.80), and design (“2x2”) are automatically
    employed.

<!-- end list -->

``` r
PowerTOST:::sampleN.TOST(CV = 0.20)
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

3.  Sample size for a four-period replicate study intended for the EMA’s
    average bioequivalence with expanding limits (ABEL), assumed
    heteroscedasticity (*CV<sub>wT</sub>* 0.40, *CV<sub>wR</sub>* 0.50).
    The defaults for *α* (0.05), *θ*<sub>0</sub> (0.90 for highly
    variable drugs / drug products), and targetpower (0.80) are
    automatically employed. The expanded limits {*θ*<sub>1</sub>,
    *θ*<sub>2</sub>} are calculated based on *CV<sub>wR</sub>*.

<!-- end list -->

``` r
PowerTOST:::sampleN.scABEL(CV = c(0.40, 0.50), design = "2x2x4")
#> 
#> +++++++++++ scaled (widened) ABEL +++++++++++
#>             Sample size estimation
#>    (simulation based on ANOVA evaluation)
#> ---------------------------------------------
#> Study design:  2x2x4 (full replicate) 
#> log-transformed data (multiplicative model)
#> 1e+05 studies for each step simulated.
#> 
#> alpha  = 0.05, target power = 0.8
#> CVw(T) = 0.4; CVw(R) = 0.5
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
#> 18   0.6963 
#> 20   0.7427 
#> 22   0.7829 
#> 24   0.8172
```

### Installation

``` r
package <- "PowerTOST"
inst    <- package %in% installed.packages()
if (length(package[!inst]) > 0) install.packages(package[!inst])
```
