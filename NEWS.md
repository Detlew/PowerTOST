# PowerTOST 1.5-6

On CRAN submitted 2024-03-18.

## Bug fixes
  * Bug fixed in `CI.BE()` regarding cases with only `Ntotal` given.
  
## Minor changes
  * Function `shadowtext()` lifted from package `TeachingDemos` and recoded since
    `TeachingDemos`was orphaned (request of Brian Ripley). Orphaned status later 
    on abandoned but nevertheless the recoded function used furthermore. 
  * Updated the RSABE-vignette and the Man page of `sampleN.scABEL()` to reflect 
    the GCC guidelines (Version 3.1 of 10 August 2022).
  * Example results in man pages of functions `power.2TOST()` and 
    `SampleN.2TOST()` corrected.

## Major changes
  * Argument `p.level` in function `power.TOST.sds()` introduced to make the significance level
    of the group-by-treatment interaction test variable. Was former hard coded to 0.1.
  * The name of functions dealing with the FDA methods for NTID are now without the acronym FDA
    in their names since the evaluation methods requested by the FDA are also required by China CDE. 
    These functions are `power.NTID()`, `sampleN.NTID()`, and `pa.NTID()`.  
    See NEWS under `PowerTOST 1.5-4`.

# PowerTOST 1.5-5

Never on CRAN, only a development version.

# PowerTOST 1.5-4

On CRAN 2022-02-21.

## Bug fixes
  * Incomplete HTML-entities in 3 man pages corrected (request of Kurt Hornik).
  * Fix in `power.TOST()` w.r.t vectorization of `CV` or `theta0`.

## Major changes
  * Functions dealing with the FDA method for NTID have now aliases without FDA in their names since the evaluation methods requested by the FDA are also required by China CDE.
    The aliases are `power.NTID()`, `sampleN.NTID()`, and `pa.NTID()`.
    The functions `power.NTIDFDA()`, `sampleN.NTIDFDA()`, and `pa.NTIDFDA()` are deprecated and will be removed in the next release.
    
## Minor changes
  * Example with vectorized `sampleN.TOST()` added to vignette ABE; vignette PA adapted to `pa.NTID()`. 
  * URLs in README and man pages updated.
  * Coefficients with more decimals (15) for the 10-point Gauss quadrature in function `tfn()`. Thanks to PharmCat for contributing them.
  * Clarification of the calculations with `gmodel = 1` in the man page section Details of function `power.TOST.sds()`.

# PowerTOST 1.5-3

On CRAN 2021-01-18

## Major changes

  * `scABEL.ad()` and `sampleN.scABEL.ad()` adapted to deal with `regulator = "GCC"`. 
    Man pages, README, and RSABE-vignette updated accordingly.
  * Regulator "GCC" introduced (GCC = Co-operation Council for the Arab States of the Gulf). The GCC evaluation framework for HVD / HVDP is/can be
    treated as a special case of ABEL, *i.e.*, use the regulatory settings with `power.scABEL()` or `sampleN.scABEL()`.

# PowerTOST 1.5-2

On CRAN 2020-10-27.

## Bug fixes
  * `stringsAsFactors = FALSE` in example of main vignette for R <4.0.0.

# PowerTOST 1.5-1

On CRAN 2020-10-22.

## Bug fixes
  * Check of arg `theta0` within range `theta1 ... theta2` fixed.
  * Fix of the default settings of `theta1, theta2` if missing in various functions.
  * Added `emmeans` to `Suggests`. Otherwise `NOTE` on r-devel-linux-x86_64-fedora-clang. 
  THX to Gábor Csárdi and Duncan Murdoch @r-pkg-devel. Will be required in R4.1.0 released next year.

## Minor changes

  * Added example for statistical assurance to the ABE-Vignette. 
  * More checks of CV and theta0 vectorized input to power.TOST(). Only one argument may be a vector.
  * Unified man pages of `sampleN.*`: Always *estimation* instead of *calculation*. References added for Fieller’s CI.
  * More examples in `README`.
  * Changed http(s) locations of References according to R CMD check.
  * Cosmetics in Vignettes. Added examples for `logscale = FALSE` to the ABE-Vignette.

# PowerTOST 1.5-0 

On CRAN 2020-08-09. (Maintenance release dedicated to 70 birthday of Detluuu)

## Bug fixes
  
  * Bug fix in `scABEL.ad` w.r.t. coercion.
  * Fix of the link to package `emmeans` in man pages of expected power.
  * Fix in Example 1 of `ABE.Rmd` (not a good idea to assign a variable with the same name as a function).

## Major changes

  * Pre-compiled RSABE-Vignette.

## Minor changes

  * Highlight clarification about *total* sample size in vignettes.
  * Removed links to man pages in vignettes (they work in the library but not in the public ones on CRAN).
  * Moved `tufte` from Imports to Suggests (Duncan Murdoch @r-pkg-devel).
  * Harmonize default value of `theta1` in `pvalue.TOST` in case of `logscale = FALSE`
  * Clarification of argument `CV` (and `theta0`, `theta1`, `theta2`) in case of `logscale = FALSE`.
  * Clarify in man pages and vignettes that all functions return the *total* sample size (and not subjects/sequence in crossovers and subjects/group in parallel designs -- like in some other software packages). Suggested by Amandine Schmutz.

# PowerTOST 1.4-9

On CRAN 2019-12-19. (Xmas gift)

## Bug fixes

  * Fix in `power.dp()` w.r.t. setting argument `CVb` if missing in case of `design="IBD"`.
  * Fix in `scABEL.ad()`: `reg$name` instead of `regulator`.
  * CV in 4th panel of `pwrA_S3methods.R` to the same precison like in the others.
  * Fix in `power.TOST.sds()` for `gmodel=1`, case of `gmodel=3` with data of the largest group (group by treatment interaction significant at p=0.1)
  
## Major changes

  * (Pre-compiled) vignettes.

## Minor changes

  * Imports package `tufte` for nice quotes in Rmarkdown.
  * Slightly enhanced man pages of `power.dp()` and `sampleN.dp()` w.r.t. the value of `CVb` in case of `design="IBD"`.
  * Cosmetics in output of `sampleN.noninf` based on `margin`.
  * Reworked minimum sample sizes in `pa.scABE.R()` according to guidances. Generally 12 (as before) but 24 for the FDA. Also 24 for the EMA if `2x2x3` design (Q&A document). Changed `N` to `n` in the S3-methods in conformity with other functions.

# PowerTOST 1.4-8

On CRAN 2019-08-29.

## Bug fixes

  * `scABEL.ad()`, `sampleN.scABEL.ad()`: `CVwT` was not given in output if `CV`was a vector.
  * In `scABEL.ad.R` `regulator`was `NULL`. Corrected to `reg$name`.
  * Smaller font in 4<sup>th</sup> screen of `pwrA_S3methods.R` of `pa.NTID()` (more lines required). Plain font instead of bold in main.
  * Broken FDA URL in many man pages corrected.
  * In `expsampleN.noninf()` wrt operator `&&` with vector arguments (new check in R 3.6.0).
  * `CI.RatioF()` fixed. Thanks to Michael (mittyri).
  
## Major changes

  * New function `sampleN.RSABE2L.sdsims()` for sample size estimation via subject simulations for the ‘exact’ method of Tóthfalusi & Endrényi “*[Algorithms for evaluating reference scaled average bioequivalence: power, bias, and consumer risk](https://doi.org/10.1002/sim.7440)*”.
  * New function `sampleN.scABEL.sdsims()` for sample size estimation for ABEL via subject simulations. Useful in case of assumed heteroscedasticity in the partial replicate design (TRT|RTR|RRT).
  * `NEWS.md` instead of `NEWS`.
  * `README.md` (knitted from `README.Rmd`).
  * Function `power.TOST.sds()` to simulate subject data & evaluate via models with group effect implemented. 
  * Unified Code base for the scaled ABE functions.

## Minor changes

  * LaTeX-builder on CRAN laments about UTF-8 characters in man-pages. No problem rendering the PDF-manual locally! Replaced all UTF-8 characters by `\enc{foo}{bar}`.
  * In `power.scABEL()`: If `nsims` not given, defaults to 1e5 (like before). If `theta0` equals one of the expanded limits, defaults to 1e6. Man-page updated.
  * Output of `sampleN.TOST()`: Same name of replicate designs like in the scaled functions.
  * In `power.scABEL()` name of scaled component `"p(BE-ABEL)"` instead of `"p(BE-wABEL)"`.
  * Updated `test_ABEL.R` in `inst/tests` to support subject simulations.
  * `power.scABEL()`: More informative warning about heteroscedasticity in the partial replicate design (use of `power.scABEL.sdsims()` suggested).
  * Add/subtract `.Machine$double.eps` if `rho` is -1 or +1 given in `power_type1_2TOST.R` (similar to `sampleN_2TOST_sim.R`). Removed warning in the latter function. Less confusing for users and the example in the man-page looks nicer.
  * Man pages reworked.

# PowerTOST 1.4-7

On CRAN 2018-04-12 (dedicated to the birthday of Alfie Schütz. :-)

## Bug fixes

  * Correction in functions using expected power.

## Major changes

  * Function `type1error.2TOST()` no longer available since it suffers from insufficient precision to obtain the type 1 error (TIE) via simulations. Due to the intersection-union principle the TIE is always upper bounded to alpha by theory.
  * Function `power.2TOST()` based on simulations to obtain the power of 2 TOSTs (statistical flaw in 4-dimensional *t*-distribution approach). `sampleN.2TOST()` and man page adapted according to this change.
  * Argument `regulator="FDA"` implemented in `scABEL.ad()`.
  * New function `power.RSABE2L.sds()` which implements the 'exact' based method for RSABE (ncTOST) of the two Lászlós. Documentation included.

# PowerTOST 1.4-6

On CRAN 2017-08-17.

## Bug fixes

  * Functions relying on simulations if `nsims > 1e7`.

## Major changes

  * Functions `CVwRfromU()` / `U2CVwR()` to calculate CV<sub>wR</sub> from the upper expanded  limit of an ABEL study according to the EMA’s or Health Canada’s rules.
  * Power and sample size for TOST: Argument `alpha` restricted to scalar. Internal functions now allow length = 2 (different alphas for the two null hypotheses).
  * Deprecated argument `dfCV` in expected power functions completely removed.

## Minor changes

  * Man-pages for non-inferiority again updated.
  * Data lazy loading (allows access by the name).

# PowerTOST 1.4-5

On CRAN 2017-05-19.

## Bug fixes

  * `power.scABEL.sdsims()` fixed which gave `power=NA` if `alpha=0`.
  * Correction in `sampleN.scABEL.ad()` if pre-specified alpha gives TIE <0.05.
  * Small bug in the S3 method plot for class `'pwrA'` corrected (text in the 4<sup>th</sup> screen was truncated on top if saved as a PDF). Replaced `text()` by `TeachingDemos::shadowtext()` to enhance legibility (interfered with underlying grid before).

## Major changes

  * Two new test scripts in `/test` subdirectory added: `test_ABEL.R` and `test_RSABE.R` which recalculate the sample size tables given in Tóthfalusi & Endrényi “*[Sample Sizes for Designing Bioequivalence Studies for Highly Variable Drugs](https://ejournals.library.ualberta.ca/index.php/JPPS/article/download/11612/9489)*”. Contributed by Helmut.
  * `power.scABEL.sdsims()` now has an argument `design_dta` to specify the design via a data.frame. May be useful for considering missing data.  
    **Attention!** This feature is experimental because the data.frame is currently not checked.
  * Updated functions `scABEL.ad()` and `sampleN.scABEL.ad()` to allow subject data simulations via `power.scABEL.sdsims()` if `regulator = "EMA"`. Removed `regulator = "ANVISA"`. Changed the order of sequences to be consistent with the other functions of PowerTOST.
  * Removed `regulator = "ANVISA"` from `pa.scABE()`.

## Minor changes

  * Man-pages for non-inferiority functions updated to include a description of the underlying hypotheses.
  * Remaining man-pages updated where term `'Null (true) ratio'` was mentioned instead of `'True ratio'`.
  * Internal change of coding in design helpers (repeated creation of data.frame with design characteristics avoided).
  * Typo corrected in `Expected_Power_for_TOST.pdf`.
  * Renamed the variable `'adj. alpha'` in `sampleN.scABEL.ad()` to `'alpha.adj'`for consistency with `scABEL.ad()`.

# PowerTOST 1.4-4

On CRAN 2017-03-15.

## Bug fixes

  * Use of removed functions `power.scABEL2()` and `sampleN.scABEL2()` in `pa.scABE()` corrected.

## Major changes

  * Added function `power.scABEL.sdsims()` to calculate power for the BE decision  via scaled (widened) BE acceptance limits (EMA recommended) based on subject data simulations.
  * Removed default `CV = 0.3` from `scABEL.ad()` and `sampleN.scABEL.ad()`. Stops execution if CV not specified.

## Minor changes

  * Administrative change for expected power: `cubature::adaptIntegrate` replaced by `cubature::hcubature` to reflect change of function name within package cubature.
  * Low CV (cosmetic) correction in `power.scABEL()` introduced.
  * DOI for references added in many man-pages. THX to Helmut.
  * Ben’s description of expected power added in `/doc` subdirectory.
  * `BE_power_sample_size_excerpt.pdf` updated to reflect the changes in computation of OwensQ. 

# PowerTOST 1.4-3

On CRAN 2016-11-01.

## Bug fixes

  * Little bug in `power.TOST()` removed which caused power <0 for large degrees of freedom.
  * Little inconsistency resolved (by Helmut) in output of power analysis  functions w.r.t. ratio or percent. Thanks to user myttiri of the BEBA forum for pointing out this inconsistency.

## Major changes

  * Functions for expected power of the TOST procedure and sample size based on expected power reworked to deal also with uncertainty of theta0 or dealing with both uncertainties of CV and theta0. Contributed mainly by Benjamin Lang.
  * Deprecated functions `power.scABEL2()` and `sampleN.scABEL2()` removed.
  * `regulator = "ANVISA"` no longer allowed in the scaled ABEL functions.
  * `OwensQ()` simpler/faster implemented. Now based solely on numerical integration in combination with non-central *t*.
  * Owen’s T-function, used in `OwensQOwen()` now based on algorithm AS76 and remarks to that algorithm to avoid numeric errors of the implementation via `integrate()`.

## Minor changes

  * Data.frame `ct9.6.6` in `data("data2x2x3")` added which was missing since a long time ago.

# PowerTOST 1.4-2

On CRAN 2016-07-14.

## Bug fixes

  * Minor bug in `sampleN.2TOST()` fixed.

## Major changes

  * Deprecated argument `point` in functions `CI2CV()`/`CVfromCI()` removed.

## Minor changes

  * `Implementation_scaledABE_sims.pdf` in `/doc` subdirectory updated to reflect changes in code of the scaled ABE functions.
  * The S3 method plot for class `'pwrA'` now has an argument `ratiolabel` for labeling the axis concerning theta0. Wish of Benjamin Lang.
  * Various enhancements in man-pages.
  * Power and sample size for FDA RSABE now take into account a CVcap if defined as finite.
  * Misleading term `'Null (true) ratio'` in output of sample size functions changed to `'True ratio'`.

# PowerTOST 1.4-1

Published on GitHub 2016-06-14.

## Major changes

  * Objects of class `'regSet'` have an additional component `'est_method'` which controls the simulations via key statistics of the evaluation using the EMA’s ANOVA or the FDA’s recommended ISC in `xyz.scABEL()` functions.
  * `power.scABEL()`/`power.scABEL2()` as well as `sampleN.scABEL()`/`sampleN.scABEL2()` are now unified. The regulator component `'est_method'` is used for switching between simulations based on the EMA’s ANOVA evaluation or ISC evaluation, respectively.  
    `power.scABEL2()`/`sampleN.scABEL2()` are therefore deprecated and will be removed in future versions. A corresponding warning is thrown.

## Minor changes

  * URL of PowerTOST on GitHub added, URL for bug reports added to `DESCRIPTION`.

# PowerTOST 1.3-7

Published on GitHub 2016-06-10.

## Bug fixes

  * In functions `scABEL.ad()` and `sampleN.scABEL.ad()` if former argument `regulator="ANVISA"` settings are used.

# PowerTOST 1.3-6

On CRAN 2016-06-06.  
Released to beta-testers 2016-05-04.

## Major changes

  * Functions for power and sample size of BE decision of highly variable drugs or drug products via ABEL (average BE with expanding limits) updated to incorporate the regulatory settings of Health Canada (new functions `power.scABEL2()` and `sampleN.scABEL2()` based on simulations of intra-subject contrasts evaluation).
  * The same functions now accept as regulator argument an object of class `'regSet'` and allow via this way User definitions of regulatory settings.
  * New function `reg_const` made visible to define objects of `'regSet'`. Class `'regSet'` has an S3 print method.

## Minor changes

  * Default value of argument `theta0` in functions `power.scABEL()`, `power.scABEL2()` and `power.RSABE()` changed to 0.90 to be in agreement with the setting in the corresponding sample size functions.
  * Function `scABEL()` (calculation of widened acceptance limits) no longer accepts `regulator="USER"`. This case may be handled via an object with class `'regSet'`, as defined by help of function `reg_const()`.
  * Functions `CVfromCI()`/`CI2CV` now use `pe` instead of `point` as argument due to more consistency with their dual `CI.BE()`. For backward compatibility `point` may be used also but then a warning is thrown. Argument point will be removed in future versions.
  * Functions for expected power for TOST and non-inferiority updated to avoid numeric dificulties with `integrate()` if `method="exact"`.
  * Documentation (typos in man-pages) rigorously enhanced. Thanks to Helmut and Ben.

# PowerTOST 1.3-5

On CRAN 2016-04-12.

## Bug fixes

  * plot method of `'pwrA'` objects reworked by Helmut Schütz to avoid some ugly overlays.

## Major changes

  * Functions for expected power for TOST and non-inferiority updated (including corresponding sample size functions). Functions now include an exact method which calculates the expected value of the power with respect to the (prior) distribution of σ<sup>2</sup> (inverse Γ distribution). Formerly only an approximation according to Julious/Owen was implemented. Contributed by Benjamin Lang.

## Minor changes

  * The functions for iteratively adjusting alpha for the EMA recommended ABEL method now have a new argument `tol` for convergence tolerance of the internally used `uniroot()`. It was formerly hardcoded as `tol=1e-5`. Contributed by Helmut Schütz.

# PowerTOST 1.3-4

On CRAN 2016-03-09.

## Bug fixes

  * In all sample size estimation functions for scaled ABE which resulted in `NA` for sample size for very low variability. Thanks to Shuanghe for detecting this bug.

## Minor changes

  * All functions for scaled ABE now have an argument `imax=100` for the maximum number of steps for sample size search. Wish of Helmut Schütz.

# PowerTOST 1.3-3

On CRAN 2016-01-15.  
Released 2016-01-04 to alpha testers.

## Major changes

  * New functions `scABEL.ad()` and `sampleN.scABEL.ad()` to iteratively adjust α in order to control the consumer risk and adapt the sample size to compensate for a potential loss in power with the EMA method of scaled ABE (ABEL). Contributed by Helmut Schütz.

## Minor changes

  * Default `theta0` changed to 0.9 in sample size estimation for scaled ABE (as recommended by the two Lászlós).

# PowerTOST 1.3-2

On CRAN 2015-12-02.

## Major changes

  * New functions for calculation of power, type 1 errror and sample size for 2 TOST instances. Contributed by Benjamin Lang.

# PowerTOST 1.3-1

On CRAN 2015-09-30.  
Released 2015-09-23 to alpha testers.

## Major changes

  * New functions `power.HVNTID()` and `sampleN.HVNTID()` introduced to calculate power and sample size for the BE decision via the FDA procedure for highly variable NTIDs (see the FDA’s dabigatran / rivaroxaban guidances).

## Minor changes

  * `power.HVNTID()` and `power.NTIDFDA()` now return the power (`p(BE)`) and the components for the scaled ABE criterion, the conventional ABE test and the test for the ratio s<sub>wT</sub>/s<sub>wR</sub> <= 2.5 if the argument `details=TRUE` (wish of Helmut Schütz).
  * `power.RSABE()` now returns the power (`p(BE)`) and the components for the scaled ABE criterion, for the point estimate criterion and for the conventional ABE test alone if the argument `details=TRUE`. Analogous changes were made in `power.scABEL()`. See man-pages.

# PowerTOST 1.2-9

On CRAN 2015-08-26.

## Bug fixes

  * `power.TOST()` with `method="exact"` reworked to give correct values in case `alpha>0.5`. `OwensQOwen()` adapted to deal with upper integration limit `R==Inf`.

## Major changes

  * New power calculation method introduced: Direct integration of the bivariate non-central *t*-distribution (`pmvt()` of package `mvtnorm`), also an exact calculation method but with somewhat lower precision and longer run-time. Contributed by Benjamin Lang.
  * New function `power.TOST.sim()` to obtain the power via simulations of the TOST. Only intended for checking purposes.

# PowerTOST 1.2-8

On CRAN 2015-07-10.

## Minor changes

  * Maintenance release with `NAMESPACE` adapted to import the necessary functions from base R installation (CRAN request).
  * `power.NTIDFDA()` and `sampleN.NTIDFDA()` now have an design argument to choose between `"2x2x4"` (full replicate 4-period) and `"2x2x3"` (full replicate 3-period design).
  * Some minor documentation improvements.

# PowerTOST 1.2-7

On CRAN 2015-06-03.

## Major changes

  * Functions `power.scABEL()` and `sampleN.scABEL()` now allow the power and sample size calculations for scaled ABEL according to the (inofficial) regulatory settings of Brasilian ANVISA (argument `regulator="ANVISA"`).
  * Function pa.scABE() also allows `regulator="ANVISA"`.
  * New helper function `scABEL` for calculation of the (widened) ABE acceptance limits.
  * Deprecated function `power2.TOST()` removed.

# PowerTOST 1.2-6

On CRAN 2015-01-23.

## Major changes

  * Function `power.TOST()` now handles balanced as well as unbalanced studies. Function `power2.TOST()` is therefore deprecated and will be removed in later versions.

## Minor changes

  * Internal change in interface to the (hidden) 'raw' power functions (sem used in function calls instead of se, n, bk).
  * `BE_power_sample_size_excerpt.pdf` in subdirectory `/doc` changed to reflect the internal changes.

# PowerTOST 1.2-5

On CRAN 2015-01-07.

## Bug fixes

  * Bug in `OwensQOwen()` removed (debug code left over), THX to Helmut for discovering this.

## Minor changes

  * Function `CVfromCI()` and alias `CI2CV()` now handle unbalanced studies. Contributed by Benjamin Lang.
  * Typos in `BE_power_sample_size_excerpt.pdf` (discovered by Ben) corrected.

# PowerTOST 1.2-4

On CRAN 2014-12-19.

## Bug fixes

  * Bug in `OwensT(h, a)` for the case of `a=-Inf` removed. Thanks to Benjamin Lang for discovering that nasty bug.

## Major changes

  * Function `pvalue.TOST()` introduced to calculate the two p-values of the TOST procedure. Contributed by Benjamin Lang.

## Minor changes

  * Functions `power.noninf()`, `exppower.TOST()`, `exppower.noninf()` adapted to deal with unbalanced studies.

# PowerTOST 1.2-3

On CRAN 2014-11-13.

## Bug fixes

  * Minimum sample size in `sampleN.dp()` introduced to avoid errors with small CV (<= 0.1) for the crossover design.

## Major changes

  * Functions `power.dp()` / `sampleN.dp()` cover "incomplete block designs".
  * Liu’s 2x2x2 design with 2 repeated measurements in each period added (see `?known.designs` under Notes).

## Minor changes

  * Partitition of n(total) to (sequence) groups reworked in most power functions if design is unbalanced.

# PowerTOST 1.2-2

On CRAN 2014-10-06.

## Minor changes

  * Functions `pa.XYZ()` adapted so that they also work under R < 3.1.1 (request of Uwe Ligges of CRAN) although with `minpower >= 0.5`.

# PowerTOST 1.2-1

On CRAN 2014-09-30.

## Bug fixes

  * Some minor bugs removed and spelling in man-pages corrected.

## Minor changes

  * Print method for class `pwrA` now calls `plot()`.
  * Power analysis in case of n=12 for the plan now drops n below 12 subjects.

# PowerTOST 1.2-0

Build 2014-09-19, released to alpha/beta testers only.

## Major changes

  * Functions for power analysis of a sample size plan for ABE (`pa.ABE()`), scaled ABE (`pa.sABE()`) and scaled ABE for NTIDs (`pa.NTIDFDA()`) analyzing power for deviations from assumptions for the sample size estimation.
  * Experimental functions for power calculations / sample size estimation for dose proportionality studies using the power model.

# PowerTOST 1.1-13

On CRAN 2014-08-12.

## Bug fixes

  * Design constant `bk(ni)` for `design="paired"` corrected. THX to Helmut Schütz who detected this bug.

# PowerTOST 1.1-12

On CRAN 2014-07-02.

## Minor changes

  * Improvements in documentation.
  * Function `power.TOST()` now throws a warning if used with imbalanced designs (n not an even multiple of the number of sequences).
  * Internal change of start value of sample size search in `sampleN.TOST()` to avoid failed searches if variability is high and theta0 is close to 1.

# PowerTOST 1.1-11

On CRAN 2014-04-30.

## Major changes

  * Utility function added which calculates 1&ndash;2α confidence interval(s) given point est., CV and n using log-tansformed evaluation.
  * Utility function added which calculates 1&ndash;2α Fieller confidence interval(s) given point est., CV (, CVb) and n for the ratio of untransformed means.

# PowerTOST 1.1-10

On CRAN 2014-01-31.

## Bug fixes

  * Bug removed where noncentrality parameters in the power calculations became `NaN` (not a number).

## Major changes

  * The package is now pre-compiled to 'byte code' via the compiler package for speed reasons. F.i. `OwensQOwen()` gains a sixfold speed boost.
  * Some bulk code changes to make the power calculations for extreme cases more bullet proof without computation time burden.
  * `OwensQ()` now tries to return a value via nct-approximation if its value is due to numeric difficulties falsely equal zero. This approximation is up to 6 decimals correct as far as tested. `OwensQ()` issues a warning if the nct-approximation is used.

## Minor changes

  * Improvements in documentation (man-pages and \*.pdf) to reflect the code changes.

# PowerTOST 1.1-9

Build 2014-01-03, not released to the public.

## Major changes

  * `OwensQ()` now uses `OwensQOwen()` in case of high delta and/or high upper integration limit. Thus extreme cases can be handled properly where the former implementation via `integrate()` was prone to fail. Thanks to Jiři Hofmann and Helmut Schütz for pointing me to such extreme cases.

# PowerTOST 1.1-8

On CRAN 2013-12-27.

## Bug fixes

  * Date typos in the history (`NEWS`) corrected from PowerTOST 1.1-00 on. THX to Julien Grassot.

## Minor changes

  * `CVfromCI()` now accepts either both CLs or one CL and the point estimate. Contributed by Helmut Schütz.
  * `power.RatioF()` and `sampleN.RatioF()` now have an argument `setseed=TRUE` which avoids the dependence of the power from the state of the random number generator (due to the calculation method of `pmvt()` of package `mvtnorm`). Thanks to Benjamin Lang.

# PowerTOST 1.1-7

On CRAN 2013-09-02.

## Major changes

  * `design="2x2x3"` (TRT|RTR) implemented in `power.scABEL()`, `sampleN.scABEL()` and in `power.RSABE()`, `sampleN.RSABE()`.

## Minor changes

  * `power.scABEL()` now throws a warning if CV<sub>wT</sub> ≠ CV<sub>wR</sub> in the design `"2x3x3"` (partial replicate).
  * Simulation details for the full replicate design slightly changed to obtain better numeric agreement of power to subject data sims.

# PowerTOST 1.1-6

On CRAN 2013-06-21.

## Bug fixes

  * Fat bug corrected in `sampleN.NTIDFDA()` which lead to false sample size for cases where the test of equal variabilities of Test vs. Reference comes into effect.

# PowerTOST 1.1-5

On CRAN 2013-06-17.

## Bug fixes

  * Flaw in sample size search corrected which does not indicate a failed search (`n=NA`) if started with a too high n. Thanks to Helmut Schütz.

## Major changes

  * Functions `sampleN.scABEL()`, `sampleN.RSABE()` and `sampleN.NTIDFDA()` return a data.frame with the input and the sample size result. The `"Sample size"` column contains the total sample size. The `"nlast"` column contains the last n value handled. Might be useful for restarting.

## Minor changes

  * Start values for sample size search in functions `sampleN.scABEL()`, `sampleN.RSABE()` reworked. Failed sample size searches are now more seldom observed.

# PowerTOST 1.1-4

Build 2013-05-15, relased to alpha testers only.

## Major changes

  * Functions `power.NTIDFDA()` and `sampleN.NTIDFDA()` introduced to calculate power and sample size for the BE decision via the method of the FDA for narrow therapeutic index drugs (NTIDs, for details see the FDA’s warfarin guidance). Power and sample size are based on simulations.

# PowerTOST 1.1-3

On CRAN 2013-05-03.

## Bug fixes

  * Flaw in implementing simulations of `"2x3x3"` design with different intra-subject variabilities in functions `power.RSABE()` and `sampleN.RSABE()` as well as in functions `power.scABEL()` and `sampleN.scABEL()` corrected.

## Minor changes

  * Methods and implementation details of the simulations for scaled ABE documented in a PDF file in the `/doc` subdirectory of the package.
  * Simulation method for the EMA scaled ABEL changed to conform better with power values from simulations via subject data.
  * Warning section in the help file of `power.scABEL()` introduced to reflect the fact that simulations via subject data and via the methods implemented in `power.scABEL()` gave empirical power values that are only approximately in agreement.

# PowerTOST 1.1-2

On CRAN 2013-02-28.

## Major changes

  * Functions `power.RSABE()` and `sampleN.RSABE(`) introduced to calculate power and sample size for the BE decision via linearized scaled ABE criterion as favored method of the FDA. Power and sample size are based on simulations.  
    Default of `nsims` changed to 1E5 (suggested by Helmut Schütz).

## Minor changes

  * Argument `setseed` introduced in the scaled ABE functions to avoid different outcomes depending on the state of the (pseudo) random number generator. If `setseed=TRUE` a `set.seed(123456)` is issued prior to each call of the simulation functions.
  * Documentation improved.

# PowerTOST 1.1-0

On CRAN 2013-02-08.

## Major changes

  * Functions `power.scABEL()` and `sampleN.scABEL()` introduced to calculate power and sample size for the BE decision via scaled (widened) BE acceptance limits based on simulations. Thanks to Helmut Schütz on the [BEBA Forum](https://forum.bebac.at/mix_entry.php?id=9997) for power testing these functions.

# PowerTOST 1.1-01

Not released, integrated in 1.1-00

## Major changes

  * Upper one-sided CL of the CV and therefore argument `alpha2` removed from `expsampleN.noninf()` and `expsampleN.TOST()` because it lead to some confusion in users thinking this had to do with the algorithm via expected power.

## Minor changes

  * `CI2CV()` as alias to `CVfromCI()` introduced because myself always typed this name if aimed to calculate the CV from a given CI.
  * Function `CVCL()` now returns a 2-element vector also if an one sided interval is requested.

# PowerTOST 1.0-0

On CRAN 2012-10-26.

## Bug fixes

  * Bugfix in `power.noninf()` to get the correct power if `theta0` is below `margin` (if `margin <1`) or `theta0` is above `margin` (if `margin >1`). `power.noninf()` calculated up to now the power of an inferiority test. Thanks to Helmut Schütz.

## Major changes

  * New helper functions `CV2mse()` and `mse2CV()`.
  * New function `CVCL()` to calculate a confidence interval of a CV.

# PowerTOST 0.9-11

On CRAN 2012-08-07.

## Bug fixes

  * In `power.TOST()`, `power2.TOST()` and `power.noninf()` to use the correct degrees of freedom depending on argument `robust`. `robust=FALSE` wrongly used the robust dfs. Thanks to Ben.

# PowerTOST 0.9-10

On CRAN 2012-07-20.

## Bug fixes

  * In internal function `.Q.integrand()`.
  * In `nmin` - has to be a multiple of steps to assure balance in sequence groups.  
    Again thanks to Helmut Schütz for detecting both.

# PowerTOST 0.9-9

On CRAN 2012-07-18.

## Bug fixes

  * Workaround introduced to handle numeric problems in `integrate()` if `CV<5.3E-6`. Thanks to Helmut Schütz.

## Minor changes

  * Minimum sample size adapted to design used (f.i. `n=2` if `design="paired"`).

# PowerTOST 0.9-8

On CRAN 2012-04-05 (Easter egg).

## Major changes

  * Functions added for 'expected' power and sample size calculations based on it for the non-inferiority test for sake of completeness.
  * Functions `CVpooled()`, `exppower.TOST()` and `expsampleN.TOST()` now also implemented for `logscale=FALSE`, *i.e.*, contain this argument in their calls.
  * Function `OwensQOwen()` made public. This is an implementation of the algorithm given by Owen in the original paper (Biometrica 1965) via repeated integration by parts.  
    This function is only for comparative purposes.
  * Function `OwensT()` made public. It is needed internally in `OwensQOwen()` but may be useful for other purposes.

# PowerTOST 0.9-6/7

On CRAN 2012-03-26.

## Bug fixes

  * PowerTOST 0.9-7 is a small bug fix of PowerTOST required by CRAN. It contains the old release number 0.9-6 in the package man-page.

## Major changes

  * Functions added for power and sample size calculations based on non-inferiority t-test. This is not a TOST procedure but eventually useful if the question of 'non-superiority' within a BE study must be evaluated.  
    Hint: Evaluation of Fluctuation in the EMA MR NfG (1999) between modified release formulation and immediate release product.

# PowerTOST 0.9-4

On CRAN 2012-03-05.

## Bug fixes

  * Little bug in `sampleN.TOST()` removed which causes extra doubled output of n and power if `n=4`. Thanks to Ben on the [BEBA Forum](https://forum.bebac.at/mix_entry.php?id=8206).

# PowerTOST 0.9-3

On CRAN 2012-02-13.

## Bug fixes

  * Bug in `power.TOST()` removed which prevented calculation of power according to `method = "exact"`.

## Minor changes

  * Sample size tables for replicate design `"2x2x3"` in `/data` sub-directory added.
  * Sample size tables for replicate design `"2x4x4"` in `/data` sub-directory added.
  * Scripts in the `/test` sub-directory made public.

# PowerTOST 0.9-2

On CRAN 2011-12-24.

## Major changes

  * Function `power2.TOST()` added to allow power calculations for studies with unbalanced (sequence) groups.
  * Argument `exact` replaced by `method` in `power.TOST()`, `sampleN.TOST()`. See help for these functions.
  * **Attention!** The sample size for the parallel group design is now the *total* sample size (to be consistent across all functions).

## Minor changes

  * Sample size tables added for the `"2x2"` crossover and for the `"parallel"` group design to alleviate validation/qualification of the package.  
    See `data(package="PowerTOST")`.
  * Scripts added in the `/test` sub-directory that create the sample size tables from the data section.
  * Updated `BE_power_sample_size_excerpt.pdf` in the `/doc` sub-directory.

# PowerTOST 0.9-0

On CRAN 2011-12-15.

## Major changes

  * Paired means design introduced.
  * `'robust'` argument added to nearly all functions.  
    With `robust=TRUE` the degrees of freedom for the so-called robust evaluation (`df2` in `known.designs()` output) will be used.  
    This may be helpful if planning is done for higher order designs evaluated via mixed model or via intra-subject contrasts (aka Senn’s basic estimator).
  * Due to the necessary `NAMESPACE` from R1.4.0 on the internal functions (names starting with `.`) are no longer exported.

# PowerTOST 0.8-7

On CRAN 2011-10-20.

## Bug fixes

  * Problem with slash in \name field of manual resolved (requested by Brian Ripley).

# PowerTOST 0.8-6

On CRAN 2011-05-18.

## Bug fixes

  * Bug removed which gave incorrect exact power values in case of `alpha>=0.5` (very unusual setting). Thanks again to Craig Zupke.  
    Cross checked results of power at equivalence margins against SAS `Proc Power`.

# PowerTOST 0.8-5

On CRAN 2011-05-16.

## Bug fixes

  * Code in Owen’s Q adapted to account for large delta or large `b` leading to integrand function almost zero over the whole range which then gave an error in `integrate()`. Thanks to Craig Zupke.


# PowerTOST 0.8-4

On CRAN 2011-03-11.

## Minor changes

  * Number of maximal steps in sample size search in the sample size functions, formerly hard coded as 50, made accessible to users via argument `imax`. Needs to be adapted only in rare extreme cases.
  * Start value for sample size search improved around `theta0=1` (logscale) or `theta0=0` (untransformed).

# PowerTOST 0.8-3

On CRAN 2011-01-18.

## Bug fixes

  * In `known.designs()`.

# PowerTOST 0.8-2

On CRAN 2011-01-10.

## Bug fixes

  * Error in df for `"3x3"` and `"4x4"` crossover designs removed.

## Major changes

  * Function for pooling CVs of different studies made public. See `?CVpooled`.

# PowerTOST 0.8-1

On CRAN 2010-11-25.

## Major changes

  * Helper function `CVfromCI()` added to estimate the CV from a confidence interval. Useful if no CV but the CI was given in literature.

# PowerTOST 0.7-3

On CRAN 2010-11-25.

## Major changes

  * Input argument `diff` removed from `sampleN.TOST()`, `expsampleN.TOST()`, `power.TOST()`, and `exppower.TOST()`.

## Minor changes

  * Documentation improved.

# PowerTOST 0.7-2

On CRAN 2010-08-27.

  * Little bug causing warnings in case of the `"2x2"` alias `"2x2x2"` design corrected.

# PowerTOST 0.7-1

On CRAN 2010-08-12.

## Major changes

  * Functions added for the power and sample size for the ratio of two means with normally distributed data on the original scale (based on Fieller’s confidence (‘fiducial’) interval).  
    AFAIK, until now only implemented in the commercial nQuery.

## Minor changes

  * The argument `diff` (Null ratio / Null diff.) is now named `theta0` since it was annoying for users to call it `diff` in case of ratios (`logscale=TRUE`). The parameter `diff` is still supported but will be removed in the next release. Therefore a warning is thrown if `diff` is used.

# PowerTOST 0.6-2

On CRAN 2010-07-21.

## Major changes

  * Previously internal hidden functions `.CV2se()` and `.se2CV()` made public.

## Minor changes

  * Some internal code consolidation
  * Minor enhancements in help pages, more examples.
  * Short documentation of used statistical apparatus `BE_power_sample_size_excerpt.pdf` for classical power / sample size in directory `/doc` added.

# PowerTOST 0.5-1

On CRAN 2010-05-07, first public release.
  
