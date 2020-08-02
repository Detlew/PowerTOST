# precompute the vignette with long run-time according to
# https://ropensci.org/technotes/2019/12/08/precompute-vignettes/

# Execute the code from the vignette rmd renamed to .Rmd.orig
knitr::knit("vignettes/RSABE.Rmd.orig", output = "vignettes/RSABE.Rmd")
