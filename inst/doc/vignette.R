## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#"
)


## -----------------------------------------------------------------------------
package <- "PowerTOST"
inst    <- package %in% installed.packages()
if (length(package[!inst]) > 0) install.packages(package[!inst])


## ---- sessioninfo-------------------------------------------------------------
options(width = 80)
devtools::session_info()

