# R/zzz.R

.onLoad <- function(libname, pkgname) {
  if (getRversion() >= "2.15.1") {
    utils::globalVariables(c("_MLFApackage_MLFA", "res_MLFA2", "variable", "value", "Iteration"))
  }
}


