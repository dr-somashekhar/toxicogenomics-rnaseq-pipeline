#!/usr/bin/env Rscript
# One-time setup: installs every R/Bioconductor package the pipeline needs, plus the Python
# TensorFlow/Keras backend. Safe to re-run; already-installed packages are skipped.
#   Rscript install_dependencies.R

if (getRversion() < "4.2.0") stop("R >= 4.2 is required (found ", getRversion(), ").")

options(repos = c(CRAN = "https://cloud.r-project.org"))

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

cran_pkgs <- c("optparse", "futile.logger", "WGCNA", "keras", "tensorflow", "reticulate",
               "xgboost", "caret", "SHAPforxgboost", "ParBayesianOptimization", "tidyverse",
               "doParallel", "testthat", "yaml", "matrixStats")
bioc_pkgs <- c("DESeq2", "sva", "GSVA", "GSEABase", "apeglm", "BiocParallel")

missing_pkgs <- function(pkgs) pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]

todo_cran <- missing_pkgs(cran_pkgs)
if (length(todo_cran)) install.packages(todo_cran)
todo_bioc <- missing_pkgs(bioc_pkgs)
if (length(todo_bioc)) BiocManager::install(todo_bioc, update = FALSE, ask = FALSE)

# Python backend for the autoencoder (creates an isolated virtualenv managed by reticulate).
if (!reticulate::py_module_available("tensorflow")) {
  message("Installing TensorFlow/Keras Python backend ...")
  keras::install_keras()
}

still_missing <- missing_pkgs(c(cran_pkgs, bioc_pkgs))
if (length(still_missing)) {
  stop("Could not install: ", paste(still_missing, collapse = ", "))
}
message("All dependencies installed. TensorFlow available: ",
        reticulate::py_module_available("tensorflow"))
