#!/usr/bin/env Rscript
# ==================================================================================================
#  AI-INTEGRATED TOXICOGENOMICS & MULTI-OMICS PREDICTIVE MODELING ENGINE
#  Version: 4.1.0  |  DILI (Drug-Induced Liver Injury) Transcriptomic Risk Modeling
# ==================================================================================================
#
#  ARCHITECTURE
#  ------------
#  01_config_and_setup   -- CLI args, config object, logging, parallel backend, reproducibility
#  02_data_simulation     -- parameterized synthetic multi-omics cohort generator
#  03_batch_correction    -- ComBat-seq + PCA-based QC diagnostics pre/post correction
#  04_differential_expr   -- DESeq2 LRT model, independent filtering, apeglm shrinkage
#  05_network_biology      -- WGCNA soft-thresholding, module detection, module-trait correlation,
#                              hub gene identification via intramodular connectivity
#  06_pathway_enrichment   -- GSVA single-sample pathway enrichment as an engineered feature block
#  07_representation_learn -- Keras autoencoder (real fit/predict, checkpointed, early-stopped)
#  08_predictive_model     -- Nested CV + Bayesian hyperparameter search for XGBoost
#  09_model_interpretation -- SHAP value decomposition of the trained model
#  10_diagnostics          -- Bootstrap CI for RMSE/R2, residual + calibration diagnostics
#  11_persistence          -- Model/artifact serialization + session provenance
#  12_unit_tests           -- testthat coverage for the simulation & preprocessing layer
#  13_main                 -- Orchestration entry point (guarded, so file is sourceable OR runnable)
#
#  DESIGN NOTES
#  ------------
#  - Every pipeline stage is a pure-ish function that takes explicit inputs and returns explicit
#    outputs (no reliance on floating global state), so stages can be unit tested in isolation.
#  - The autoencoder is ACTUALLY trained here (no commented-out placeholder branch). If keras/
#    tensorflow are unavailable in the runtime environment, `build_and_train_autoencoder()`
#    fails loudly rather than silently substituting random noise for a "latent space."
#  - Model selection never touches the held-out test set: hyperparameters AND the number of boosting
#    rounds are chosen by inner CV on the training split only.
#  - Logging replaces ad hoc `cat()` calls so the pipeline produces an audit trail suitable for
#    a lab notebook or CI log.
# ==================================================================================================

# --------------------------------------------------------------------------------------------------
# 01. CONFIG, LOGGING, CLI, AND ENVIRONMENT SETUP
# --------------------------------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(optparse)     # CLI argument parsing
  library(futile.logger) # Structured logging
  library(DESeq2)
  library(WGCNA)
  library(keras)
  library(tensorflow)
  library(xgboost)
  library(caret)
  library(sva)
  library(GSVA)
  library(GSEABase)
  library(SHAPforxgboost)
  library(ParBayesianOptimization)
  library(tidyverse)
  library(doParallel)
  library(testthat)
  library(yaml)
  library(BiocParallel)
})

# Be explicit about masked generics (AnnotationDbi / WGCNA mask dplyr::select and stats::cor).
select <- dplyr::select
cor    <- stats::cor

# ---- 1.1 CLI ARGUMENTS -----------------------------------------------------------------------
# Allows `Rscript ai_toxicogenomics_master_pipeline.R --n_samples 500 --seed 2026 --run_tests TRUE`
option_list <- list(
  make_option("--n_samples", type = "integer", default = 500,
              help = "Number of simulated hepatocyte samples [default %default]"),
  make_option("--n_genes", type = "integer", default = 20000,
              help = "Number of simulated genes [default %default]"),
  make_option("--seed", type = "integer", default = 2026,
              help = "Global random seed for full reproducibility [default %default]"),
  make_option("--latent_dim", type = "integer", default = 32,
              help = "Autoencoder bottleneck dimensionality [default %default]"),
  make_option("--bo_iterations", type = "integer", default = 15,
              help = "Bayesian optimization iterations for XGBoost tuning [default %default]"),
  make_option("--output_dir", type = "character", default = "pipeline_artifacts",
              help = "Directory for models, logs, and figures [default %default]"),
  make_option("--run_tests", type = "logical", default = TRUE,
              help = "Run the testthat unit-test suite before the pipeline executes [default %default]"),
  make_option("--config_file", type = "character", default = NA,
              help = "Optional YAML file overriding any of the above defaults")
)

parse_cli_args <- function(args = commandArgs(trailingOnly = TRUE)) {
  parser <- OptionParser(option_list = option_list)
  opts <- parse_args(parser, args = args)

  if (!is.na(opts$config_file) && file.exists(opts$config_file)) {
    yaml_overrides <- yaml::read_yaml(opts$config_file)
    opts <- modifyList(opts, yaml_overrides)
  } else if (!is.na(opts$config_file)) {
    stop("config_file not found: ", opts$config_file)
  }
  validate_opts(opts)
  opts
}

validate_opts <- function(opts) {
  stopifnot(
    "n_samples must be >= 30 (needed for CV, WGCNA and the autoencoder)" = opts$n_samples >= 30,
    "n_genes must be >= 2000 (vst() and WGCNA need >1000 genes after filtering)" = opts$n_genes >= 2000,
    "latent_dim must be >= 2" = opts$latent_dim >= 2,
    "bo_iterations must be >= 1" = opts$bo_iterations >= 1
  )
  invisible(opts)
}

# ---- 1.2 LOGGING -----------------------------------------------------------------------------
init_logger <- function(output_dir) {
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  log_path <- file.path(output_dir, sprintf("pipeline_run_%s.log", format(Sys.time(), "%Y%m%d_%H%M%S")))
  flog.appender(appender.tee(log_path))
  flog.threshold(INFO)
  flog.info("Logger initialized. Writing to %s", log_path)
  invisible(log_path)
}

# ---- 1.3 PARALLEL BACKEND ---------------------------------------------------------------------
init_parallel_backend <- function() {
  detected <- parallel::detectCores()
  if (is.na(detected)) detected <- 2
  num_cores <- max(1, detected - 1)
  registerDoParallel(cores = num_cores)
  allowWGCNAThreads(nThreads = num_cores)
  # BiocParallel (used by DESeq2 and GSVA): fork on Linux/macOS, sockets on Windows.
  bp <- make_bpparam(num_cores)
  BiocParallel::register(bp, default = TRUE)
  flog.info("Parallel backend registered across %d cores.", num_cores)
  num_cores
}

make_bpparam <- function(num_cores) {
  if (num_cores <= 1) return(BiocParallel::SerialParam())
  if (.Platform$OS.type == "windows") BiocParallel::SnowParam(workers = num_cores)
  else BiocParallel::MulticoreParam(workers = num_cores)
}

# ---- 1.4 REPRODUCIBILITY SNAPSHOT ---------------------------------------------------------------
capture_provenance <- function(output_dir, config) {
  prov <- list(
    timestamp = as.character(Sys.time()),
    session_info = capture.output(sessionInfo()),
    config = config
  )
  saveRDS(prov, file.path(output_dir, "provenance.rds"))
  yaml::write_yaml(config, file.path(output_dir, "resolved_config.yaml"))
  flog.info("Provenance and resolved config snapshot written to %s", output_dir)
}



# --------------------------------------------------------------------------------------------------
# 02. SYNTHETIC MULTI-OMICS COHORT SIMULATION
# --------------------------------------------------------------------------------------------------
BATCH_LEVELS     <- c("Lab_A", "Lab_B", "Lab_C")
TREATMENT_LEVELS <- c("Vehicle", "Low_Dose", "High_Dose")

#' Simulate a dose-response hepatotoxicity RNA-seq cohort
#'
#' Counts are drawn from a negative binomial whose mean is the product of a per-gene baseline,
#' a dose-dependent module effect, and a per-batch (lab) gene-wise multiplicative factor, so the
#' batch-correction stage has a real technical effect to remove.
#'
#' @param n_samples number of simulated samples
#' @param n_genes number of simulated genes (module sizes scale with it)
#' @param seed RNG seed
#' @param batch_effect_sd SD of the log-normal gene-wise batch factor (0 disables batch effects)
#' @return list(counts = integer matrix [genes x samples], pheno = data.frame,
#'              ground_truth_modules = list of gene row indices)
simulate_multiomics_cohort <- function(n_samples = 500, n_genes = 20000, seed = 2026,
                                       batch_effect_sd = 0.5) {
  set.seed(seed)
  flog.info("Simulating cohort: n_samples=%d, n_genes=%d", n_samples, n_genes)

  pheno_data <- data.frame(
    Sample_ID = paste0("Hepatocyte_", seq_len(n_samples)),
    Batch = factor(sample(BATCH_LEVELS, n_samples, replace = TRUE), levels = BATCH_LEVELS),
    Treatment_Group = factor(
      sample(TREATMENT_LEVELS, n_samples, replace = TRUE, prob = c(0.4, 0.3, 0.3)),
      levels = TREATMENT_LEVELS
    )
  )
  rownames(pheno_data) <- pheno_data$Sample_ID

  pheno_data$Toxicity_Score <- with(pheno_data, {
    base <- ifelse(Treatment_Group == "Vehicle", rnorm(n_samples, 10, 5),
             ifelse(Treatment_Group == "Low_Dose", rnorm(n_samples, 40, 15),
                    rnorm(n_samples, 85, 10)))
    pmin(pmax(base, 0), 100)
  })

  # Per-gene baseline mean, expanded to a genes x samples matrix of NB means.
  base_means <- pmax(rnbinom(n_genes, mu = 800, size = 1.2), 1)
  mu <- matrix(base_means, nrow = n_genes, ncol = n_samples)

  # Ground-truth mechanistic gene modules (non-overlapping), sized as a fraction of n_genes.
  module_sizes <- pmax(1L, floor(n_genes * c(0.02, 0.02, 0.015)))
  mod_ox_stress <- sample(seq_len(n_genes), module_sizes[1])
  mod_apoptosis <- sample(setdiff(seq_len(n_genes), mod_ox_stress), module_sizes[2])
  mod_necrosis  <- sample(setdiff(seq_len(n_genes), c(mod_ox_stress, mod_apoptosis)), module_sizes[3])

  scale_module <- function(mu, genes, samples, lo, hi) {
    if (length(samples) == 0) return(mu)
    mu[genes, samples] <- mu[genes, samples] *
      matrix(runif(length(genes) * length(samples), lo, hi), nrow = length(genes))
    mu
  }
  high_idx <- which(pheno_data$Treatment_Group == "High_Dose")
  low_idx  <- which(pheno_data$Treatment_Group == "Low_Dose")
  mu <- scale_module(mu, mod_ox_stress, high_idx, 3, 10)
  mu <- scale_module(mu, mod_apoptosis, high_idx, 0.1, 0.5)
  mu <- scale_module(mu, mod_necrosis,  high_idx, 4, 12)
  mu <- scale_module(mu, mod_ox_stress, low_idx,  1.5, 3)

  # Technical batch effect: Lab_A is the reference; other labs get gene-wise log-normal factors.
  batch_factor <- matrix(1, nrow = n_genes, ncol = length(BATCH_LEVELS))
  if (batch_effect_sd > 0) {
    for (b in 2:length(BATCH_LEVELS)) batch_factor[, b] <- exp(rnorm(n_genes, 0, batch_effect_sd))
  }
  mu <- mu * batch_factor[, as.integer(pheno_data$Batch), drop = FALSE]

  raw_counts <- matrix(rnbinom(length(mu), mu = mu, size = 1.2), nrow = n_genes)
  stopifnot(max(raw_counts) < .Machine$integer.max)
  storage.mode(raw_counts) <- "integer"
  rownames(raw_counts) <- paste0("ENSG", sprintf("%06d", seq_len(n_genes)))
  colnames(raw_counts) <- pheno_data$Sample_ID

  list(
    counts = raw_counts,
    pheno = pheno_data,
    ground_truth_modules = list(
      oxidative_stress = mod_ox_stress,
      apoptosis = mod_apoptosis,
      necrosis = mod_necrosis
    )
  )
}


# --------------------------------------------------------------------------------------------------
# 03. BATCH CORRECTION + QC DIAGNOSTICS
# --------------------------------------------------------------------------------------------------
#' Largest R^2 between batch label and any of the first `n_pcs` principal components
batch_pc_r2 <- function(counts, batch, n_pcs = 5) {
  n_pcs <- min(n_pcs, ncol(counts) - 1)
  pca <- prcomp(t(log1p(counts)), center = TRUE, scale. = FALSE, rank. = n_pcs)
  r2 <- apply(pca$x, 2, function(pc) summary(lm(pc ~ batch))$r.squared)
  max(r2)
}

#' Run ComBat-seq batch correction and report a PCA-based sanity check
run_batch_correction <- function(raw_counts, pheno_data, output_dir) {
  flog.info("Running ComBat-seq batch correction across %d batches.", nlevels(pheno_data$Batch))

  pre_batch_r2 <- batch_pc_r2(raw_counts, pheno_data$Batch)

  adjusted_counts <- ComBat_seq(
    counts = as.matrix(raw_counts),
    batch = pheno_data$Batch,
    group = pheno_data$Treatment_Group
  )

  post_batch_r2 <- batch_pc_r2(adjusted_counts, pheno_data$Batch)

  flog.info("Max batch R2 over PC1-5 -- before: %.3f | after: %.3f", pre_batch_r2, post_batch_r2)

  qc_report <- data.frame(
    stage = c("pre_correction", "post_correction"),
    max_pc1to5_batch_r2 = c(pre_batch_r2, post_batch_r2)
  )
  write_csv(qc_report, file.path(output_dir, "batch_correction_qc.csv"))

  if (post_batch_r2 >= pre_batch_r2) {
    flog.warn("ComBat-seq did not reduce PC-batch association -- inspect batch design further.")
  }

  adjusted_counts
}


# --------------------------------------------------------------------------------------------------
# 04. DIFFERENTIAL EXPRESSION (DESeq2, LRT, apeglm SHRINKAGE)
# --------------------------------------------------------------------------------------------------
#' Fit DESeq2 with a likelihood ratio test across the full dose-response design
run_deseq2_lrt <- function(adjusted_counts, pheno_data, num_cores, output_dir = NULL, alpha = 0.05) {
  flog.info("Fitting DESeq2 LRT model (full: ~Treatment_Group, reduced: ~1).")

  dds <- DESeqDataSetFromMatrix(
    countData = adjusted_counts,
    colData = pheno_data,
    design = ~ Treatment_Group
  )
  dds <- dds[rowSums(counts(dds) >= 10) >= 0.2 * ncol(dds), ]  # low-count gene filtering
  flog.info("%d genes retained after low-count filtering.", nrow(dds))
  dds <- DESeq(dds, test = "LRT", reduced = ~1, parallel = (num_cores > 1))

  res <- results(dds, alpha = alpha, independentFiltering = TRUE)

  # apeglm shrinkage for the High_Dose vs Vehicle contrast, used downstream for effect-size ranking
  shrunk <- tryCatch({
    lfcShrink(dds, coef = "Treatment_Group_High_Dose_vs_Vehicle", type = "apeglm")
  }, error = function(e) {
    flog.warn("apeglm shrinkage failed (%s); falling back to unshrunk LFCs.", conditionMessage(e))
    res
  })

  n_sig <- sum(res$padj < alpha, na.rm = TRUE)
  flog.info("DESeq2 LRT complete: %d genes significant at padj < %.2f", n_sig, alpha)

  # vst() subsamples 1000 genes; fall back to the full transform for small gene sets.
  vsd <- if (nrow(dds) >= 1000) vst(dds, blind = FALSE) else varianceStabilizingTransformation(dds, blind = FALSE)

  if (!is.null(output_dir)) {
    de_table <- as.data.frame(res)
    de_table$shrunk_log2FC_high_vs_vehicle <- shrunk$log2FoldChange[match(rownames(de_table), rownames(shrunk))]
    de_table <- tibble::rownames_to_column(de_table, "gene_id")
    write_csv(de_table[order(de_table$padj), ], file.path(output_dir, "deseq2_lrt_results.csv"))
  }

  list(dds = dds, results = res, shrunk_lfc = shrunk, vst_matrix = assay(vsd))
}


# --------------------------------------------------------------------------------------------------
# 05. NETWORK BIOLOGY (WGCNA)
# --------------------------------------------------------------------------------------------------
#' Build a signed weighted gene co-expression network and identify trait-correlated modules
run_wgcna_analysis <- function(vst_matrix, pheno_data, num_cores, n_top_genes = 5000) {
  n_top_genes <- min(n_top_genes, nrow(vst_matrix))
  flog.info("Selecting top %d variable genes for WGCNA network construction.", n_top_genes)

  gene_vars <- matrixStats::rowVars(vst_matrix)
  top_genes <- head(order(gene_vars, decreasing = TRUE), n_top_genes)
  datExpr <- t(vst_matrix[top_genes, ])

  gsg <- goodSamplesGenes(datExpr, verbose = 0)
  if (!gsg$allOK) {
    flog.warn("Removing %d low-quality genes / %d samples flagged by goodSamplesGenes().",
              sum(!gsg$goodGenes), sum(!gsg$goodSamples))
    datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes]
  }

  powers <- c(1:10, seq(12, 20, by = 2))
  sft <- pickSoftThreshold(datExpr, powerVector = powers, verbose = 0, networkType = "signed")
  optimal_power <- sft$powerEstimate
  if (is.na(optimal_power)) {
    flog.warn("No power reached the scale-free R^2 threshold; defaulting to power=6.")
    optimal_power <- 6
  }
  flog.info("Soft-thresholding power selected: %d (scale-free fit R^2 = %.2f)",
            as.integer(optimal_power),
            sft$fitIndices$SFT.R.sq[sft$fitIndices$Power == optimal_power])

  net <- blockwiseModules(
    datExpr, power = optimal_power,
    networkType = "signed", TOMType = "signed", minModuleSize = 50,
    maxBlockSize = ncol(datExpr),
    reassignThreshold = 0, mergeCutHeight = 0.25,
    numericLabels = TRUE, pamRespectsDendro = FALSE,
    saveTOMs = FALSE, verbose = 0, nThreads = num_cores
  )

  module_eigengenes <- net$MEs
  rownames(module_eigengenes) <- rownames(datExpr)
  module_labels <- net$colors

  # Module-trait correlation against continuous toxicity score
  toxicity_aligned <- pheno_data$Toxicity_Score[match(rownames(datExpr), pheno_data$Sample_ID)]
  module_trait_cor <- cor(module_eigengenes, toxicity_aligned, use = "pairwise.complete.obs")
  module_trait_p <- corPvalueStudent(module_trait_cor, nrow(datExpr))

  # Hub gene identification via intramodular connectivity (module 0 = unassigned, excluded)
  connectivity <- intramodularConnectivity(
    adjacency(datExpr, power = optimal_power, type = "signed"), module_labels
  )
  assigned <- connectivity[module_labels != 0, , drop = FALSE]
  hub_genes <- head(rownames(assigned)[order(-assigned$kWithin)], 20)

  flog.info("WGCNA complete: %d modules (incl. unassigned ME0) detected; top module-trait |r| = %.2f",
            ncol(module_eigengenes), max(abs(module_trait_cor)))

  list(
    module_eigengenes = module_eigengenes,
    module_labels = module_labels,
    module_trait_cor = module_trait_cor,
    module_trait_p = module_trait_p,
    hub_genes = hub_genes,
    soft_power = optimal_power,
    sample_ids = rownames(datExpr)
  )
}


# --------------------------------------------------------------------------------------------------
# 06. PATHWAY-LEVEL ENRICHMENT (GSVA) -- used as an engineered feature block, not decorative import
# --------------------------------------------------------------------------------------------------
#' Build an illustrative hepatotoxicity gene-set panel for the simulated data.
#'
#' Each "mechanistic" set is up to ~70% true members of the corresponding injected module plus ~30%
#' random decoys (mimicking an imperfect curated annotation). Only gene IDs are used -- no sample
#' labels -- so this does not leak the outcome. In production replace with
#' `msigdbr::msigdbr(species = "Homo sapiens", category = "H")`.
#' @param module_genes list(oxidative_stress, apoptosis, necrosis) of gene IDs (optional)
build_illustrative_gene_sets <- function(all_genes, module_genes = NULL, seed = 1) {
  set.seed(seed)
  make_set <- function(true_genes, size) {
    true_genes <- intersect(true_genes, all_genes)
    true_genes <- if (length(true_genes) > 0)
      sample(true_genes, min(floor(0.7 * size), length(true_genes))) else character(0)
    decoys <- sample(setdiff(all_genes, true_genes), size - length(true_genes))
    c(true_genes, decoys)
  }
  gm <- module_genes
  size <- function(n) min(n, length(all_genes))
  list(
    HALLMARK_OXIDATIVE_STRESS_LIKE = make_set(gm$oxidative_stress, size(150)),
    HALLMARK_APOPTOSIS_LIKE        = make_set(gm$apoptosis, size(150)),
    HALLMARK_NECROSIS_LIKE         = make_set(gm$necrosis, size(120)),
    HALLMARK_XENOBIOTIC_METAB_LIKE = make_set(NULL, size(100))
  )
}

#' Compute single-sample GSVA enrichment scores for a hepatotoxicity gene-set panel
run_gsva_enrichment <- function(vst_matrix, num_cores, gene_sets = NULL) {
  flog.info("Computing GSVA single-sample enrichment scores for curated hepatotoxicity gene sets.")

  if (is.null(gene_sets)) gene_sets <- build_illustrative_gene_sets(rownames(vst_matrix))
  gene_sets <- lapply(gene_sets, intersect, rownames(vst_matrix))

  gsva_param <- GSVA::gsvaParam(as.matrix(vst_matrix), gene_sets, kcdf = "Gaussian")
  gsva_scores <- GSVA::gsva(gsva_param, BPPARAM = make_bpparam(num_cores), verbose = FALSE)

  flog.info("GSVA complete: %d pathway-level features generated per sample.", nrow(gsva_scores))
  t(gsva_scores)
}


# --------------------------------------------------------------------------------------------------
# 07. REPRESENTATION LEARNING -- DEEP AUTOENCODER (ACTUALLY TRAINED)
# --------------------------------------------------------------------------------------------------
#' Min-max scale each gene (column) to [0, 1]; constant genes become 0.
minmax_scale_columns <- function(x) {
  mins <- apply(x, 2, min)
  rngs <- apply(x, 2, max) - mins
  rngs[rngs == 0] <- 1
  sweep(sweep(x, 2, mins, "-"), 2, rngs, "/")
}

#' Train a deep autoencoder and return the encoder-side latent representation
#'
#' Unlike a placeholder implementation, this function performs a real fit() call and a real
#' predict() call. If TensorFlow/Keras are not installed, it fails explicitly rather than
#' fabricating a latent space with random noise. The autoencoder is unsupervised (it never sees
#' the toxicity label) and is fit on all samples.
build_and_train_autoencoder <- function(vst_matrix, latent_dim = 32, output_dir, seed = 2026,
                                         epochs = 100, batch_size = 32, val_split = 0.2) {
  n_genes <- nrow(vst_matrix)
  flog.info("Building autoencoder: input_dim=%d -> 512 -> 128 -> %d (bottleneck)",
            n_genes, latent_dim)

  if (!reticulate::py_module_available("tensorflow")) {
    stop("TensorFlow backend not available in this environment. Install it first, e.g. ",
         "`keras::install_keras()` (keras 2.x) or `keras3::install_keras()` (keras3), then re-run.")
  }
  tensorflow::set_random_seed(seed)

  x_train <- minmax_scale_columns(t(vst_matrix))  # samples x genes, each gene scaled to 0-1

  input_layer <- layer_input(shape = c(n_genes))
  encoded <- input_layer %>%
    layer_dense(units = 512, activation = "relu") %>%
    layer_batch_normalization() %>%
    layer_dropout(rate = 0.2) %>%
    layer_dense(units = 128, activation = "relu") %>%
    layer_dense(units = latent_dim, activation = "linear", name = "bottleneck_latent_space")

  decoded <- encoded %>%
    layer_dense(units = 128, activation = "relu") %>%
    layer_dense(units = 512, activation = "relu") %>%
    layer_dense(units = n_genes, activation = "sigmoid")

  autoencoder <- keras_model(inputs = input_layer, outputs = decoded)
  autoencoder %>% compile(optimizer = optimizer_adam(learning_rate = 0.001), loss = "mse")

  checkpoint_path <- file.path(output_dir, "autoencoder_best.keras")
  callbacks <- list(
    callback_early_stopping(monitor = "val_loss", patience = 10, restore_best_weights = TRUE),
    callback_model_checkpoint(checkpoint_path, monitor = "val_loss", save_best_only = TRUE),
    callback_reduce_lr_on_plateau(monitor = "val_loss", factor = 0.5, patience = 5)
  )

  flog.info("Training autoencoder for up to %d epochs (batch_size=%d, val_split=%.2f).",
            epochs, batch_size, val_split)

  history <- autoencoder %>% fit(
    x_train, x_train,
    epochs = epochs,
    batch_size = batch_size,
    validation_split = val_split,
    callbacks = callbacks,
    verbose = 0
  )

  final_val_loss <- min(unlist(history$metrics$val_loss))  # best epoch (weights were restored)
  flog.info("Autoencoder training complete after %d epochs. Best val_loss = %.5f",
            length(history$metrics$val_loss), final_val_loss)

  encoder_model <- keras_model(
    inputs = autoencoder$input,
    outputs = get_layer(autoencoder, "bottleneck_latent_space")$output
  )
  latent_features <- predict(encoder_model, x_train, verbose = 0)
  colnames(latent_features) <- paste0("DL_Latent_Dim_", seq_len(latent_dim))
  rownames(latent_features) <- colnames(vst_matrix)

  list(
    latent_features = latent_features,
    encoder = encoder_model,
    history = history,
    final_val_loss = final_val_loss,
    checkpoint_path = checkpoint_path
  )
}


# --------------------------------------------------------------------------------------------------
# 08. PREDICTIVE MODEL -- HOLDOUT + INNER-CV BAYESIAN HYPERPARAMETER OPTIMIZATION (XGBOOST)
# --------------------------------------------------------------------------------------------------
#' Coefficient of determination (1 - SSE/SST); unlike cor()^2 it penalizes biased predictions.
r2_score <- function(pred, actual) {
  1 - sum((actual - pred)^2) / sum((actual - mean(actual))^2)
}

#' Train an XGBoost regressor with an outer holdout split and an inner Bayesian-optimized
#' hyperparameter search. Hyperparameters AND the number of boosting rounds are chosen by
#' 5-fold CV on the training split only; the test split is touched exactly once, for scoring.
train_xgboost_model <- function(ai_dataset, bo_iterations = 15, seed = 2026, num_cores = 1) {
  flog.info("Preparing train/test split and inner CV scoring function for XGBoost.")

  set.seed(seed)
  train_idx <- createDataPartition(ai_dataset$Toxicity, p = 0.8, list = FALSE)
  train_data <- ai_dataset[train_idx, ]
  test_data  <- ai_dataset[-train_idx, ]

  x_train <- as.matrix(select(train_data, -Toxicity))
  y_train <- train_data$Toxicity
  x_test  <- as.matrix(select(test_data, -Toxicity))
  y_test  <- test_data$Toxicity

  dtrain <- xgb.DMatrix(data = x_train, label = y_train)
  dtest  <- xgb.DMatrix(data = x_test, label = y_test)

  cv_folds <- createFolds(y_train, k = 5)  # fixed folds -> comparable scores across BO steps

  make_params <- function(eta, max_depth, subsample, colsample_bytree, min_child_weight) {
    list(
      objective = "reg:squarederror",
      eta = eta,
      max_depth = as.integer(round(max_depth)),
      subsample = subsample,
      colsample_bytree = colsample_bytree,
      min_child_weight = min_child_weight,
      nthread = num_cores
    )
  }

  run_cv <- function(params) {
    cv <- xgb.cv(
      params = params, data = dtrain, nrounds = 1000, folds = cv_folds,
      early_stopping_rounds = 20, verbose = 0, metrics = "rmse"
    )
    log <- as.data.frame(cv$evaluation_log)
    best <- which.min(log$test_rmse_mean)
    list(rmse = log$test_rmse_mean[best], nrounds = best)
  }

  cv_score <- function(eta, max_depth, subsample, colsample_bytree, min_child_weight) {
    cv <- run_cv(make_params(eta, max_depth, subsample, colsample_bytree, min_child_weight))
    list(Score = -cv$rmse)  # ParBayesianOptimization maximizes
  }

  flog.info("Running Bayesian hyperparameter optimization (%d iterations).", bo_iterations)
  bounds <- list(
    eta = c(0.01, 0.3),
    max_depth = c(3, 10),
    subsample = c(0.5, 1.0),
    colsample_bytree = c(0.5, 1.0),
    min_child_weight = c(1, 10)
  )

  set.seed(seed)
  opt_result <- bayesOpt(
    FUN = cv_score, bounds = bounds,
    initPoints = max(6, length(bounds) + 2),
    iters.n = bo_iterations, verbose = 0
  )

  best_params <- getBestPars(opt_result)
  flog.info("Best hyperparameters found: %s",
            paste(sprintf("%s=%.4f", names(best_params), unlist(best_params)), collapse = ", "))

  final_params <- do.call(make_params, as.list(best_params[names(bounds)]))
  final_cv <- run_cv(final_params)
  flog.info("Inner-CV RMSE at best params: %.2f (%d boosting rounds).", final_cv$rmse, final_cv$nrounds)

  xgb_model <- xgb.train(params = final_params, data = dtrain, nrounds = final_cv$nrounds, verbose = 0)

  preds <- predict(xgb_model, dtest)
  rmse <- sqrt(mean((preds - y_test)^2))
  r_squared <- r2_score(preds, y_test)

  flog.info("Final XGBoost model -- Test RMSE: %.2f | Test R^2: %.3f", rmse, r_squared)

  list(
    model = xgb_model,
    best_params = final_params,
    cv_rmse = final_cv$rmse,
    optimization_trace = opt_result,
    test_predictions = preds,
    test_actual = y_test,
    metrics = list(rmse = rmse, r_squared = r_squared, cv_rmse = final_cv$rmse),
    train_data = train_data,
    test_data = test_data
  )
}


# --------------------------------------------------------------------------------------------------
# 09. MODEL INTERPRETATION -- SHAP DECOMPOSITION
# --------------------------------------------------------------------------------------------------
#' Decompose model predictions into per-feature SHAP contributions (replaces raw gain/cover
#' importance with a game-theoretically grounded, per-sample attribution).
compute_shap_values <- function(xgb_model, train_data) {
  flog.info("Computing SHAP values for model interpretation.")

  x_train <- as.matrix(select(train_data, -Toxicity))
  shap_values <- shap.values(xgb_model = xgb_model, X_train = x_train)
  shap_long <- shap.prep(shap_contrib = shap_values$shap_score, X_train = x_train)

  top_features <- shap_values$mean_shap_score %>%
    sort(decreasing = TRUE) %>%
    head(10)

  flog.info("Top %d SHAP-ranked features: %s", min(3, length(top_features)),
            paste(head(names(top_features), 3), collapse = ", "))

  list(shap_values = shap_values, shap_long = shap_long, top_features = top_features)
}


# --------------------------------------------------------------------------------------------------
# 10. STATISTICAL DIAGNOSTICS -- BOOTSTRAP CI, RESIDUALS, CALIBRATION
# --------------------------------------------------------------------------------------------------
#' Bootstrap a confidence interval for test-set RMSE and R^2, and compute residual diagnostics
bootstrap_model_diagnostics <- function(preds, actual, n_boot = 2000, seed = 2026) {
  set.seed(seed)
  n <- length(actual)
  rmse_boot <- numeric(n_boot)
  r2_boot <- numeric(n_boot)

  for (b in seq_len(n_boot)) {
    idx <- sample(seq_len(n), n, replace = TRUE)
    rmse_boot[b] <- sqrt(mean((preds[idx] - actual[idx])^2))
    r2_boot[b] <- r2_score(preds[idx], actual[idx])
  }

  residuals <- preds - actual
  shapiro_p <- if (n >= 3 && n <= 5000) shapiro.test(residuals)$p.value else NA

  list(
    rmse_ci = quantile(rmse_boot, c(0.025, 0.5, 0.975)),
    r2_ci = quantile(r2_boot, c(0.025, 0.5, 0.975), na.rm = TRUE),
    residual_normality_p = shapiro_p,
    residual_mean = mean(residuals),
    residual_sd = sd(residuals)
  )
}


# --------------------------------------------------------------------------------------------------
# 11. PERSISTENCE -- ARTIFACT SERIALIZATION
# --------------------------------------------------------------------------------------------------
#' Write human-readable results (metrics, predictions, hub genes, plots) next to the models.
write_reports <- function(output_dir, xgb_result, shap_result, wgcna_result, diagnostics) {
  metrics <- data.frame(
    metric = c("test_rmse", "test_r2", "inner_cv_rmse", "rmse_ci_low", "rmse_ci_high",
               "r2_ci_low", "r2_ci_high", "residual_mean", "residual_sd", "shapiro_p"),
    value = c(xgb_result$metrics$rmse, xgb_result$metrics$r_squared, xgb_result$metrics$cv_rmse,
              diagnostics$rmse_ci[[1]], diagnostics$rmse_ci[[3]],
              diagnostics$r2_ci[[1]], diagnostics$r2_ci[[3]],
              diagnostics$residual_mean, diagnostics$residual_sd, diagnostics$residual_normality_p)
  )
  write_csv(metrics, file.path(output_dir, "model_metrics.csv"))
  write_csv(data.frame(actual = xgb_result$test_actual, predicted = xgb_result$test_predictions),
            file.path(output_dir, "test_predictions.csv"))
  write_csv(data.frame(rank = seq_along(wgcna_result$hub_genes), gene_id = wgcna_result$hub_genes),
            file.path(output_dir, "wgcna_hub_genes.csv"))
  write_csv(tibble::tibble(feature = names(shap_result$top_features),
                           mean_abs_shap = as.numeric(shap_result$top_features)),
            file.path(output_dir, "shap_top_features.csv"))

  tryCatch({
    p <- ggplot(data.frame(actual = xgb_result$test_actual, predicted = xgb_result$test_predictions),
                aes(actual, predicted)) +
      geom_point(alpha = 0.6) + geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
      labs(title = "Held-out test set: predicted vs. actual toxicity score",
           x = "Actual Toxicity_Score", y = "Predicted Toxicity_Score") +
      theme_minimal()
    ggsave(file.path(output_dir, "predicted_vs_actual.png"), p, width = 6, height = 5, dpi = 150)
    ggsave(file.path(output_dir, "shap_summary.png"),
           SHAPforxgboost::shap.plot.summary(shap_result$shap_long), width = 7, height = 6, dpi = 150)
  }, error = function(e) flog.warn("Plot generation failed: %s", conditionMessage(e)))
}

persist_pipeline_artifacts <- function(output_dir, xgb_result, ae_result, wgcna_result, shap_result) {
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  xgb.save(xgb_result$model, file.path(output_dir, "xgboost_dili_model.xgb"))
  saveRDS(wgcna_result, file.path(output_dir, "wgcna_result.rds"))
  saveRDS(shap_result, file.path(output_dir, "shap_result.rds"))
  saveRDS(ae_result$latent_features, file.path(output_dir, "autoencoder_latent_features.rds"))
  tryCatch(
    ae_result$encoder$save(file.path(output_dir, "autoencoder_encoder.keras")),
    error = function(e) flog.warn("Could not save the encoder model: %s", conditionMessage(e))
  )

  flog.info("All model artifacts persisted to %s", output_dir)
}


# --------------------------------------------------------------------------------------------------
# 12. UNIT TESTS (testthat) -- run inline against the simulation & preprocessing layer
# --------------------------------------------------------------------------------------------------
run_unit_tests <- function() {
  flog.info("Running unit-test suite for simulation and preprocessing functions.")

  reporter <- testthat::ListReporter$new()
  testthat::with_reporter(reporter, {

    test_that("simulate_multiomics_cohort returns correctly shaped objects", {
      sim <- simulate_multiomics_cohort(n_samples = 20, n_genes = 500, seed = 1)
      expect_equal(dim(sim$counts), c(500, 20))
      expect_equal(nrow(sim$pheno), 20)
      expect_identical(rownames(sim$pheno), colnames(sim$counts))
      expect_true(all(sim$pheno$Toxicity_Score >= 0 & sim$pheno$Toxicity_Score <= 100))
      expect_true(all(sapply(sim$ground_truth_modules, length) > 0))
    })

    test_that("ground-truth modules are disjoint and within range", {
      sim <- simulate_multiomics_cohort(n_samples = 20, n_genes = 2000, seed = 4)
      all_idx <- unlist(sim$ground_truth_modules)
      expect_equal(anyDuplicated(all_idx), 0)
      expect_true(all(all_idx >= 1 & all_idx <= 2000))
    })

    test_that("simulated counts are non-negative integers", {
      sim <- simulate_multiomics_cohort(n_samples = 10, n_genes = 200, seed = 2)
      expect_true(all(sim$counts >= 0))
      expect_true(is.integer(sim$counts))
    })

    test_that("simulation is reproducible for a fixed seed", {
      a <- simulate_multiomics_cohort(n_samples = 10, n_genes = 100, seed = 7)
      b <- simulate_multiomics_cohort(n_samples = 10, n_genes = 100, seed = 7)
      expect_identical(a$counts, b$counts)
    })

    test_that("high-dose samples show elevated toxicity relative to vehicle", {
      sim <- simulate_multiomics_cohort(n_samples = 200, n_genes = 300, seed = 3)
      veh_mean <- mean(sim$pheno$Toxicity_Score[sim$pheno$Treatment_Group == "Vehicle"])
      high_mean <- mean(sim$pheno$Toxicity_Score[sim$pheno$Treatment_Group == "High_Dose"])
      expect_gt(high_mean, veh_mean)
    })

    test_that("the simulated batch effect is detectable (and absent when disabled)", {
      with_b <- simulate_multiomics_cohort(n_samples = 90, n_genes = 1000, seed = 5, batch_effect_sd = 0.8)
      no_b   <- simulate_multiomics_cohort(n_samples = 90, n_genes = 1000, seed = 5, batch_effect_sd = 0)
      expect_gt(batch_pc_r2(with_b$counts, with_b$pheno$Batch),
                batch_pc_r2(no_b$counts, no_b$pheno$Batch))
    })

    test_that("autoencoder input scaling maps every gene into [0, 1]", {
      x <- matrix(c(1, 5, 3, 2, 2, 2), nrow = 3)  # second gene is constant
      scaled <- minmax_scale_columns(x)
      expect_equal(dim(scaled), dim(x))
      expect_true(all(scaled >= 0 & scaled <= 1))
      expect_equal(scaled[, 2], c(0, 0, 0))
    })

    test_that("r2_score and bootstrap diagnostics behave sensibly", {
      y <- c(10, 20, 30, 40, 50)
      expect_equal(r2_score(y, y), 1)
      expect_lt(r2_score(rep(mean(y), 5), y), 1e-12)
      d <- bootstrap_model_diagnostics(y + c(1, -1, 1, -1, 0), y, n_boot = 50, seed = 1)
      expect_lte(d$rmse_ci[[1]], d$rmse_ci[[3]])
    })

    test_that("illustrative gene sets are enriched for the injected modules", {
      sim <- simulate_multiomics_cohort(n_samples = 10, n_genes = 2000, seed = 6)
      module_genes <- lapply(sim$ground_truth_modules, function(idx) rownames(sim$counts)[idx])
      gs <- build_illustrative_gene_sets(rownames(sim$counts), module_genes)
      truth <- module_genes$oxidative_stress
      # far above the chance rate of length(truth) / n_genes
      expect_gt(mean(gs$HALLMARK_OXIDATIVE_STRESS_LIKE %in% truth),
                5 * length(truth) / nrow(sim$counts))
      expect_true(all(unlist(gs) %in% rownames(sim$counts)))
    })

    test_that("CLI parsing applies defaults and overrides", {
      o <- parse_cli_args(c("--n_samples", "60", "--n_genes", "3000"))
      expect_equal(o$n_samples, 60)
      expect_equal(o$seed, 2026)
      expect_error(parse_cli_args(c("--n_samples", "5")))
    })
  })

  results <- as.data.frame(reporter$get_results())
  bad <- results[results$failed > 0 | results$error, , drop = FALSE]
  if (nrow(bad) > 0) {
    stop("Unit tests failed: ", paste(bad$test, collapse = "; "))
  }
  flog.info("Unit-test suite passed (%d test blocks).", nrow(results))
}


# --------------------------------------------------------------------------------------------------
# 13. MAIN ORCHESTRATION
# --------------------------------------------------------------------------------------------------
run_pipeline <- function(opts) {
  init_logger(opts$output_dir)
  num_cores <- init_parallel_backend()
  capture_provenance(opts$output_dir, opts)

  if (isTRUE(opts$run_tests)) run_unit_tests()

  flog.info("=== STAGE 1/7: Data simulation ===")
  sim <- simulate_multiomics_cohort(opts$n_samples, opts$n_genes, opts$seed)

  flog.info("=== STAGE 2/7: Batch correction ===")
  adjusted_counts <- run_batch_correction(sim$counts, sim$pheno, opts$output_dir)

  flog.info("=== STAGE 3/7: Differential expression (DESeq2 LRT) ===")
  de_result <- run_deseq2_lrt(adjusted_counts, sim$pheno, num_cores, opts$output_dir)

  flog.info("=== STAGE 4/7: Network biology (WGCNA) ===")
  wgcna_result <- run_wgcna_analysis(de_result$vst_matrix, sim$pheno, num_cores)

  flog.info("=== STAGE 5/7: Pathway enrichment (GSVA) ===")
  module_genes <- lapply(sim$ground_truth_modules, function(idx) rownames(sim$counts)[idx])
  gene_sets <- build_illustrative_gene_sets(rownames(de_result$vst_matrix), module_genes,
                                            seed = opts$seed)
  gsva_features <- run_gsva_enrichment(de_result$vst_matrix, num_cores, gene_sets)

  flog.info("=== STAGE 6/7: Representation learning (autoencoder) ===")
  ae_result <- build_and_train_autoencoder(de_result$vst_matrix, opts$latent_dim, opts$output_dir,
                                           seed = opts$seed)

  flog.info("=== STAGE 7/7: Predictive modeling (XGBoost + Bayesian tuning) ===")
  common_samples <- Reduce(intersect, list(
    rownames(ae_result$latent_features),
    rownames(wgcna_result$module_eigengenes),
    rownames(gsva_features),
    sim$pheno$Sample_ID
  ))

  ai_dataset <- data.frame(
    Toxicity = sim$pheno$Toxicity_Score[match(common_samples, sim$pheno$Sample_ID)],
    ae_result$latent_features[common_samples, , drop = FALSE],
    wgcna_result$module_eigengenes[common_samples, , drop = FALSE],
    gsva_features[common_samples, , drop = FALSE]
  )

  xgb_result <- train_xgboost_model(ai_dataset, opts$bo_iterations, opts$seed, num_cores)
  shap_result <- compute_shap_values(xgb_result$model, xgb_result$train_data)
  diagnostics <- bootstrap_model_diagnostics(xgb_result$test_predictions, xgb_result$test_actual,
                                             seed = opts$seed)

  flog.info("Bootstrap RMSE 95%% CI: [%.2f, %.2f] | R^2 95%% CI: [%.3f, %.3f]",
            diagnostics$rmse_ci[1], diagnostics$rmse_ci[3],
            diagnostics$r2_ci[1], diagnostics$r2_ci[3])

  persist_pipeline_artifacts(opts$output_dir, xgb_result, ae_result, wgcna_result, shap_result)
  write_reports(opts$output_dir, xgb_result, shap_result, wgcna_result, diagnostics)

  flog.info("Pipeline execution complete. Artifacts written to %s", opts$output_dir)

  invisible(list(
    simulation = sim,
    differential_expression = de_result,
    wgcna = wgcna_result,
    gsva = gsva_features,
    autoencoder = ae_result,
    xgboost = xgb_result,
    shap = shap_result,
    diagnostics = diagnostics
  ))
}

# Guard so the file can be `source()`-d for interactive/unit-test use, or run directly via Rscript.
if (identical(environment(), globalenv()) && sys.nframe() == 0) {
  opts <- parse_cli_args()
  results <- run_pipeline(opts)
}
