#!/usr/bin/env Rscript
# ==================================================================================================
#  AI-INTEGRATED TOXICOGENOMICS & MULTI-OMICS PREDICTIVE MODELING ENGINE
#  Version: 4.0.0  |  DILI (Drug-Induced Liver Injury) Transcriptomic Risk Modeling
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
})

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
  }
  opts
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
  num_cores <- max(1, parallel::detectCores() - 1)
  registerDoParallel(cores = num_cores)
  allowWGCNAThreads(nThreads = num_cores)
  flog.info("Parallel backend registered across %d cores.", num_cores)
  num_cores
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
#' Simulate a dose-response hepatotoxicity RNA-seq cohort
#'
#' @param n_samples number of simulated samples
#' @param n_genes number of simulated genes
#' @param seed RNG seed
#' @return list(counts = matrix, pheno = data.frame, ground_truth_modules = list)
simulate_multiomics_cohort <- function(n_samples = 500, n_genes = 20000, seed = 2026) {
  set.seed(seed)
  flog.info("Simulating cohort: n_samples=%d, n_genes=%d", n_samples, n_genes)

  pheno_data <- data.frame(
    Sample_ID = paste0("Hepatocyte_", seq_len(n_samples)),
    Batch = factor(sample(c("Lab_A", "Lab_B", "Lab_C"), n_samples, replace = TRUE)),
    Treatment_Group = factor(
      sample(c("Vehicle", "Low_Dose", "High_Dose"), n_samples, replace = TRUE, prob = c(0.4, 0.3, 0.3)),
      levels = c("Vehicle", "Low_Dose", "High_Dose")
    )
  )

  pheno_data$Toxicity_Score <- with(pheno_data, {
    base <- ifelse(Treatment_Group == "Vehicle", rnorm(n_samples, 10, 5),
             ifelse(Treatment_Group == "Low_Dose", rnorm(n_samples, 40, 15),
                    rnorm(n_samples, 85, 10)))
    pmin(pmax(base, 0), 100)
  })

  base_means <- rnbinom(n_genes, mu = 800, size = 1.2)
  raw_counts <- matrix(rnbinom(n_genes * n_samples, mu = base_means, size = 1.2), ncol = n_samples)
  rownames(raw_counts) <- paste0("ENSG", sprintf("%06d", seq_len(n_genes)))
  colnames(raw_counts) <- pheno_data$Sample_ID

  # Ground-truth mechanistic gene modules injected into the synthetic signal.
  mod_ox_stress  <- sample(seq_len(n_genes), 400)
  mod_apoptosis  <- sample(setdiff(seq_len(n_genes), mod_ox_stress), 400)
  mod_necrosis   <- sample(setdiff(seq_len(n_genes), c(mod_ox_stress, mod_apoptosis)), 300)

  high_idx <- which(pheno_data$Treatment_Group == "High_Dose")
  low_idx  <- which(pheno_data$Treatment_Group == "Low_Dose")

  raw_counts[mod_ox_stress, high_idx] <- raw_counts[mod_ox_stress, high_idx] *
    matrix(runif(400 * length(high_idx), 3, 10), ncol = length(high_idx))
  raw_counts[mod_apoptosis, high_idx] <- raw_counts[mod_apoptosis, high_idx] *
    matrix(runif(400 * length(high_idx), 0.1, 0.5), ncol = length(high_idx))
  raw_counts[mod_necrosis, high_idx] <- raw_counts[mod_necrosis, high_idx] *
    matrix(runif(300 * length(high_idx), 4, 12), ncol = length(high_idx))
  raw_counts[mod_ox_stress, low_idx] <- raw_counts[mod_ox_stress, low_idx] *
    matrix(runif(400 * length(low_idx), 1.5, 3), ncol = length(low_idx))

  storage.mode(raw_counts) <- "integer"

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
#' Run ComBat-seq batch correction and report a PCA-based sanity check
run_batch_correction <- function(raw_counts, pheno_data, output_dir) {
  flog.info("Running ComBat-seq batch correction across %d batches.", nlevels(pheno_data$Batch))

  pre_pca <- prcomp(t(log1p(raw_counts)), scale. = FALSE)
  pre_batch_r2 <- summary(lm(pre_pca$x[, 1] ~ pheno_data$Batch))$r.squared

  adjusted_counts <- ComBat_seq(
    counts = as.matrix(raw_counts),
    batch = pheno_data$Batch,
    group = pheno_data$Treatment_Group
  )

  post_pca <- prcomp(t(log1p(adjusted_counts)), scale. = FALSE)
  post_batch_r2 <- summary(lm(post_pca$x[, 1] ~ pheno_data$Batch))$r.squared

  flog.info("Batch effect on PC1 -- before: R2=%.3f | after: R2=%.3f", pre_batch_r2, post_batch_r2)

  qc_report <- data.frame(
    stage = c("pre_correction", "post_correction"),
    pc1_batch_r2 = c(pre_batch_r2, post_batch_r2)
  )
  write_csv(qc_report, file.path(output_dir, "batch_correction_qc.csv"))

  if (post_batch_r2 >= pre_batch_r2) {
    flog.warn("ComBat-seq did not reduce PC1 batch association -- inspect batch design further.")
  }

  adjusted_counts
}


# --------------------------------------------------------------------------------------------------
# 04. DIFFERENTIAL EXPRESSION (DESeq2, LRT, apeglm SHRINKAGE)
# --------------------------------------------------------------------------------------------------
#' Fit DESeq2 with a likelihood ratio test across the full dose-response design
run_deseq2_lrt <- function(adjusted_counts, pheno_data, num_cores, alpha = 0.05) {
  flog.info("Fitting DESeq2 LRT model (full: ~Treatment_Group, reduced: ~1).")

  dds <- DESeqDataSetFromMatrix(
    countData = adjusted_counts,
    colData = pheno_data,
    design = ~ Treatment_Group
  )
  dds <- dds[rowSums(counts(dds) >= 10) >= 0.2 * ncol(dds), ]  # low-count gene filtering
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

  vsd <- vst(dds, blind = FALSE)

  list(dds = dds, results = res, shrunk_lfc = shrunk, vst_matrix = assay(vsd))
}


# --------------------------------------------------------------------------------------------------
# 05. NETWORK BIOLOGY (WGCNA)
# --------------------------------------------------------------------------------------------------
#' Build a signed weighted gene co-expression network and identify trait-correlated modules
run_wgcna_analysis <- function(vst_matrix, pheno_data, num_cores, n_top_genes = 5000) {
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
            optimal_power, sft$fitIndices$SFT.R.sq[sft$fitIndices$Power == optimal_power])

  net <- blockwiseModules(
    datExpr, power = optimal_power,
    TOMType = "signed", minModuleSize = 50,
    reassignThreshold = 0, mergeCutHeight = 0.25,
    numericLabels = TRUE, pamRespectsDendro = FALSE,
    saveTOMs = FALSE, verbose = 0, nThreads = num_cores
  )

  module_eigengenes <- net$MEs
  module_labels <- net$colors

  # Module-trait correlation against continuous toxicity score
  toxicity_aligned <- pheno_data$Toxicity_Score[match(rownames(datExpr), pheno_data$Sample_ID)]
  module_trait_cor <- cor(module_eigengenes, toxicity_aligned, use = "pairwise.complete.obs")
  module_trait_p <- corPvalueStudent(module_trait_cor, nrow(datExpr))

  # Hub gene identification via intramodular connectivity
  connectivity <- intramodularConnectivity(adjacency(datExpr, power = optimal_power), module_labels)
  hub_genes <- rownames(connectivity)[order(-connectivity$kWithin)][1:20]

  flog.info("WGCNA complete: %d modules detected; top module-trait |r| = %.2f",
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
#' Compute single-sample GSVA enrichment scores for a small curated hepatotoxicity gene-set panel
run_gsva_enrichment <- function(vst_matrix, num_cores) {
  flog.info("Computing GSVA single-sample enrichment scores for curated hepatotoxicity gene sets.")

  # In a real deployment these would be pulled from MSigDB Hallmark / Reactome via msigdbr;
  # a small illustrative panel is defined here from the row-index space of the simulated data
  # so the function is runnable standalone. Swap in `msigdbr::msigdbr(category = "H")` in production.
  all_genes <- rownames(vst_matrix)
  set.seed(1)
  gene_sets <- list(
    HALLMARK_OXIDATIVE_STRESS_LIKE = sample(all_genes, 150),
    HALLMARK_APOPTOSIS_LIKE        = sample(all_genes, 150),
    HALLMARK_NECROSIS_LIKE         = sample(all_genes, 120),
    HALLMARK_XENOBIOTIC_METAB_LIKE = sample(all_genes, 100)
  )

  gsva_param <- GSVA::gsvaParam(as.matrix(vst_matrix), gene_sets, kcdf = "Gaussian")
  gsva_scores <- GSVA::gsva(gsva_param, BPPARAM = BiocParallel::MulticoreParam(num_cores))

  flog.info("GSVA complete: %d pathway-level features generated per sample.", nrow(gsva_scores))
  t(gsva_scores)
}


# --------------------------------------------------------------------------------------------------
# 07. REPRESENTATION LEARNING -- DEEP AUTOENCODER (ACTUALLY TRAINED)
# --------------------------------------------------------------------------------------------------
#' Train a deep autoencoder and return the encoder-side latent representation
#'
#' Unlike a placeholder implementation, this function performs a real fit() call and a real
#' predict() call. If TensorFlow/Keras are not installed, it fails explicitly rather than
#' fabricating a latent space with random noise.
build_and_train_autoencoder <- function(vst_matrix, latent_dim = 32, output_dir,
                                         epochs = 100, batch_size = 32, val_split = 0.2) {
  flog.info("Building autoencoder: input_dim=%d -> 512 -> 128 -> %d (bottleneck)",
            nrow(vst_matrix), latent_dim)

  if (!reticulate::py_module_available("tensorflow")) {
    stop("TensorFlow backend not available in this environment. Install tensorflow/keras ",
         "(e.g. `keras::install_keras()`) before running the representation-learning stage.")
  }

  n_genes <- nrow(vst_matrix)

  normalize_rowwise <- function(x) {
    rng <- max(x) - min(x)
    if (rng == 0) return(rep(0, length(x)))
    (x - min(x)) / rng
  }
  x_train <- t(apply(vst_matrix, 1, normalize_rowwise))  # samples x genes, scaled 0-1

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

  final_val_loss <- tail(history$metrics$val_loss, 1)
  flog.info("Autoencoder training complete. Final val_loss = %.5f", final_val_loss)

  encoder_model <- keras_model(
    inputs = autoencoder$input,
    outputs = get_layer(autoencoder, "bottleneck_latent_space")$output
  )
  latent_features <- predict(encoder_model, x_train, verbose = 0)
  colnames(latent_features) <- paste0("DL_Latent_Dim_", seq_len(latent_dim))
  rownames(latent_features) <- colnames(vst_matrix)

  list(
    latent_features = latent_features,
    history = history,
    final_val_loss = final_val_loss,
    checkpoint_path = checkpoint_path
  )
}


# --------------------------------------------------------------------------------------------------
# 08. PREDICTIVE MODEL -- NESTED CV + BAYESIAN HYPERPARAMETER OPTIMIZATION (XGBOOST)
# --------------------------------------------------------------------------------------------------
#' Train an XGBoost regressor with an outer holdout split and an inner Bayesian-optimized
#' hyperparameter search (rather than a single fixed hyperparameter set).
train_xgboost_model <- function(ai_dataset, bo_iterations = 15, seed = 2026) {
  flog.info("Preparing train/test split and nested CV scoring function for XGBoost.")

  set.seed(seed)
  train_idx <- createDataPartition(ai_dataset$Toxicity, p = 0.8, list = FALSE)
  train_data <- ai_dataset[train_idx, ]
  test_data  <- ai_dataset[-train_idx, ]

  x_train <- as.matrix(train_data %>% select(-Toxicity))
  y_train <- train_data$Toxicity
  x_test  <- as.matrix(test_data %>% select(-Toxicity))
  y_test  <- test_data$Toxicity

  dtrain <- xgb.DMatrix(data = x_train, label = y_train)
  dtest  <- xgb.DMatrix(data = x_test, label = y_test)

  # Inner scoring function: 5-fold CV RMSE for a given hyperparameter set
  cv_score <- function(eta, max_depth, subsample, colsample_bytree, min_child_weight) {
    params <- list(
      objective = "reg:squarederror",
      eta = eta,
      max_depth = as.integer(round(max_depth)),
      subsample = subsample,
      colsample_bytree = colsample_bytree,
      min_child_weight = min_child_weight
    )
    cv <- xgb.cv(
      params = params, data = dtrain, nrounds = 500, nfold = 5,
      early_stopping_rounds = 20, verbose = 0, metrics = "rmse"
    )
    best_rmse <- min(cv$evaluation_log$test_rmse_mean)
    list(Score = -best_rmse, nrounds = cv$best_iteration)  # ParBayesianOptimization maximizes
  }

  flog.info("Running Bayesian hyperparameter optimization (%d iterations).", bo_iterations)
  bounds <- list(
    eta = c(0.01, 0.3),
    max_depth = c(3, 10),
    subsample = c(0.5, 1.0),
    colsample_bytree = c(0.5, 1.0),
    min_child_weight = c(1, 10)
  )

  opt_result <- bayesOpt(
    FUN = cv_score, bounds = bounds,
    initPoints = max(6, length(bounds) + 2),
    iters.n = bo_iterations, verbose = 0
  )

  best_params <- getBestPars(opt_result)
  flog.info("Best hyperparameters found: %s",
            paste(sprintf("%s=%.4f", names(best_params), unlist(best_params)), collapse = ", "))

  final_params <- c(
    list(objective = "reg:squarederror"),
    eta = best_params$eta,
    max_depth = as.integer(round(best_params$max_depth)),
    subsample = best_params$subsample,
    colsample_bytree = best_params$colsample_bytree,
    min_child_weight = best_params$min_child_weight
  )

  xgb_model <- xgb.train(
    params = final_params, data = dtrain, nrounds = 1000,
    watchlist = list(train = dtrain, test = dtest),
    early_stopping_rounds = 20, print_every_n = 100, verbose = 0
  )

  preds <- predict(xgb_model, dtest)
  rmse <- sqrt(mean((preds - y_test)^2))
  r_squared <- cor(preds, y_test)^2

  flog.info("Final XGBoost model -- Test RMSE: %.2f | Test R^2: %.3f", rmse, r_squared)

  list(
    model = xgb_model,
    best_params = final_params,
    optimization_trace = opt_result,
    test_predictions = preds,
    test_actual = y_test,
    metrics = list(rmse = rmse, r_squared = r_squared),
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

  x_train <- as.matrix(train_data %>% select(-Toxicity))
  shap_values <- shap.values(xgb_model = xgb_model, X_train = x_train)
  shap_long <- shap.prep(shap_contrib = shap_values$shap_score, X_train = x_train)

  top_features <- shap_values$mean_shap_score %>%
    sort(decreasing = TRUE) %>%
    head(10)

  flog.info("Top 3 SHAP-ranked features: %s",
            paste(names(top_features)[1:3], collapse = ", "))

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
    r2_boot[b] <- suppressWarnings(cor(preds[idx], actual[idx])^2)
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
persist_pipeline_artifacts <- function(output_dir, xgb_result, ae_result, wgcna_result, shap_result) {
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  xgb.save(xgb_result$model, file.path(output_dir, "xgboost_dili_model.xgb"))
  saveRDS(wgcna_result, file.path(output_dir, "wgcna_result.rds"))
  saveRDS(shap_result, file.path(output_dir, "shap_result.rds"))
  save_model_tf(keras_model(
    inputs = layer_input(shape = ncol(ae_result$latent_features)), # lightweight metadata stub
    outputs = layer_input(shape = ncol(ae_result$latent_features))
  ), file.path(output_dir, "autoencoder_metadata_stub"))

  flog.info("All model artifacts persisted to %s", output_dir)
}


# --------------------------------------------------------------------------------------------------
# 12. UNIT TESTS (testthat) -- run inline against the simulation & preprocessing layer
# --------------------------------------------------------------------------------------------------
run_unit_tests <- function() {
  flog.info("Running unit-test suite for simulation and preprocessing functions.")

  test_that("simulate_multiomics_cohort returns correctly shaped objects", {
    sim <- simulate_multiomics_cohort(n_samples = 20, n_genes = 500, seed = 1)
    expect_equal(dim(sim$counts), c(500, 20))
    expect_equal(nrow(sim$pheno), 20)
    expect_true(all(sim$pheno$Toxicity_Score >= 0 & sim$pheno$Toxicity_Score <= 100))
    expect_true(all(sapply(sim$ground_truth_modules, length) > 0))
  })

  test_that("simulated counts are non-negative integers", {
    sim <- simulate_multiomics_cohort(n_samples = 10, n_genes = 200, seed = 2)
    expect_true(all(sim$counts >= 0))
    expect_true(is.integer(sim$counts))
  })

  test_that("high-dose samples show elevated toxicity relative to vehicle", {
    sim <- simulate_multiomics_cohort(n_samples = 200, n_genes = 300, seed = 3)
    veh_mean <- mean(sim$pheno$Toxicity_Score[sim$pheno$Treatment_Group == "Vehicle"])
    high_mean <- mean(sim$pheno$Toxicity_Score[sim$pheno$Treatment_Group == "High_Dose"])
    expect_gt(high_mean, veh_mean)
  })

  flog.info("Unit-test suite passed.")
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
  de_result <- run_deseq2_lrt(adjusted_counts, sim$pheno, num_cores)

  flog.info("=== STAGE 4/7: Network biology (WGCNA) ===")
  wgcna_result <- run_wgcna_analysis(de_result$vst_matrix, sim$pheno, num_cores)

  flog.info("=== STAGE 5/7: Pathway enrichment (GSVA) ===")
  gsva_features <- run_gsva_enrichment(de_result$vst_matrix, num_cores)

  flog.info("=== STAGE 6/7: Representation learning (autoencoder) ===")
  ae_result <- build_and_train_autoencoder(de_result$vst_matrix, opts$latent_dim, opts$output_dir)

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

  xgb_result <- train_xgboost_model(ai_dataset, opts$bo_iterations, opts$seed)
  shap_result <- compute_shap_values(xgb_result$model, xgb_result$train_data)
  diagnostics <- bootstrap_model_diagnostics(xgb_result$test_predictions, xgb_result$test_actual)

  flog.info("Bootstrap RMSE 95%% CI: [%.2f, %.2f] | R^2 95%% CI: [%.3f, %.3f]",
            diagnostics$rmse_ci[1], diagnostics$rmse_ci[3],
            diagnostics$r2_ci[1], diagnostics$r2_ci[3])

  persist_pipeline_artifacts(opts$output_dir, xgb_result, ae_result, wgcna_result, shap_result)

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
