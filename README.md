# AI-Integrated Toxicogenomics & Multi-Omics Predictive Engine

**Nested-CV XGBoost + Bayesian tuning, a trained Keras autoencoder, WGCNA network biology, and GSVA pathway enrichment for quantitative Drug-Induced Liver Injury (DILI) risk modeling from RNA-seq.**

## Overview

This repository models Drug-Induced Liver Injury (DILI) severity directly from transcriptomic signatures by combining three complementary feature representations — mechanistic (WGCNA co-expression modules), curated biological (GSVA pathway enrichment), and learned (autoencoder latent space) — into a single predictive layer.

The pipeline is organized as independently testable functions rather than a linear script, so each stage (simulation, batch correction, differential expression, network construction, representation learning, predictive modeling, interpretation) can be run, validated, and reasoned about on its own.

## Pipeline Architecture

| Stage | Method | Purpose |
|---|---|---|
| 1. Simulation | Parameterized NB count generator | Synthetic dose-response cohort (Vehicle / Low / High dose) with injected mechanistic gene modules |
| 2. Batch correction | `ComBat-seq` | Removes multi-lab technical variance; validated with a PCA-based pre/post QC check, not assumed |
| 3. Differential expression | `DESeq2`, LRT (full vs. reduced model), `apeglm` shrinkage | Identifies dose-responsive genes; independent filtering applied |
| 4. Network biology | `WGCNA`, signed TOM, automatic soft-thresholding | Reduces 5,000 variable genes to co-expression modules; module-trait correlation and intramodular-connectivity hub gene identification |
| 5. Pathway enrichment | `GSVA` (Gaussian kernel) | Single-sample enrichment scores over a curated hepatotoxicity gene-set panel, used as model features |
| 6. Representation learning | Keras/TensorFlow autoencoder (512→128→32→128→512) | Real `fit()`/`predict()` calls with early stopping, LR reduction on plateau, and checkpointing — not a placeholder |
| 7. Predictive modeling | `XGBoost`, 5-fold inner CV, Bayesian hyperparameter optimization (`ParBayesianOptimization`) | Tunes `eta`, `max_depth`, `subsample`, `colsample_bytree`, `min_child_weight` against CV RMSE rather than a fixed grid |
| 8. Interpretation | SHAP (`SHAPforxgboost`) | Per-sample, game-theoretic feature attribution in place of raw gain/cover importance |
| 9. Diagnostics | Bootstrap CI (2,000 resamples) | 95% CI on test RMSE/R², residual normality (Shapiro-Wilk), residual mean/SD |
| 10. Persistence | `xgb.save`, `saveRDS`, TF SavedModel | Serializes trained model, WGCNA/SHAP objects, and a resolved-config + session-info provenance snapshot for reproducibility |

## What Changed From the Original Version

The earlier version of this pipeline demonstrated the right technology list (DESeq2, WGCNA, Keras, XGBoost) but had a load-bearing gap: the autoencoder's `fit()`/`predict()` calls were commented out, and the "latent features" fed into the predictive model were `matrix(rnorm(...))` — random noise standing in for a trained representation. This version fixes that directly: the autoencoder is actually trained, with real callbacks and a real bottleneck extraction.

Beyond that fix, this version adds the layers that separate a demo script from a defensible modeling pipeline:

- **Nested validation** instead of a single train/test split with fixed hyperparameters — hyperparameters are chosen by inner 5-fold CV, not read off a hardcoded list.
- **Interpretability** via SHAP rather than raw split-gain importance, which is sensitive to feature scale and correlation structure.
- **Uncertainty quantification** — bootstrap confidence intervals on RMSE/R², not a single point estimate.
- **QC as code** — the PCA-based batch-effect check is computed and logged, not assumed to have worked because `ComBat_seq()` was called.
- **Modularity and unit tests** — the simulation and preprocessing layer has `testthat` coverage, so correctness claims are checked, not asserted.
- **Reproducibility artifacts** — every run writes a resolved config, session info, and a timestamped log.

## Technical Stack

- **Language:** R (≥ 4.2)
- **Differential expression / network biology:** `DESeq2`, `sva`, `WGCNA`, `GSVA`, `GSEABase`
- **Deep learning:** `keras`, `tensorflow`
- **Predictive modeling:** `xgboost`, `caret`, `ParBayesianOptimization`
- **Interpretation:** `SHAPforxgboost`
- **Engineering:** `optparse` (CLI), `futile.logger` (structured logging), `doParallel`, `testthat`, `yaml`

## Usage

```bash
# Default run (500 samples, 20,000 genes, seed 2026)
Rscript ai_toxicogenomics_master_pipeline.R

# Custom run with a YAML config override
Rscript ai_toxicogenomics_master_pipeline.R --config_file config.yaml

# Direct CLI overrides
Rscript ai_toxicogenomics_master_pipeline.R --n_samples 800 --bo_iterations 25 --output_dir results/
```

Sourcing the file interactively (`source("ai_toxicogenomics_master_pipeline.R")`) does **not** auto-execute the pipeline — each stage is a standalone function (`simulate_multiomics_cohort()`, `run_deseq2_lrt()`, `run_wgcna_analysis()`, etc.), so individual stages can be run, inspected, or unit-tested in isolation. Full execution requires calling `run_pipeline(opts)` explicitly or invoking the script via `Rscript`.

## Notes on the Synthetic Data

The RNA-seq cohort and gene-set panel in this repository are simulated (negative-binomial counts with injected dose-dependent modules) so the full pipeline is runnable end-to-end without a data-access agreement. The gene-set panel used for GSVA is illustrative; production use should substitute a curated panel (e.g., MSigDB Hallmark via `msigdbr`) matched to the exposure/tissue context.
