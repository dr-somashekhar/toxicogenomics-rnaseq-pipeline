# AI-Integrated Toxicogenomics & Multi-Omics Predictive Engine

**Nested-CV XGBoost + Bayesian tuning, a trained Keras autoencoder, WGCNA network biology, and GSVA pathway enrichment for quantitative Drug-Induced Liver Injury (DILI) risk modeling from RNA-seq.**

## Overview

This repository models Drug-Induced Liver Injury (DILI) severity directly from transcriptomic signatures by combining three complementary feature representations — mechanistic (WGCNA co-expression modules), curated biological (GSVA pathway enrichment), and learned (autoencoder latent space) — into a single predictive layer.

The pipeline is organized as independently testable functions rather than a linear script, so each stage (simulation, batch correction, differential expression, network construction, representation learning, predictive modeling, interpretation) can be run, validated, and reasoned about on its own.

## Pipeline Architecture

| Stage | Method | Purpose |
|---|---|---|
| 1. Simulation | Parameterized NB count generator | Synthetic dose-response cohort (Vehicle / Low / High dose) with injected mechanistic gene modules |
| 2. Batch correction | `ComBat-seq` | Removes the simulated multi-lab technical variance; validated with a PCA-based pre/post QC check (max batch R² over PC1–5), not assumed |
| 3. Differential expression | `DESeq2`, LRT (full vs. reduced model), `apeglm` shrinkage | Identifies dose-responsive genes; independent filtering applied |
| 4. Network biology | `WGCNA`, signed TOM, automatic soft-thresholding | Reduces 5,000 variable genes to co-expression modules; module-trait correlation and intramodular-connectivity hub gene identification |
| 5. Pathway enrichment | `GSVA` (Gaussian kernel) | Single-sample enrichment scores over a curated hepatotoxicity gene-set panel, used as model features |
| 6. Representation learning | Keras/TensorFlow autoencoder (512→128→32→128→512) | Real `fit()`/`predict()` calls with early stopping, LR reduction on plateau, and checkpointing — not a placeholder |
| 7. Predictive modeling | `XGBoost`, 80/20 holdout + 5-fold inner CV, Bayesian hyperparameter optimization (`ParBayesianOptimization`) | Tunes `eta`, `max_depth`, `subsample`, `colsample_bytree`, `min_child_weight` and the number of boosting rounds against CV RMSE on the training split only; the test split is scored once |
| 8. Interpretation | SHAP (`SHAPforxgboost`) | Per-sample, game-theoretic feature attribution in place of raw gain/cover importance |
| 9. Diagnostics | Bootstrap CI (2,000 resamples) | 95% CI on test RMSE and R² (1 − SSE/SST), residual normality (Shapiro-Wilk), residual mean/SD |
| 10. Persistence | `xgb.save`, `saveRDS`, TF SavedModel | Serializes the trained model, WGCNA/SHAP objects, the encoder, CSV/PNG reports, and a resolved-config + session-info provenance snapshot for reproducibility |

## What Changed From the Original Version

The earlier version of this pipeline demonstrated the right technology list (DESeq2, WGCNA, Keras, XGBoost) but had a load-bearing gap: the autoencoder's `fit()`/`predict()` calls were commented out, and the "latent features" fed into the predictive model were `matrix(rnorm(...))` — random noise standing in for a trained representation. This version fixes that directly: the autoencoder is actually trained, with real callbacks and a real bottleneck extraction.

Beyond that fix, this version adds the layers that separate a demo script from a defensible modeling pipeline:

- **Leak-free model selection** — hyperparameters and the number of boosting rounds are chosen by inner 5-fold CV on the training split; the held-out test split is never used for early stopping or tuning.
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
- **Engineering:** `BiocParallel`, `optparse` (CLI), `futile.logger` (structured logging), `doParallel`, `testthat`, `yaml`

## Running it on your desktop

1. **Install R ≥ 4.2** from <https://cran.r-project.org> (Windows users also need [Rtools](https://cran.r-project.org/bin/windows/Rtools/); macOS users need the Xcode command-line tools: `xcode-select --install`). Python 3.9–3.12 must be available for TensorFlow (reticulate can fetch one).
2. **Clone the repo:** `git clone https://github.com/dr-somashekhar/toxicogenomics-rnaseq-pipeline.git && cd toxicogenomics-rnaseq-pipeline`
3. **Install dependencies (once, ~20–40 min):** `Rscript install_dependencies.R`
4. **Smoke test (a few minutes):** `Rscript ai_toxicogenomics_master_pipeline.R --config_file config.example.yaml`
5. **Full run (defaults: 500 samples × 20,000 genes; needs ~16 GB RAM, 1–3 h):** `Rscript ai_toxicogenomics_master_pipeline.R`
6. **Inspect results** in the output directory (`results_demo/` or `pipeline_artifacts/`) — see *Outputs* below.

## Outputs

| File | Content |
|---|---|
| `pipeline_run_*.log` | Timestamped audit trail |
| `resolved_config.yaml`, `provenance.rds` | Exact settings + `sessionInfo()` |
| `batch_correction_qc.csv` | Batch R² on PC1–5 before/after ComBat-seq |
| `deseq2_lrt_results.csv` | LRT results per gene with apeglm-shrunk log2FC (High vs Vehicle) |
| `wgcna_hub_genes.csv`, `wgcna_result.rds` | Hub genes (unassigned module excluded) and full network results |
| `model_metrics.csv` | Test RMSE/R², inner-CV RMSE, bootstrap CIs, residual diagnostics |
| `test_predictions.csv`, `predicted_vs_actual.png` | Held-out predictions and plot |
| `shap_top_features.csv`, `shap_summary.png`, `shap_result.rds` | SHAP interpretation |
| `xgboost_dili_model.xgb`, `autoencoder_encoder.keras`, `autoencoder_best.keras` | Trained models |

## Usage

```bash
# Default run (500 samples, 20,000 genes, seed 2026)
Rscript ai_toxicogenomics_master_pipeline.R

# Custom run with a YAML config override
Rscript ai_toxicogenomics_master_pipeline.R --config_file config.yaml

# Direct CLI overrides
Rscript ai_toxicogenomics_master_pipeline.R --n_samples 800 --bo_iterations 25 --output_dir results/

# Skip the inline unit tests
Rscript ai_toxicogenomics_master_pipeline.R --run_tests FALSE
```

Command-line values are overridden by the YAML file when both are given. `--n_samples` must be ≥ 30 and `--n_genes` ≥ 2000.

Sourcing the file interactively (`source("ai_toxicogenomics_master_pipeline.R")`) does **not** auto-execute the pipeline — each stage is a standalone function (`simulate_multiomics_cohort()`, `run_deseq2_lrt()`, `run_wgcna_analysis()`, etc.), so individual stages can be run, inspected, or unit-tested in isolation. Full execution requires calling `run_pipeline(opts)` explicitly or invoking the script via `Rscript`.

## Notes on the Synthetic Data

The RNA-seq cohort and gene-set panel in this repository are simulated (negative-binomial counts with injected dose-dependent modules) so the full pipeline is runnable end-to-end without a data-access agreement. The gene-set panel used for GSVA is illustrative; production use should substitute a curated panel (e.g., MSigDB Hallmark via `msigdbr`) matched to the exposure/tissue context.

Because the autoencoder, WGCNA and GSVA are unsupervised and fit on all samples (before the train/test split), the held-out R² is a slightly optimistic estimate of performance on a truly independent cohort; for real data, refit these stages inside each training fold.
