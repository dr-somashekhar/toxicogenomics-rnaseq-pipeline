#  AI-Integrated Toxicogenomics & Multi-Omics Engine 
**Deep Learning (Keras/TensorFlow), WGCNA, and XGBoost for Predictive Toxicology**

##  Project Overview
This repository contains a high-dimensional, enterprise-grade programmatic architecture designed to model **Drug-Induced Liver Injury (DILI)**. Going far beyond traditional differential expression, this pipeline integrates Artificial Intelligence (AI) and Systems Biology to predict hepatotoxicity severity directly from raw RNA-seq transcriptomic signatures.

By bridging deep biological knowledge (Pharm.D.) with advanced machine learning, this pipeline demonstrates PhD-level capabilities in systems toxicology, biomarker discovery, and computational biology.

---

##  Advanced Computational Architecture

### 1. Batch Correction & Likelihood Ratio Testing (DESeq2 / ComBat-seq)
Raw sequencing data is intrinsically noisy. This pipeline utilizes `ComBat-seq` for robust removal of batch effects (adjusting for multi-center laboratory variance). Differential expression is modeled via `DESeq2` utilizing a **Likelihood Ratio Test (LRT)** to capture complex, dose-response transcriptomic perturbations across Vehicle, Low-Dose, and High-Dose exposure groups.

### 2. Weighted Gene Co-expression Network Analysis (WGCNA)
To move from individual genes to biological systems, the pipeline constructs a scale-free topological network. By computing a Topological Overlap Matrix (TOM), it clusters 20,000 genes into functional modules (e.g., Oxidative Stress, Apoptosis). The Module Eigengenes (MEs) are extracted to represent entire molecular pathways mathematically.

### 3. Deep Learning Autoencoder (Keras / TensorFlow)
High-dimensional omics data (20,000+ features) suffers from the "curse of dimensionality" when fed into predictive models. This script builds a Deep Neural Network Autoencoder to compress the transcriptomic data into a dense, 32-dimensional **Latent Space**. This performs non-linear feature extraction, capturing the underlying biological geometry of hepatotoxicity.

### 4. Predictive Machine Learning (XGBoost)
The pipeline merges the Systems Biology metrics (WGCNA Eigengenes) and the AI metrics (Deep Learning Latent Features) into a unified dataset. An **Extreme Gradient Boosting (XGBoost)** regression model is trained to predict the exact quantitative severity of liver injury. The model utilizes cross-validation, hyperparameter tuning, and early stopping, outputting feature importance metrics to identify the strongest predictors of necrosis.

---

##  Technical Stack
* **Language:** R (Version 4.2+)
* **Deep Learning:** `keras`, `tensorflow`
* **Machine Learning:** `xgboost`, `caret`
* **Systems Biology:** `WGCNA`, `DESeq2`, `sva`, `GSVA`
* **Infrastructure:** `doParallel` (Multi-core tensor processing)

##  Execution
To run this multi-omics AI pipeline:
```R
source("ai_toxicogenomics_master_pipeline.R")
