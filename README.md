#  Toxicogenomics RNA-Seq Pipeline: Drug-Induced Liver Injury (DILI)
**Advanced Bioinformatics Workflow for Differential Gene Expression and Pathway Enrichment**

##  Project Overview
This repository contains an enterprise-grade computational biology pipeline designed to process and analyze high-dimensional transcriptomic data. Specifically, it models **Drug-Induced Liver Injury (DILI)** by evaluating the molecular response of human hepatocytes exposed to a known hepatotoxin (e.g., Acetaminophen/APAP).

Bridging clinical pharmacotherapy with molecular biology, this pipeline transforms raw RNA sequencing counts into actionable toxicological insights, identifying transcriptomic biomarkers of cellular stress, apoptosis, and metabolic dysfunction.

---

## Computational & Molecular Workflow

### 1. High-Dimensional Transcriptomic Simulation
To demonstrate robust data wrangling capabilities, the pipeline simulates a massive RNA-seq count matrix containing **20,000 genes across 12 biological replicates** (6 Controls vs. 6 Exposed). The baseline read counts are modeled using a negative binomial distribution, with true biological signals injected to simulate transcriptomic alterations (upregulation of oxidative stress and downregulation of basal metabolic pathways).

### 2. Differential Gene Expression (DESeq2)
The core statistical engine utilizes the Bioconductor package `DESeq2`. 
* **Pre-filtering:** Removes low-count genes to optimize memory and statistical power.
* **Dispersion Estimation & Modeling:** Fits a generalized linear model (GLM) using negative binomial distribution assumptions.
* **FDR Control & Log-Fold Shrinkage:** Applies the Benjamini-Hochberg procedure for strict False Discovery Rate (FDR < 0.05) control, and utilizes `apeglm` to shrink extreme Log2 Fold Changes (LFC) associated with low-expression genes, reducing false positives.

### 3. Dimensionality Reduction & Quality Control
Before executing statistical tests, the pipeline verifies data integrity and batch effects:
* Applies **Variance Stabilizing Transformation (VST)** to homoskedasticize the count data.
* Executes **Principal Component Analysis (PCA)** to visualize sample clustering, ensuring separation is driven by hepatotoxin exposure rather than technical artifacts.

### 4. Functional Pathway Enrichment Analysis
Identifying a list of differentially expressed genes (DEGs) is only the first step. The pipeline utilizes `clusterProfiler` to translate gene lists into biological meaning:
* **Over-Representation Analysis (ORA):** Maps significantly perturbed genes to the **Gene Ontology (GO)** database.
* **Mechanism Identification:** Programmatically identifies whether the hepatotoxin triggers specific biological processes (BP) such as cytochrome P450 (CYP) induction, reactive oxygen species (ROS) metabolism, or caspase-mediated apoptosis.

---

##  Translational Relevance
This project directly reflects the methodologies utilized in modern translational toxicology. By moving beyond traditional serum biomarkers (like ALT/AST) and evaluating toxicity at the transcriptomic level, researchers can identify the exact molecular mechanisms of adverse drug reactions long before macroscopic cellular necrosis occurs. This computational approach is foundational for predictive toxicology and the development of safer pharmaceutical compounds.

---

## Technical Stack
* **Language:** R (Version 4.2+)
* **Core Bioinformatics:** `Bioconductor`
* **Differential Expression:** `DESeq2`
* **Data Visualization:** `EnhancedVolcano`, `pheatmap`, `ggplot2`
* **Gene Set Enrichment:** `clusterProfiler`, `org.Hs.eg.db`
* **Data Wrangling:** `tidyverse`

##  Execution
To run the full simulation and analysis pipeline locally:
```R
# Ensure BiocManager is installed, then run:
source("dili_transcriptomics_pipeline.R")
