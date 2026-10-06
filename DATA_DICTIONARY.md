#  Data Dictionary & Multi-Omics Codebook

This document describes the high-dimensional variables, simulated transcriptomic signatures, and machine learning features utilized within the Toxicogenomics & AI pipeline.

## 1. Experimental Metadata & Phenotype
| Variable Name | Data Type | Description | Unit / Coding |
| :--- | :--- | :--- | :--- |
| `Sample_ID` | String | Unique identifier for the simulated hepatocyte culture | `Hepatocyte_1` to `N` |
| `Batch` | Categorical | Technical sequencing batch (used for ComBat-seq correction); a gene-wise multiplicative batch effect is simulated for `Lab_B`/`Lab_C` | `Lab_A`, `Lab_B`, `Lab_C` |
| `Treatment_Group`| Categorical | Experimental hepatotoxin (APAP) exposure level | `Vehicle`, `Low_Dose`, `High_Dose` |

## 2. Transcriptomic Features (RNA-Seq)
| Variable Name | Data Type | Description | Unit / Coding |
| :--- | :--- | :--- | :--- |
| `ENSG000...` | Numeric | Ensembl Gene ID. Contains negative binomial read counts | Raw Integer Counts |
| `ME0` to `ME_N` | Numeric | WGCNA Module Eigengenes (first principal component of co-expressed gene clusters); `ME0` is the unassigned (grey) set | Continuous |
| `HALLMARK_*_LIKE` | Numeric | GSVA single-sample enrichment score for an illustrative hepatotoxicity gene set (oxidative stress, apoptosis, necrosis, xenobiotic metabolism) | Continuous, −1 to 1 |

## 3. Artificial Intelligence / Deep Learning Features
| Variable Name | Data Type | Description | Unit / Coding |
| :--- | :--- | :--- | :--- |
| `DL_Latent_Dim_X`| Numeric | Compressed 32-dimensional latent space features extracted via Keras Autoencoder | Scaled Continuous |

## 4. Toxicological Ground Truth (Target Variable)
| Variable Name | Data Type | Description | Unit / Coding |
| :--- | :--- | :--- | :--- |
| `Toxicity_Score` | Numeric | Simulated severity of hepatic necrosis and injury (model target, called `Toxicity` in the modeling table) | Scale `0` - `100` |

*Note: All datasets are synthetically generated for pipeline validation and computational demonstration.*
