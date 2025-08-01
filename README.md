# ASD-module
The code implements matrix factorization and machine learning to identify ASD-related APA gene modules from human brain snRNA-seq data.
# About
This study presents a systematic and integrative analytical pipeline for identifying autism spectrum disorder (ASD)-associated alternative polyadenylation (APA) gene modules and evaluating their predictive utility. 
Cell-type-specific RUD matrices are first randomly partitioned into training and test sets at a 7:3 ratio, with the partitioning process repeated 10 times to ensure robustness. Sparse Module Activity Factorization (SMAF) is then applied to the training data to decompose APA profiles into a gene-module matrix (U) and a module-cell activity matrix (W), enabling the identification of latent APA co-regulation patterns. A statistical filtering strategy is subsequently employed to detect modules that are significantly dysregulated in ASD and strongly associated with phenotypic outcomes. Core genes and distinctive APA features are extracted from each significant module. To enhance reliability, module stability is assessed using "cross-split reproducibility," retaining only those modules consistently identified across independent data partitions. Finally, ASD classification models are constructed using XGBoost (single-modality) and PSVM-2V (multi-modality fusion) based on the derived APA modules, allowing evaluation of their diagnostic potential.
By integrating unsupervised module discovery with supervised prediction, this pipeline provides an interpretable and biologically grounded framework for elucidating the role of APA in ASD pathogenesis.
# Requirements
You'll need to install the following packages in order to run the codes.
  * python==3.8
  * R==4.3
# Example codes
## preprocessing.R
The data were cleaned by checking the quality of nuclei and genes, and removing a few nuclei from different cell cycle stages.
## dataset_splitting.R
The expression data for each cell type were split into training and test sets.
## run_smaf.py
The SMAF algorithm was applied to decompose the RUD training set matrix for each cell type to identify APA gene modules.
## identify_modules.R
The core code identifies APA gene modules through differential analysis and correlation testing, and evaluates module robustness by combining gene selection strategies with cross-data-block recurrence rates.
## trained_model.R
Core code for model comparison.