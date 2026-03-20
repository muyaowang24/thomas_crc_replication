# Replication of Thomas et al.on cross-cohort metagenomic signatures of colorectal cancer

This repository contains a public-data replication of the analytical framework from Thomas et al., *Nature Medicine* (2019), focusing on cross-cohort metagenomic signatures of colorectal cancer.

## Project goal

The original paper asked a simple but important question:

Can gut microbial signals linked to colorectal cancer be reproduced across independent cohorts, rather than appearing only in a single dataset?

This replication focuses on that core idea using publicly accessible processed metagenomic profiles from `curatedMetagenomicData`.

**This replication project is also intended as a learning-oriented research exercise. In addition to reproducing the main analytical framework of the original paper, it helps build practical understanding of cross-cohort microbiome analysis, meta-analysis, external validation, and reproducible research workflow design**.

## What this project does

This project reproduces three main parts of the original analytical logic:

1. Load and harmonize public colorectal cancer stool metagenomic cohorts
2. Identify species that repeatedly differ between CRC cases and controls across studies
3. Test whether species-level microbial signals generalize across studies using leave-one-dataset-out classification

## Data source

Data were accessed through:

- `curatedMetagenomicData`
- `curatedMetagenomicAnalyses`

These packages provide standardized microbiome abundance profiles and metadata for public colorectal cancer cohorts used in the replication.

This repository does not redistribute raw sequencing files.

## Main results

Using 11 public CRC cohorts and 1395 total stool samples, this replication found:

- Several species repeatedly enriched in CRC across studies, including  
  `Parvimonas micra`, `Gemella morbillorum`, `Peptostreptococcus stomatis`, `Dialister pneumosintes`, and `Fusobacterium nucleatum`
- Random-effects meta-analysis showed positive pooled effects for these recurrent taxa
- Leave-one-dataset-out classification showed moderate to strong external performance across most held-out cohorts

These results support the main idea of the original paper: some CRC-associated microbial signals are reproducible across independent studies.

## Project structure

```text
.
├── README.md
├── .gitignore
├── data_raw/
├── data_processed/
├── output/
├── results/
├── R/
└── report/
