# Explainable-AI-for-Variant-Prioritization
# CFTR Variant Prioritization

## Project Overview

This project aims to develop an interpretable machine learning framework for prioritizing disease-causing genetic variants using ClinVar knowledge integration and Whole Exome Sequencing (WES/WXS) data. The initial focus is on the CFTR gene as a proof-of-concept case study for rare disease variant prioritization.

## Current Progress

- Downloaded ClinVar variant_summary dataset
- Loaded ClinVar data in R
- Extracted CFTR variants
- Analyzed clinical significance distribution

## Dataset

Source:
ClinVar Variant Summary Dataset

Gene of Interest:
CFTR

Number of CFTR Variants:
12421

Assembly: GRCh38

## Clinical Labels Retained:

Pathogenic
Likely pathogenic
Benign
Likely benign
Repository Structure
scripts/
    cftr_preprocessing.R

results/
    cftr_variant_count.txt
    clinical_significance_summary.csv

## Next Steps
Create machine-learning-ready labels
Perform exploratory data analysis
Engineer predictive features
Train baseline machine learning models
Evaluate classification performance
Implement explainable AI techniques
Extend framework toward WES variant prioritization

## Long-Term Goal
Develop an explainable machine learning framework for prioritizing disease-causing variants from Whole Exome Sequencing data using ClinVar knowledge integration for rare disease diagnosis.
