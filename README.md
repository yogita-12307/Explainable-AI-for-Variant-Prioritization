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

## Clinical Labels Retained

- Pathogenic
- Likely pathogenic
- Benign
- Likely benign

## Feature Engineering

Selected machine learning features:

- Type
- OriginSimple
- ReviewStatus
- NumberSubmitters

Generated binary labels:

- 1 = Pathogenic + Likely pathogenic
- 0 = Benign + Likely benign

Dataset Summary:

- Total CFTR variants: 12,421
- GRCh38 CFTR variants: 6,094
- High-confidence variants: 2,800
- Pathogenic class: 1,212
- Benign class: 1,588

## Repository Structure

scripts/
    cftr_preprocessing.R

results/
    cftr_variant_count.txt
    cftr_clean_variant_count.txt
    clinical_significance_summary.csv
    cftr_ml_dataset.csv
    Random Forest Performance

Accuracy: 84.79%
Kappa: 0.6849
Sensitivity: 92.11%
Specificity: 75.21%
OOB Error Rate: 14.23%

Top Predictive Features:
1. Type
2. ReviewStatus
3. NumberSubmitters
4. OriginSimple

## Next Steps

- Split dataset into training and testing sets
- Train Random Forest classifier
- Evaluate model performance
- Perform feature importance analysis
- Implement explainable AI methods
- Integrate WES-derived variants for prioritization

## Long-Term Goal

Develop an explainable machine learning framework for prioritizing disease-causing variants from Whole Exome Sequencing (WES) data through integration of ClinVar clinical knowledge, variant annotation, and interpretable machine learning models.
Implement explainable AI techniques
Extend framework toward WES variant prioritization

Develop an explainable machine learning framework for prioritizing disease-causing variants from Whole Exome Sequencing data using ClinVar knowledge integration for rare disease diagnosis.
