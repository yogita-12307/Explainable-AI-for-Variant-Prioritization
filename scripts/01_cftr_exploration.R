# --------------------------------------------------
# CFTR Variant Extraction and Cleaning from ClinVar
# Project: Rare Disease Variant Prioritization
# --------------------------------------------------

library(data.table)

# Set working directory containing ClinVar files
setwd("E:/Private/variant_summary.txt")
getwd()

# Load ClinVar variant summary dataset
clinVar <- fread("variant_summary.txt")

# Extract variants belonging to the CFTR gene
cftr <- clinVar[GeneSymbol == "CFTR"]

# Display total number of CFTR variants
cat(
  "Total CFTR variants identified:",
  nrow(cftr),
  "\n"
)

# Summarize clinical significance categories
clinical_summary <- as.data.frame(
  table(cftr$ClinicalSignificance)
)

print(clinical_summary)

# Retain only variants mapped to the GRCh38 assembly
cftr <- cftr[Assembly == "GRCh38"]
nrow(cftr)

# Keep only high-confidence pathogenicity classes
good_labels <- c(
  "Pathogenic",
  "Likely pathogenic",
  "Benign",
  "Likely benign"
)

cftr_clean <- cftr[
  ClinicalSignificance %in% good_labels
]

# Examine class distribution after filtering
table(cftr_clean$ClinicalSignificance)

# Review ClinVar evidence levels for retained variants
table(cftr_clean$ReviewStatus)

# Create results directory if it does not exist
dir.create("results", showWarnings = FALSE)

# Save variant count summary
writeLines(
  paste(
    "Total CFTR variants identified:",
    nrow(cftr)
  ),
  "results/cftr_variant_count.txt"
)

# Save clinical significance summary table
write.csv(
  clinical_summary,
  "results/clinical_significance_summary.csv",
  row.names = FALSE
)
nrow(cftr_clean)
nrow(cftr)
head(cftr$GeneSymbol)
# Save cleaned dataset for downstream machine learning
write.csv(
  cftr_clean,
  "results/cftr_clean_dataset.csv",
  row.names = FALSE
)
