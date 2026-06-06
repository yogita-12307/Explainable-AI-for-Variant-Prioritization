getwd()
setwd("E:/Private/variant_summary.txt")

library(data.table)

# Load ClinVar dataset
clinVar <- fread("variant_summary.txt")

# Extract CFTR variants
cftr <- clinVar[GeneSymbol == "CFTR"]

# Count variants
cat(
  "Total CFTR variants identified:",
  nrow(cftr),
  "\n"
)

# Clinical significance distribution
clinical_summary <- as.data.frame(
  table(cftr$ClinicalSignificance)
)

print(clinical_summary)

# Create result folder
dir.create("results")

# Save outputs
writeLines(
  paste(
    "Total CFTR variants identified:",
    nrow(cftr)
  ),
  "results/cftr_variant_count.txt"
)

write.csv(
  clinical_summary,
  "results/clinical_significance_summary.csv",
  row.names = FALSE
)
