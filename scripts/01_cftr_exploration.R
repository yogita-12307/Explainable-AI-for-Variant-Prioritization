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

# Creates binary labels
cftr_clean$labels <- ifelse(
  cftr_clean$ClinicalSignificance %in% c("Pathogenic", "Likely pathogenic" ),
  1,0
)

table(cftr_clean$labels)

# Feature Exploration

# Label distribution
table(cftr_clean$labels)

# Variant type distribution
table(cftr_clean$Type)

# Variant origin distribution
table(cftr_clean$OriginSimple)

# Number of submitters summary
summary(cftr_clean$NumberSubmitters)

# Review status distribution
table(cftr_clean$ReviewStatus)

# Create Machine Learning Dataset
ml_data <- cftr_clean[, .(
  Type,
  OriginSimple,
  ReviewStatus,
  NumberSubmitters,
  labels
)]

# Inspect dataset
str(ml_data)
summary(ml_data)


# Convert Variables to Factors
ml_data$Type <- as.factor(ml_data$Type)
ml_data$OriginSimple <- as.factor(ml_data$OriginSimple)
ml_data$ReviewStatus <- as.factor(ml_data$ReviewStatus)
ml_data$labels <- as.factor(ml_data$labels)

# Verify structure
str(ml_data)

# Save ML Dataset
write.csv(
  ml_data,
  "results/cftr_ml_dataset.csv",
  row.names = FALSE
)

# Split data into training and testing sets

library(caret)

set.seed(123)

trainIndex <- createDataPartition(
  ml_data$labels,
  p = 0.8,
  list = FALSE
)

train_data <- ml_data[trainIndex, ]
test_data <- ml_data[-trainIndex, ]

dim(train_data)
dim(test_data)

write.csv(
  train_data,
  "results/train_data.csv",
  row.names = FALSE
)

write.csv(
  test_data,
  "results/test_data.csv",
  row.names = FALSE
)

# Train Random Forest classifier

library(randomForest)

rf_model <- randomForest(
  labels ~ .,
  data = train_data,
  ntree = 500,
  importance = TRUE
)

print(rf_model)

# Generate predictions

predictions <- predict(
  rf_model,
  test_data
)

# Evaluate model performance

rf_results <- confusionMatrix(
  predictions,
  test_data$labels
)

print(rf_results)

# Feature importance analysis

importance(rf_model)

varImpPlot(rf_model)

# Save model

saveRDS(
  rf_model,
  "results/cftr_random_forest_model.rds"
)

# Save evaluation results

capture.output(
  rf_results,
  file = "results/random_forest_results.txt"
)
