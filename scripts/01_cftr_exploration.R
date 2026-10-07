# ============================================================
# CFTR PATHOGENICITY PROJECT
# MODEL A vs MODEL B — CLEAN REBUILD
# ============================================================

# ------------------------------------------------------------
# 0. Packages
# ------------------------------------------------------------

library(data.table)
library(caret)
library(randomForest)



set.seed(123)

# ------------------------------------------------------------
# 1. Load / prepare original CFTR dataset
# ------------------------------------------------------------

# Start from the original ClinVar dataset

setwd("E:/Private")
clinVar <- fread("E:/Private/variant_summary.txt/variant_summary.txt")
cftr <- clinVar[GeneSymbol == "CFTR"]

# Keep only GRCh38 variants and the four target classifications
cftr_clean <- cftr[
  Assembly == "GRCh38" &
    ClinicalSignificance %in% c(
      "Pathogenic",
      "Likely pathogenic",
      "Benign",
      "Likely benign"
    )
]

# Create binary outcome:
# 1 = pathogenic / likely pathogenic
# 0 = benign / likely benign
cftr_clean[, labels := ifelse(
  ClinicalSignificance %in% c(
    "Pathogenic",
    "Likely pathogenic"
  ),
  1,
  0
)]

# Make outcome a factor for caret
cftr_clean[, labels_factor := factor(
  labels,
  levels = c(0, 1),
  labels = c("Benign", "Pathogenic")
)]

cat("\n================ DATASET =================\n")
cat("Number of variants:", nrow(cftr_clean), "\n")
print(table(cftr_clean$labels_factor))


# ------------------------------------------------------------
# 2. Create ONE train/test split
# ------------------------------------------------------------
# IMPORTANT:
# Both Model A and Model B use exactly the same variants
# in training and testing.

set.seed(123)

train_index <- createDataPartition(
  cftr_clean$labels_factor,
  p = 0.80,
  list = FALSE
)

train_data <- cftr_clean[train_index]
test_data  <- cftr_clean[-train_index]

cat("\n================ TRAIN / TEST =================\n")
cat("Training variants:", nrow(train_data), "\n")
cat("Testing variants:", nrow(test_data), "\n")

cat("\nTraining class distribution:\n")
print(table(train_data$labels_factor))

cat("\nTesting class distribution:\n")
print(table(test_data$labels_factor))


# ============================================================
# MODEL A
# CURATION-DEPENDENT BASELINE
# ============================================================

# Model A deliberately includes ClinVar curation metadata.
#
# These variables describe aspects of the ClinVar evidence/
# curation process and therefore make this a curation-dependent
# baseline rather than a purely variant-level model.

# ------------------------------------------------------------
# 3. Prepare Model A variables
# ------------------------------------------------------------

modelA_vars <- c(
  "labels_factor",
  "Type",
  "OriginSimple",
  "ReviewStatus",
  "NumberSubmitters"
)

modelA_train <- train_data[, ..modelA_vars]
modelA_test  <- test_data[, ..modelA_vars]

# Convert categorical variables to factors
modelA_train[, Type := as.factor(Type)]
modelA_test[, Type := factor(Type, levels = levels(modelA_train$Type))]

modelA_train[, OriginSimple := as.factor(OriginSimple)]
modelA_test[, OriginSimple := factor(
  OriginSimple,
  levels = levels(modelA_train$OriginSimple)
)]

modelA_train[, ReviewStatus := as.factor(ReviewStatus)]
modelA_test[, ReviewStatus := factor(
  ReviewStatus,
  levels = levels(modelA_train$ReviewStatus)
)]

# Replace missing categorical values with explicit level
for (v in c("Type", "OriginSimple", "ReviewStatus")) {
  
  if ("Missing" %in% levels(modelA_train[[v]])) {
    next
  }
  
  levels(modelA_train[[v]]) <- c(
    levels(modelA_train[[v]]),
    "Missing"
  )
  
  modelA_train[is.na(get(v)), (v) := "Missing"]
  
  modelA_test[is.na(get(v)), (v) := "Missing"]
}

# Make sure NumberSubmitters is numeric
modelA_train[, NumberSubmitters := as.numeric(NumberSubmitters)]
modelA_test[, NumberSubmitters := as.numeric(NumberSubmitters)]

# Replace missing numeric values with training median
median_submitters <- median(
  modelA_train$NumberSubmitters,
  na.rm = TRUE
)

modelA_train[
  is.na(NumberSubmitters),
  NumberSubmitters := median_submitters
]

modelA_test[
  is.na(NumberSubmitters),
  NumberSubmitters := median_submitters
]


# ------------------------------------------------------------
# 4. Train Model A
# ------------------------------------------------------------

set.seed(123)

modelA_rf <- randomForest(
  labels_factor ~ Type +
    OriginSimple +
    ReviewStatus +
    NumberSubmitters,
  data = modelA_train,
  ntree = 500,
  importance = TRUE
)

cat("\n================ MODEL A =================\n")
print(modelA_rf)

cat("\nModel A variable importance:\n")
print(importance(modelA_rf))


# ------------------------------------------------------------
# 5. Evaluate Model A
# ------------------------------------------------------------

modelA_pred <- predict(
  modelA_rf,
  newdata = modelA_test
)

cat("\nModel A confusion matrix:\n")

cm_A <- confusionMatrix(
  modelA_pred,
  modelA_test$labels_factor,
  positive = "Pathogenic"
)

print(cm_A)

# ============================================================
# FIX MODEL A TEST DATA
# ============================================================

# Rebuild the test data directly from the original train/test
# split so that factor types match the model exactly.

modelA_test <- test_data[, ..modelA_vars]

# Match factor levels exactly to the training data
modelA_test[, Type := factor(
  Type,
  levels = levels(modelA_train$Type)
)]

modelA_test[, OriginSimple := factor(
  OriginSimple,
  levels = levels(modelA_train$OriginSimple)
)]

modelA_test[, ReviewStatus := factor(
  ReviewStatus,
  levels = levels(modelA_train$ReviewStatus)
)]

# NumberSubmitters must be numeric
modelA_test[, NumberSubmitters := as.numeric(NumberSubmitters)]

# Replace missing values using the training median
modelA_test[
  is.na(NumberSubmitters),
  NumberSubmitters := median_submitters
]

# Predict Model A
modelA_pred <- predict(
  modelA_rf,
  newdata = modelA_test
)

# Evaluate
cat("\n================ MODEL A TEST RESULTS =================\n")

cm_A <- confusionMatrix(
  modelA_pred,
  modelA_test$labels_factor,
  positive = "Pathogenic"
)

print(cm_A)

cat("\nKey metrics:\n")
cat("Accuracy:",
    round(cm_A$overall["Accuracy"], 4), "\n")

cat("Kappa:",
    round(cm_A$overall["Kappa"], 4), "\n")

cat("Pathogenic sensitivity:",
    round(cm_A$byClass["Sensitivity"], 4), "\n")

cat("Pathogenic specificity:",
    round(cm_A$byClass["Specificity"], 4), "\n")

cat("Balanced accuracy:",
    round(cm_A$byClass["Balanced Accuracy"], 4), "\n")


# ============================================================
# MODEL B
# REDUCED-CIRCULARITY VARIANT-LEVEL MODEL
# ============================================================

# IMPORTANT:
# We DO NOT filter out "other" variants.
#
# The previous parser caused severe selection bias because
# many benign variants were classified as "other".
#
# Instead, every one of the original 2,800 variants remains
# in the dataset.


# ------------------------------------------------------------
# 6. Create broader consequence categories
# ------------------------------------------------------------

# Work from the HGVS/name field already used in your project.
#
# The rules are deliberately broad.
# Unrecognized variants remain "other" rather than being removed.

cftr_clean[, consequence := "other"]

# Frameshift
cftr_clean[
  grepl("fs", Name, ignore.case = TRUE),
  consequence := "frameshift"
]

# Nonsense / stop-gain
cftr_clean[
  grepl("Ter|\\*", Name, ignore.case = TRUE),
  consequence := "nonsense"
]

# Splice-site variants
cftr_clean[
  grepl(
    "\\+[1-9][0-9]*[A-Z]>|-[1-9][0-9]*[A-Z]>",
    Name,
    ignore.case = TRUE
  ),
  consequence := "splice_site"
]

# Deletions
cftr_clean[
  grepl("del", Name, ignore.case = TRUE),
  consequence := "deletion"
]

# Duplications
cftr_clean[
  grepl("dup", Name, ignore.case = TRUE),
  consequence := "duplication"
]

# Insertions
cftr_clean[
  grepl("ins", Name, ignore.case = TRUE),
  consequence := "insertion"
]

# Missense
cftr_clean[
  grepl(
    "[A-Z][a-z]{2}[0-9]+[A-Z][a-z]{2}",
    Name
  ),
  consequence := "missense"
]

# Synonymous variants
cftr_clean[
  grepl("=", Name),
  consequence := "synonymous"
]

# ------------------------------------------------------------
# 7. Inspect consequence distribution
# ------------------------------------------------------------

cat("\n================ CONSEQUENCE DISTRIBUTION =================\n")

print(table(cftr_clean$consequence))

cat("\nConsequence by clinical class:\n")

print(
  table(
    cftr_clean$labels_factor,
    cftr_clean$consequence
  )
)

cat("\nProportions within each clinical class:\n")

print(
  round(
    prop.table(
      table(
        cftr_clean$labels_factor,
        cftr_clean$consequence
      ),
      margin = 1
    ),
    3
  )
)


# ============================================================
# MODEL B TRAINING
# ============================================================

# ------------------------------------------------------------
# 8. Recreate train/test data after adding consequence
# ------------------------------------------------------------

# The original train_index is retained, so Model A and Model B
# use exactly the same train/test variants.

train_data_B <- cftr_clean[train_index]
test_data_B  <- cftr_clean[-train_index]

train_data_B[, consequence := as.factor(consequence)]

test_data_B[, consequence := factor(
  consequence,
  levels = levels(train_data_B$consequence)
)]


# ------------------------------------------------------------
# 9. Add basic variant-level information
# ------------------------------------------------------------

# Type is a variant-level attribute rather than a ClinVar
# curation-confidence variable.

train_data_B[, Type := as.factor(Type)]

test_data_B[, Type := factor(
  Type,
  levels = levels(train_data_B$Type)
)]

# Remove unused factor levels
train_data_B[, consequence := droplevels(consequence)]
test_data_B[, consequence := factor(
  consequence,
  levels = levels(train_data_B$consequence)
)]


# ------------------------------------------------------------
# 10. Train Model B
# ------------------------------------------------------------

set.seed(123)

modelB_rf <- randomForest(
  labels_factor ~ Type + consequence,
  data = train_data_B,
  ntree = 500,
  importance = TRUE
)

cat("\n================ MODEL B =================\n")
print(modelB_rf)

cat("\nModel B variable importance:\n")
print(importance(modelB_rf))


# ------------------------------------------------------------
# 11. Evaluate Model B
# ------------------------------------------------------------

modelB_pred <- predict(
  modelB_rf,
  newdata = test_data_B
)

cat("\nModel B confusion matrix:\n")

cm_B <- confusionMatrix(
  modelB_pred,
  test_data_B$labels_factor,
  positive = "Pathogenic"
)

print(cm_B)


# ============================================================
# 12. DIRECT MODEL COMPARISON
# ============================================================

cat("\n============================================================\n")
cat("                MODEL A vs MODEL B\n")
cat("============================================================\n")

cat("\nModel A:\n")
cat("Accuracy:", round(cm_A$overall["Accuracy"], 4), "\n")
cat(
  "Kappa:",
  round(cm_A$overall["Kappa"], 4),
  "\n"
)
cat(
  "Pathogenic sensitivity:",
  round(cm_A$byClass["Sensitivity"], 4),
  "\n"
)
cat(
  "Pathogenic specificity:",
  round(cm_A$byClass["Specificity"], 4),
  "\n"
)
cat(
  "Balanced accuracy:",
  round(cm_A$byClass["Balanced Accuracy"], 4),
  "\n"
)

cat("\nModel B:\n")
cat("Accuracy:", round(cm_B$overall["Accuracy"], 4), "\n")
cat(
  "Kappa:",
  round(cm_B$overall["Kappa"], 4),
  "\n"
)
cat(
  "Pathogenic sensitivity:",
  round(cm_B$byClass["Sensitivity"], 4),
  "\n"
)
cat(
  "Pathogenic specificity:",
  round(cm_B$byClass["Specificity"], 4),
  "\n"
)
cat(
  "Balanced accuracy:",
  round(cm_B$byClass["Balanced Accuracy"], 4),
  "\n"
)


# ============================================================
# 13. Save important objects
# ============================================================

results <- list(
  dataset = cftr_clean,
  train_index = train_index,
  modelA = modelA_rf,
  modelB = modelB_rf,
  modelA_confusion = cm_A,
  modelB_confusion = cm_B
)

saveRDS(
  results,
  file = "CFTR_model_results.rds"
)

cat("\n============================================================\n")
cat("Analysis complete.\n")
cat("Results saved to: CFTR_model_results.rds\n")
cat("============================================================\n")
