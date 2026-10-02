library(data.table)
library(caret)
library(randomForest)

# 1) LOAD AND CLEAN CLINVAR DATA  (unchanged from your original)

setwd("E:/Private")  # folder containing variant_summary.txt

clinVar <- fread("variant_summary.txt")

cftr <- clinVar[GeneSymbol == "CFTR"]

cat("Total CFTR variants identified:", nrow(cftr), "\n")

clinical_summary <- as.data.frame(table(cftr$ClinicalSignificance))
print(clinical_summary)

cftr <- cftr[Assembly == "GRCh38"]

good_labels <- c("Pathogenic", "Likely pathogenic", "Benign", "Likely benign")
cftr_clean <- cftr[ClinicalSignificance %in% good_labels]

table(cftr_clean$ClinicalSignificance)
table(cftr_clean$ReviewStatus)

dir.create("results", showWarnings = FALSE)

writeLines(
  paste("Total CFTR variants identified:", nrow(cftr)),
  "results/cftr_variant_count.txt"
)
write.csv(clinical_summary, "results/clinical_significance_summary.csv", row.names = FALSE)

write.csv(cftr_clean, "results/cftr_clean_dataset.csv", row.names = FALSE)

# Binary labels
cftr_clean$labels <- ifelse(
  cftr_clean$ClinicalSignificance %in% c("Pathogenic", "Likely pathogenic"),
  1, 0
)
table(cftr_clean$labels)

# 2) BASELINE MODEL — original 4 features
#    (includes ClinVar curation-confidence metadata on purpose,
#    this is the model we are checking FOR circularity)

ml_data_baseline <- cftr_clean[, .(
  Type, OriginSimple, ReviewStatus, NumberSubmitters, labels
)]

ml_data_baseline$Type         <- as.factor(ml_data_baseline$Type)
ml_data_baseline$OriginSimple <- as.factor(ml_data_baseline$OriginSimple)
ml_data_baseline$ReviewStatus <- as.factor(ml_data_baseline$ReviewStatus)
ml_data_baseline$labels       <- as.factor(ml_data_baseline$labels)

write.csv(ml_data_baseline, "results/cftr_ml_dataset_baseline.csv", row.names = FALSE)

set.seed(123)
trainIndex_base <- createDataPartition(ml_data_baseline$labels, p = 0.8, list = FALSE)
train_base <- ml_data_baseline[trainIndex_base, ]
test_base  <- ml_data_baseline[-trainIndex_base, ]

rf_model_baseline <- randomForest(
  labels ~ ., data = train_base, ntree = 500, importance = TRUE
)
print(rf_model_baseline)

pred_base <- predict(rf_model_baseline, test_base)
results_baseline <- confusionMatrix(pred_base, test_base$labels)
print(results_baseline)

saveRDS(rf_model_baseline, "results/cftr_rf_model_baseline.rds")
capture.output(results_baseline, file = "results/random_forest_results_baseline.txt")

# 3) ANNOTATE VARIANTS — gnomAD allele frequency + CADD score
#    via myvariant.info (requires internet access)

# install.packages("BiocManager")
# BiocManager::install("myvariant")
library(myvariant)

# Build HGVS-style genomic query IDs from ClinVar's VCF-style columns.
# NOTE: adjust column names below if your variant_summary.txt uses
# different headers (ClinVar's standard columns are shown here:
# Chromosome, Start, ReferenceAlleleVCF, AlternateAlleleVCF)
cftr_clean$query_id <- paste0(
  "chr", cftr_clean$Chromosome, ":g.",
  cftr_clean$Start, cftr_clean$ReferenceAlleleVCF, ">", cftr_clean$AlternateAlleleVCF
)

# Query in batches to avoid timeouts/rate limits
annotate_batch <- function(ids, batch_size = 200) {
  results <- list()
  for (i in seq(1, length(ids), by = batch_size)) {
    batch <- ids[i:min(i + batch_size - 1, length(ids))]
    res <- tryCatch(
      getVariants(batch, fields = c("gnomad_genome.af.af", "cadd.phred")),
      error = function(e) {
        cat("Batch", i, "failed:", conditionMessage(e), "\n")
        NULL
      }
    )
    if (!is.null(res)) results[[length(results) + 1]] <- res
  }
  if (length(results) == 0) return(NULL)
  rbindlist(results, fill = TRUE)
}

annotations <- annotate_batch(cftr_clean$query_id)

if (!is.null(annotations)) {
  cftr_clean$gnomad_af   <- annotations$gnomad_genome.af.af[match(cftr_clean$query_id, annotations$query)]
  cftr_clean$cadd_phred  <- annotations$cadd.phred[match(cftr_clean$query_id, annotations$query)]
} else {
  cat("Annotation failed — gnomad_af/cadd_phred will be NA. See fallback note below.\n")
  cftr_clean$gnomad_af  <- NA
  cftr_clean$cadd_phred <- NA
}

# Variants with no population frequency record are effectively absent
# from gnomAD (i.e. not seen in large population cohorts) — treat NA as 0
cftr_clean$gnomad_af[is.na(cftr_clean$gnomad_af)] <- 0

write.csv(cftr_clean, "results/cftr_annotated_dataset.csv", row.names = FALSE)

cat("Variants missing CADD score:", sum(is.na(cftr_clean$cadd_phred)), "of", nrow(cftr_clean), "\n")

# 4) BIOLOGICAL MODEL — no ClinVar curation-confidence features
#    This is the direct test of the original project goal:
#    "identify pathogenicity from the variant's own features,
#    not because ClinVar already flagged it with high confidence"

cftr_bio <- cftr_clean[!is.na(cadd_phred)]  # drop rows with no CADD score

ml_data_bio <- cftr_bio[, .(
  Type, OriginSimple, gnomad_af, cadd_phred, labels
)]

ml_data_bio$Type         <- as.factor(ml_data_bio$Type)
ml_data_bio$OriginSimple <- as.factor(ml_data_bio$OriginSimple)
ml_data_bio$labels       <- as.factor(ml_data_bio$labels)

write.csv(ml_data_bio, "results/cftr_ml_dataset_biological.csv", row.names = FALSE)

set.seed(123)
trainIndex_bio <- createDataPartition(ml_data_bio$labels, p = 0.8, list = FALSE)
train_bio <- ml_data_bio[trainIndex_bio, ]
test_bio  <- ml_data_bio[-trainIndex_bio, ]

rf_model_bio <- randomForest(
  labels ~ ., data = train_bio, ntree = 500, importance = TRUE
)
print(rf_model_bio)

pred_bio <- predict(rf_model_bio, test_bio)
results_bio <- confusionMatrix(pred_bio, test_bio$labels)
print(results_bio)

saveRDS(rf_model_bio, "results/cftr_rf_model_biological.rds")
capture.output(results_bio, file = "results/random_forest_results_biological.txt")

# 5) COMPARE BASELINE vs BIOLOGICAL MODEL

comparison <- data.frame(
  Model = c("Baseline (incl. ReviewStatus/NumberSubmitters)",
            "Biological (Type/Origin/gnomAD/CADD)"),
  Accuracy    = c(results_baseline$overall["Accuracy"], results_bio$overall["Accuracy"]),
  Sensitivity = c(results_baseline$byClass["Sensitivity"], results_bio$byClass["Sensitivity"]),
  Specificity = c(results_baseline$byClass["Specificity"], results_bio$byClass["Specificity"])
)
print(comparison)
write.csv(comparison, "results/model_comparison.csv", row.names = FALSE)

# 6) SHAP EXPLAINABILITY — on the biological model
#    (explains individual predictions, not just global importance)

# install.packages("fastshap")
library(fastshap)
library(ggplot2)

pred_wrapper <- function(object, newdata) {
  predict(object, newdata, type = "prob")[, "1"]
}

shap_values <- explain(
  rf_model_bio,
  X = train_bio[, -"labels"],
  pred_wrapper = pred_wrapper,
  nsim = 50
)

shap_plot <- autoplot(shap_values, type = "importance")
ggsave("results/shap_importance_biological_model.png", shap_plot, width = 7, height = 5)

saveRDS(shap_values, "results/shap_values_biological.rds")

cat("\nDone. Key outputs in results/:\n",
    "- model_comparison.csv          (baseline vs biological accuracy)\n",
    "- shap_importance_biological_model.png\n",
    "- random_forest_results_baseline.txt / _biological.txt\n")
