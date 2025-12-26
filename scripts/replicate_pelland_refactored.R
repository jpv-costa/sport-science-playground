# scripts/replicate_pelland_refactored.R
# Refactored Replication of Pelland et al. (2024) Meta-Regression
#
# Uses OOP/SOLID architecture with R6 classes:
# - PellandDataLoader: Data loading and preprocessing
# - CorrelationImputer: Pre-post correlation imputation
# - MetaRegressionModel: Multi-level meta-regression fitting
#
# Original OSF: https://osf.io/7knsj/

# ==============================================================================
# Setup
# ==============================================================================

# Ensure we're in project root for box module resolution
if (file.exists("deadlift-study.Rproj")) {
  # Already in project root
} else if (file.exists("../deadlift-study.Rproj")) {
  setwd("..")
}

# Set box path explicitly
options(box.path = getwd())

box::use(
  R/loaders/pelland_data_loader[PellandDataLoader],
  R/models/meta_regression[MetaRegressionModel, MetaRegressionResult]
)

cat("=== Pelland et al. (2024) Meta-Regression Replication ===\n")
cat("Using refactored OOP/SOLID architecture\n\n")

# ==============================================================================
# Configuration
# ==============================================================================

DATA_PATH <- "data/external/pelland_meta_regression/Data__Code__and_Estimation_Materials/V2.PTF.Data.xlsx"
OUTPUT_PATH <- "data/processed/pelland_replication_results_refactored.rds"

# ==============================================================================
# Data Loading
# ==============================================================================

cat("=== Loading and Preprocessing Data ===\n")

loader <- PellandDataLoader$new(DATA_PATH)
data <- loader$load()
summary_stats <- loader$summarize(data)

cat("Data Summary:\n")
cat("  Total effect sizes:", summary_stats$total_effects, "\n")
cat("  Strength effects:", summary_stats$strength_effects, "\n")
cat("  Hypertrophy effects:", summary_stats$hypertrophy_effects, "\n")
cat("  Strength studies:", summary_stats$strength_studies, "\n")
cat("  Hypertrophy studies:", summary_stats$hypertrophy_studies, "\n")
cat("  Imputed r (Strength):", round(summary_stats$imputed_correlation_strength, 4), "\n")
cat("  Imputed r (Hypertrophy):", round(summary_stats$imputed_correlation_hypertrophy, 4), "\n\n")

# ==============================================================================
# Model Fitting
# ==============================================================================

cat("=== Fitting Meta-Regression Models ===\n\n")

model <- MetaRegressionModel$new(
  moderators = c("avg.rir", "load.set", "set.rep.equated", "weeks", "train.status"),
  random_effects = list(~1|study, ~1|group, ~1|obs),
  method = "REML",
  test = "t"
)

# Fit both outcomes
results <- model$fit_both_outcomes(data$strength, data$hypertrophy)

# ==============================================================================
# Results Extraction
# ==============================================================================

cat("=== STRENGTH RESULTS ===\n")
strength_summary <- results$strength$summarize()
cat("N effect sizes:", strength_summary$n_effects, "\n")
cat("N studies:", strength_summary$n_studies, "\n")
cat("RIR coefficient:", round(strength_summary$rir_effect, 4), "\n")
cat("RIR p-value:", round(strength_summary$rir_p_value, 4), "\n")
cat("RIR significant:", strength_summary$rir_significant, "\n\n")

cat("Coefficients Table:\n")
print(results$strength$get_coefficients_table())

cat("\n\n=== HYPERTROPHY RESULTS ===\n")
hypertrophy_summary <- results$hypertrophy$summarize()
cat("N effect sizes:", hypertrophy_summary$n_effects, "\n")
cat("N studies:", hypertrophy_summary$n_studies, "\n")
cat("RIR coefficient:", round(hypertrophy_summary$rir_effect, 4), "\n")
cat("RIR p-value:", round(hypertrophy_summary$rir_p_value, 4), "\n")
cat("RIR significant:", hypertrophy_summary$rir_significant, "\n\n")

cat("Coefficients Table:\n")
print(results$hypertrophy$get_coefficients_table())

# ==============================================================================
# Paper Comparison
# ==============================================================================

cat("\n\n=== PAPER COMPARISON ===\n")
cat("STRENGTH:\n")
cat("  Paper claim: CI contains null (not significant)\n")
cat("  Our result: p =", round(strength_summary$rir_p_value, 3), "\n")
cat("  Match:", ifelse(!strength_summary$rir_significant, "YES", "NO"), "\n")

cat("\nHYPERTROPHY:\n")
cat("  Paper claim: Negative slope, CI excludes null\n")
cat("  Our result: b =", round(hypertrophy_summary$rir_effect, 4),
    ", p =", round(hypertrophy_summary$rir_p_value, 4), "\n")
cat("  Match:", ifelse(hypertrophy_summary$rir_significant &&
                        hypertrophy_summary$rir_effect < 0, "YES", "NO"), "\n")

# ==============================================================================
# Save Results
# ==============================================================================

cat("\n\n=== Saving Results ===\n")

output <- list(
  strength_model = results$strength$model,
  hypertrophy_model = results$hypertrophy$model,
  strength_summary = strength_summary,
  hypertrophy_summary = hypertrophy_summary,
  data_strength = data$strength,
  data_hypertrophy = data$hypertrophy,
  correlation_strength = data$correlation_strength$imputed_correlation,
  correlation_hypertrophy = data$correlation_hypertrophy$imputed_correlation,
  loader_summary = summary_stats
)

saveRDS(output, OUTPUT_PATH)
cat("Results saved to:", OUTPUT_PATH, "\n")

cat("\n=== REPLICATION COMPLETE ===\n")
