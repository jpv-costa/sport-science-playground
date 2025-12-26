# scripts/replicate_rir_velocity_refactored.R
# Replication of Jukic et al. (2024) - RIR-Velocity Relationship Modeling
#
# Uses OOP/SOLID architecture with R6 classes:
# - RirVelocityDataLoader: Data loading and preprocessing
# - RirVelocityModeler: Model fitting and evaluation
#
# Original: https://doi.org/10.1113/EP091747

# ==============================================================================
# Setup
# ==============================================================================

# Ensure we're in project root for box module resolution
if (file.exists("deadlift-study.Rproj")) {
  # Already in project root
} else if (file.exists("../deadlift-study.Rproj")) {
  setwd("..")
}

options(box.path = getwd())

box::use(
  R/loaders/rir_velocity_data_loader[RirVelocityDataLoader],
  R/calculators/rir_velocity_modeler[
    RirVelocityModeler,
    RirVelocityModelResult,
    RirVelocityModelSummary
  ]
)

cat("=== Jukic et al. (2024) RIR-Velocity Relationship Replication ===\n")
cat("Using refactored OOP/SOLID architecture\n\n")

# ==============================================================================
# Configuration
# ==============================================================================

DATA_PATH <- "data/external/zourdos_rir/rir_velocity_dataset.xlsx"
OUTPUT_PATH <- "data/processed/rir_velocity_replication_results.rds"

# ==============================================================================
# Data Loading
# ==============================================================================

cat("=== Loading Data ===\n")

loader <- RirVelocityDataLoader$new(DATA_PATH)
data <- loader$load()
summary_stats <- loader$summarize(data)

cat("Data Summary:\n")
cat("  Observations:", summary_stats$n_observations, "\n")
cat("  Participants:", summary_stats$n_participants, "\n")
cat("  Load types:", paste(summary_stats$load_types, collapse = ", "), "\n")
cat("  Days:", paste(summary_stats$days, collapse = ", "), "\n")
cat("  Velocity range:", round(summary_stats$velocity_range, 3), "m/s\n")
cat("  RIR range:", summary_stats$rir_range, "\n")
cat("  Mean relative strength:", round(summary_stats$mean_relative_strength, 2), "x BW\n\n")

# ==============================================================================
# General Model Analysis
# ==============================================================================

cat("=== General (Pooled) RIR-Velocity Models ===\n\n")

modeler <- RirVelocityModeler$new()

# Fit general models by load
general_results <- modeler$fit_general_by_load(data, "set_type")

cat("General Model Results by Load:\n")
cat("-" |> rep(60) |> paste(collapse = ""), "\n")
cat(sprintf("%-10s %-12s %8s %8s\n", "Load", "Model", "R²", "RSE"))
cat("-" |> rep(60) |> paste(collapse = ""), "\n")

for (load in names(general_results)) {
  for (model_type in c("linear", "polynomial")) {
    result <- general_results[[load]][[model_type]]
    cat(sprintf("%-10s %-12s %8.3f %8.3f\n",
                load, model_type, result$r_squared, result$rse))
  }
}
cat("\n")

# ==============================================================================
# Individual Model Analysis
# ==============================================================================

cat("=== Individual RIR-Velocity Models ===\n\n")

# Fit individual models
individual_results <- modeler$fit_individual(data, "id", "set_type")

# Summarize by model type
linear_summary <- modeler$summarize_individual(individual_results, "linear")
poly_summary <- modeler$summarize_individual(individual_results, "polynomial")

cat("Individual Model Summary (across all participants and loads):\n")
cat("-" |> rep(60) |> paste(collapse = ""), "\n")
cat(sprintf("%-12s %8s %8s %8s %8s\n", "Model", "Med R²", "Min R²", "Max R²", "Med RSE"))
cat("-" |> rep(60) |> paste(collapse = ""), "\n")
cat(sprintf("%-12s %8.3f %8.3f %8.3f %8.3f\n",
            "Linear",
            linear_summary$median_r_squared,
            linear_summary$min_r_squared,
            linear_summary$max_r_squared,
            linear_summary$median_rse))
cat(sprintf("%-12s %8.3f %8.3f %8.3f %8.3f\n",
            "Polynomial",
            poly_summary$median_r_squared,
            poly_summary$min_r_squared,
            poly_summary$max_r_squared,
            poly_summary$median_rse))
cat("\n")

# ==============================================================================
# Comparison: General vs Individual
# ==============================================================================

cat("=== General vs Individual Model Comparison ===\n\n")

# Average R² for general models (polynomial)
general_r2 <- mean(sapply(general_results, function(x) x$polynomial$r_squared), na.rm = TRUE)
individual_r2 <- poly_summary$median_r_squared

cat(sprintf("General polynomial R² (average across loads): %.3f\n", general_r2))
cat(sprintf("Individual polynomial R² (median): %.3f\n", individual_r2))
cat(sprintf("Improvement factor: %.2fx\n\n", individual_r2 / general_r2))

# ==============================================================================
# Predictive Validity: Day 1 → Day 2
# ==============================================================================

cat("=== Predictive Validity (Day 1 → Day 2) ===\n\n")

day1_data <- loader$filter_by_day(data, "Day 1")
day2_data <- loader$filter_by_day(data, "Day 2")

# General prediction (using pooled Day 1 data)
general_accuracy <- modeler$calculate_prediction_accuracy(
  day1_data, day2_data, model_type = "polynomial"
)

cat("General Model Prediction (Day 1 → Day 2):\n")
cat(sprintf("  Mean absolute error: %.2f reps\n", mean(general_accuracy$absolute_error, na.rm = TRUE)))
cat(sprintf("  Median absolute error: %.2f reps\n", stats::median(general_accuracy$absolute_error, na.rm = TRUE)))
cat(sprintf("  %% within 2 reps: %.1f%%\n",
            mean(general_accuracy$absolute_error <= 2, na.rm = TRUE) * 100))
cat("\n")

# ==============================================================================
# Paper Comparison
# ==============================================================================

cat("=== PAPER COMPARISON ===\n\n")

cat("General Model R² (polynomial):\n")
cat("  Paper claim: ~0.50-0.60 across loads\n")
cat(sprintf("  Our result: %.2f (average)\n", general_r2))
cat("  Match: ", ifelse(general_r2 > 0.4 && general_r2 < 0.7, "YES", "CHECK"), "\n\n")

cat("Individual Model R² (polynomial):\n")
cat("  Paper claim: ~0.85-0.95 (median)\n")
cat(sprintf("  Our result: %.2f (median)\n", individual_r2))
cat("  Match: ", ifelse(individual_r2 > 0.7, "YES", "CHECK"), "\n\n")

cat("Improvement Factor:\n")
cat("  Paper claim: Individual ~2x better fit than general\n")
cat(sprintf("  Our result: %.2fx\n", individual_r2 / general_r2))
cat("  Match: ", ifelse(individual_r2 / general_r2 > 1.5, "YES", "CHECK"), "\n\n")

cat("Prediction Accuracy:\n")
cat("  Paper claim: <2 rep error for individual models\n")
cat(sprintf("  Our result: %.2f reps (general model)\n",
            mean(general_accuracy$absolute_error, na.rm = TRUE)))
cat("\n")

# ==============================================================================
# Save Results
# ==============================================================================

cat("=== Saving Results ===\n")

output <- list(
  data = data,
  summary = summary_stats,
  general_results = lapply(general_results, function(load) {
    lapply(load, function(m) m$to_list())
  }),
  linear_summary = linear_summary$to_list(),
  polynomial_summary = poly_summary$to_list(),
  general_prediction_accuracy = general_accuracy,
  comparison = list(
    general_r2 = general_r2,
    individual_r2 = individual_r2,
    improvement_factor = individual_r2 / general_r2
  )
)

saveRDS(output, OUTPUT_PATH)
cat("Results saved to:", OUTPUT_PATH, "\n")

cat("\n=== REPLICATION COMPLETE ===\n")
