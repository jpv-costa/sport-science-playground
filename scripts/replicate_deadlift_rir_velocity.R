# scripts/replicate_deadlift_rir_velocity.R
# Study 4: Deadlift RIR-Velocity Relationship Modeling
#
# Original thesis research on conventional deadlift RIR-velocity relationships.
# Applies the same methodology as Jukic et al. (2024) but for deadlift exercise.
#
# Uses OOP/SOLID architecture with R6 classes:
# - DeadliftRirDataLoader: Data loading and preprocessing (thesis format)
# - RirVelocityModeler: Model fitting and evaluation (reused)

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
  R/loaders/deadlift_rir_data_loader[DeadliftRirDataLoader],
  R/calculators/rir_velocity_modeler[
    RirVelocityModeler,
    RirVelocityModelResult,
    RirVelocityModelSummary
  ]
)

cat("=== Deadlift RIR-Velocity Relationship Study ===\n")
cat("Original thesis research on conventional deadlift\n\n")

# ==============================================================================
# Configuration
# ==============================================================================

DATA_PATH <- "deadlift-study/TESE - Análise estatística (Preliminar).xlsx"
OUTPUT_PATH <- "data/processed/deadlift_rir_velocity_results.rds"

# ==============================================================================
# Data Loading
# ==============================================================================

cat("=== Loading Data ===\n")

loader <- DeadliftRirDataLoader$new(DATA_PATH)
data <- loader$load()
summary_stats <- loader$summarize(data)

cat("Data Summary:\n")
cat("  Observations:", summary_stats$n_observations, "\n")
cat("  Participants:", summary_stats$n_participants, "\n")
cat("    - Male:", round(summary_stats$n_male), "\n")
cat("    - Female:", round(summary_stats$n_female), "\n")
cat("  Load types:", paste(summary_stats$load_types, collapse = ", "), "\n")
cat("  Days:", paste(summary_stats$days, collapse = ", "), "\n")
cat("  Velocity range:", round(summary_stats$velocity_range, 3), "m/s\n")
cat("  RIR range:", summary_stats$rir_range, "\n")
cat("  Weight range:", summary_stats$weight_range, "kg\n\n")

# ==============================================================================
# General Model Analysis
# ==============================================================================

cat("=== General (Pooled) RIR-Velocity Models ===\n\n")

modeler <- RirVelocityModeler$new()

# Fit general models by load
general_results <- modeler$fit_general_by_load(data, "load_percentage")

cat("General Model Results by Load:\n")
cat("-" |> rep(60) |> paste(collapse = ""), "\n")
cat(sprintf("%-10s %-12s %8s %8s %8s\n", "Load", "Model", "R²", "RSE", "N"))
cat("-" |> rep(60) |> paste(collapse = ""), "\n")

for (load in names(general_results)) {
  for (model_type in c("linear", "polynomial")) {
    result <- general_results[[load]][[model_type]]
    cat(sprintf("%-10s %-12s %8.3f %8.3f %8d\n",
                load, model_type, result$r_squared, result$rse, result$n))
  }
}
cat("\n")

# ==============================================================================
# Individual Model Analysis
# ==============================================================================

cat("=== Individual RIR-Velocity Models ===\n\n")

# Fit individual models
individual_results <- modeler$fit_individual(data, "id", "load_percentage")

# Summarize by model type
linear_summary <- modeler$summarize_individual(individual_results, "linear")
poly_summary <- modeler$summarize_individual(individual_results, "polynomial")

cat("Individual Model Summary (across all participants and loads):\n")
cat("-" |> rep(60) |> paste(collapse = ""), "\n")
cat(sprintf("%-12s %8s %8s %8s %8s %8s\n", "Model", "N", "Med R²", "Min R²", "Max R²", "Med RSE"))
cat("-" |> rep(60) |> paste(collapse = ""), "\n")
cat(sprintf("%-12s %8d %8.3f %8.3f %8.3f %8.3f\n",
            "Linear",
            linear_summary$n_models,
            linear_summary$median_r_squared,
            linear_summary$min_r_squared,
            linear_summary$max_r_squared,
            linear_summary$median_rse))
cat(sprintf("%-12s %8d %8.3f %8.3f %8.3f %8.3f\n",
            "Polynomial",
            poly_summary$n_models,
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

cat(sprintf("Day 1 observations: %d\n", nrow(day1_data)))
cat(sprintf("Day 2 observations: %d\n\n", nrow(day2_data)))

# General prediction (using pooled Day 1 data)
general_accuracy <- modeler$calculate_prediction_accuracy(
  day1_data, day2_data, model_type = "polynomial"
)

cat("General Model Prediction (Day 1 → Day 2):\n")
cat(sprintf("  Mean absolute error: %.2f reps\n", mean(general_accuracy$absolute_error, na.rm = TRUE)))
cat(sprintf("  Median absolute error: %.2f reps\n", stats::median(general_accuracy$absolute_error, na.rm = TRUE)))
cat(sprintf("  %% within 1 rep: %.1f%%\n",
            mean(general_accuracy$absolute_error <= 1, na.rm = TRUE) * 100))
cat(sprintf("  %% within 2 reps: %.1f%%\n",
            mean(general_accuracy$absolute_error <= 2, na.rm = TRUE) * 100))
cat("\n")

# ==============================================================================
# Comparison with Squat Study (Study 3)
# ==============================================================================

cat("=== Comparison: Deadlift vs Squat (Study 3) ===\n\n")

cat("                      Squat (Jukic)    Deadlift (Thesis)\n")
cat("-" |> rep(60) |> paste(collapse = ""), "\n")
cat(sprintf("Participants:          %d               %d\n", 46, summary_stats$n_participants))
cat(sprintf("Loads:                 70/80/90%%        80/90%%\n"))
cat(sprintf("General R²:            ~0.50            %.2f\n", general_r2))
cat(sprintf("Individual R² (med):   ~0.88            %.2f\n", individual_r2))
cat(sprintf("Improvement:           ~1.8x            %.2fx\n", individual_r2 / general_r2))
cat("\n")

# ==============================================================================
# Deadlift-Specific Insights
# ==============================================================================

cat("=== Deadlift-Specific Insights ===\n\n")

# Compare velocity ranges between loads
load_80_data <- loader$filter_by_load(data, "80%")
load_90_data <- loader$filter_by_load(data, "90%")

cat("Velocity by Load:\n")
cat(sprintf("  80%% 1RM: %.3f - %.3f m/s (mean: %.3f)\n",
            min(load_80_data$mean_velocity),
            max(load_80_data$mean_velocity),
            mean(load_80_data$mean_velocity)))
cat(sprintf("  90%% 1RM: %.3f - %.3f m/s (mean: %.3f)\n",
            min(load_90_data$mean_velocity),
            max(load_90_data$mean_velocity),
            mean(load_90_data$mean_velocity)))
cat("\n")

# RIR distribution
cat("RIR Distribution:\n")
rir_table <- table(data$rir)
for (r in sort(as.numeric(names(rir_table)))) {
  cat(sprintf("  RIR %d: %d observations (%.1f%%)\n",
              r, rir_table[as.character(r)],
              rir_table[as.character(r)] / nrow(data) * 100))
}
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
  ),
  day1_data = day1_data,
  day2_data = day2_data
)

saveRDS(output, OUTPUT_PATH)
cat("Results saved to:", OUTPUT_PATH, "\n")

cat("\n=== STUDY 4 COMPLETE ===\n")
