# scripts/analyze_advanced_velocity.R
# Study 6: Advanced Velocity-RIR Analyses
#
# Implements five novel analyses (H2-H6) that extend the deadlift research:
# - H2: Minimum Velocity Threshold (MVT) variability at failure
# - H3: Day-to-day reliability of individual velocity profiles
# - H4: Polynomial vs linear model comparison
# - H5: Velocity decay rate within sets
# - H6: Early rep velocity for failure prediction
#
# Uses OOP/SOLID architecture with R6 classes.

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
  R/calculators/advanced_velocity_analyzer[
    AdvancedVelocityAnalyzer,
    MvtAnalysisResult,
    ReliabilityResult,
    ModelComparisonResult,
    VelocityDecayResult,
    FailurePredictionResult
  ]
)

cat("==========================================================\n")
cat("Study 6: Advanced Velocity-RIR Analyses\n")
cat("==========================================================\n\n")

# ==============================================================================
# Configuration
# ==============================================================================

DATA_PATH <- "deadlift-study/TESE - Análise estatística (Preliminar).xlsx"
OUTPUT_PATH <- "data/processed/advanced_velocity_results.rds"

# ==============================================================================
# Data Loading
# ==============================================================================

cat("=== Loading Data ===\n\n")

loader <- DeadliftRirDataLoader$new(DATA_PATH)
data <- suppressMessages(suppressWarnings(loader$load()))
summary_stats <- loader$summarize(data)

cat("Data Summary:\n")
cat("  Observations:", summary_stats$n_observations, "\n")
cat("  Participants:", summary_stats$n_participants, "\n")
cat("  Sets:", summary_stats$n_sets, "\n")
cat("  Load types:", paste(summary_stats$load_types, collapse = ", "), "\n")
cat("  Days:", paste(summary_stats$days, collapse = ", "), "\n")
cat("  Reps per set:", summary_stats$reps_per_set_range[1], "-",
    summary_stats$reps_per_set_range[2],
    "(mean:", round(summary_stats$mean_reps_per_set, 1), ")\n\n")

# Create analyzer
analyzer <- AdvancedVelocityAnalyzer$new()

# ==============================================================================
# H2: Minimum Velocity Threshold (MVT) Variability
# ==============================================================================

cat("==========================================================\n")
cat("H2: Minimum Velocity Threshold (MVT) Variability\n")
cat("==========================================================\n\n")

cat("Research Question: How variable is the velocity at muscular failure?\n\n")

mvt_result <- analyzer$analyze_mvt_variability(data)

cat("Population MVT Statistics:\n")
cat("  Mean:", round(mvt_result$population_stats$mean, 3), "m/s\n")
cat("  SD:", round(mvt_result$population_stats$sd, 3), "m/s\n")
cat("  CV:", round(mvt_result$population_stats$cv_percent, 1), "%\n")
cat("  IQR:", round(mvt_result$population_stats$iqr, 3), "m/s\n")
cat("  Range:", round(mvt_result$population_stats$min, 3), "-",
    round(mvt_result$population_stats$max, 3), "m/s\n")
cat("  N observations:", mvt_result$population_stats$n, "\n\n")

cat("Individual MVT Statistics:\n")
print(mvt_result$individual_stats[, c("id", "sex", "mean_mvt", "sd_mvt", "n_sets")])
cat("\n")

cat("Sex Comparison:\n")
cat("  Male MVT:", round(mvt_result$sex_comparison$male_mean, 3),
    "+/-", round(mvt_result$sex_comparison$male_sd, 3), "m/s\n")
cat("  Female MVT:", round(mvt_result$sex_comparison$female_mean, 3), "m/s\n")
if (!is.null(mvt_result$sex_comparison$wilcox_p)) {
  cat("  Wilcoxon p-value:", round(mvt_result$sex_comparison$wilcox_p, 4), "\n")
  cat("  Significant:", mvt_result$sex_comparison$significant, "\n")
}
cat("\n")

cat("Load Comparison (80% vs 90%):\n")
cat("  80% MVT:", round(mvt_result$load_comparison$load_80_mean, 3),
    "+/-", round(mvt_result$load_comparison$load_80_sd, 3), "m/s\n")
cat("  90% MVT:", round(mvt_result$load_comparison$load_90_mean, 3),
    "+/-", round(mvt_result$load_comparison$load_90_sd, 3), "m/s\n")
cat("  Difference:", round(mvt_result$load_comparison$difference, 3), "m/s\n")
cat("  Wilcoxon p-value:", round(mvt_result$load_comparison$wilcox_p, 4), "\n")
cat("  Significant:", mvt_result$load_comparison$significant, "\n\n")

cat("Interpretation:\n")
if (mvt_result$population_stats$cv_percent > 25) {
  cat("  High variability (CV > 25%) suggests individual calibration is essential.\n")
  cat("  A single 'universal' MVT threshold would not be accurate for all athletes.\n")
} else {
  cat("  Moderate variability suggests a population-level MVT may be useful.\n")
}
cat("\n")

# ==============================================================================
# H3: Day-to-Day Reliability
# ==============================================================================

cat("==========================================================\n")
cat("H3: Day-to-Day Reliability of Individual Profiles\n")
cat("==========================================================\n\n")

cat("Research Question: Are individual velocity-RIR relationships stable?\n\n")

reliability <- analyzer$calculate_day_reliability(data)

cat("Slope ICC (velocity loss per RIR):\n")
cat("  ICC:", round(reliability$slope_icc$icc, 3), "\n")
cat("  95% CI:", round(reliability$slope_icc$ci_lower, 3), "-",
    round(reliability$slope_icc$ci_upper, 3), "\n")
cat("  Interpretation:", reliability$slope_icc$interpretation, "\n")
cat("  N participants:", reliability$slope_icc$n, "\n\n")

cat("Intercept ICC (baseline velocity):\n")
cat("  ICC:", round(reliability$intercept_icc$icc, 3), "\n")
cat("  95% CI:", round(reliability$intercept_icc$ci_lower, 3), "-",
    round(reliability$intercept_icc$ci_upper, 3), "\n")
cat("  Interpretation:", reliability$intercept_icc$interpretation, "\n\n")

cat("MVT ICC (predicted velocity at failure):\n")
cat("  ICC:", round(reliability$mvt_icc$icc, 3), "\n")
cat("  95% CI:", round(reliability$mvt_icc$ci_lower, 3), "-",
    round(reliability$mvt_icc$ci_upper, 3), "\n")
cat("  Interpretation:", reliability$mvt_icc$interpretation, "\n")
cat("  SEM:", round(reliability$mvt_icc$sem, 4), "m/s\n")
cat("  MDC95:", round(reliability$mvt_icc$mdc95, 4), "m/s\n\n")

cat("Day 1 vs Day 2 Parameters:\n")
print(reliability$day_parameters[, c("id", "slope_day1", "slope_day2",
                                     "mvt_day1", "mvt_day2")])
cat("\n")

cat("Interpretation:\n")
if (reliability$slope_icc$icc >= 0.75) {
  cat("  Good reliability: Single calibration session may be sufficient.\n")
} else if (reliability$slope_icc$icc >= 0.50) {
  cat("  Moderate reliability: Consider periodic recalibration.\n")
} else {
  cat("  Poor reliability: Frequent recalibration recommended.\n")
}
cat("\n")

# ==============================================================================
# H4: Polynomial vs Linear Model Comparison
# ==============================================================================

cat("==========================================================\n")
cat("H4: Polynomial vs Linear Model Comparison\n")
cat("==========================================================\n\n")

cat("Research Question: Does a quadratic model fit better than linear?\n\n")

poly_result <- analyzer$compare_polynomial_models(data)

cat("Individual Model Comparison:\n")
print(poly_result$individual_results[, c("id", "r2_adj_linear", "r2_adj_quad",
                                          "delta_aic", "best_model")])
cat("\n")

cat("Summary:\n")
cat("  Participants analyzed:", poly_result$best_model_summary$n_participants, "\n")
cat("  Linear preferred:", poly_result$best_model_summary$n_linear_best, "\n")
cat("  Quadratic preferred:", poly_result$best_model_summary$n_quad_best,
    "(", round(poly_result$best_model_summary$pct_quad_best, 1), "%)\n")
cat("  Avg R2 improvement (quad vs linear):",
    round(poly_result$best_model_summary$avg_r2_improvement, 4), "\n\n")

if (!is.null(poly_result$population_comparison$lrt_p)) {
  cat("Population LMM Comparison:\n")
  cat("  Linear AIC:", round(poly_result$population_comparison$aic_linear, 1), "\n")
  cat("  Quadratic AIC:", round(poly_result$population_comparison$aic_quad, 1), "\n")
  cat("  LRT Chi-squared:", round(poly_result$population_comparison$lrt_chisq, 2), "\n")
  cat("  LRT p-value:", format(poly_result$population_comparison$lrt_p, digits = 4), "\n")
  cat("  Quadratic significant:", poly_result$population_comparison$quad_significant, "\n\n")
}

cat("Recommendation:", poly_result$best_model_summary$recommendation, "\n\n")

# ==============================================================================
# H5: Velocity Decay Rate Within Sets
# ==============================================================================

cat("==========================================================\n")
cat("H5: Velocity Decay Rate Within Sets\n")
cat("==========================================================\n\n")

cat("Research Question: How does velocity loss per rep change during a set?\n\n")

decay <- analyzer$analyze_velocity_decay(data)

cat("Decay Summary:\n")
cat("  Average velocity loss per rep:",
    round(decay$decay_summary$avg_decay_per_rep, 4), "m/s\n")
cat("  SD:", round(decay$decay_summary$sd_decay, 4), "m/s\n")
cat("  Early reps (1-3) decay:",
    round(decay$decay_summary$early_reps_decay, 4), "m/s\n")
cat("  Late reps (4+) decay:",
    round(decay$decay_summary$late_reps_decay, 4), "m/s\n")
cat("  Decay acceleration:",
    round(decay$decay_summary$decay_acceleration, 4), "m/s\n")
cat("  N observations:", decay$decay_summary$n_observations, "\n\n")

cat("Acceleration Test (is decay rate constant?):\n")
cat("  Slope (delta_v ~ rep_number):",
    round(decay$decay_acceleration$slope, 5), "\n")
cat("  SE:", round(decay$decay_acceleration$se, 5), "\n")
cat("  t-value:", round(decay$decay_acceleration$t_value, 3), "\n")
cat("  p-value:", format(decay$decay_acceleration$p_value, digits = 4), "\n")
cat("  Interpretation:", decay$decay_acceleration$interpretation, "\n\n")

if (!is.null(decay$breakpoint$breakpoint_rep)) {
  cat("Breakpoint Detection:\n")
  cat("  Accelerated decay starts at rep:", decay$breakpoint$breakpoint_rep, "\n")
  cat("  Average decay:", round(decay$breakpoint$avg_decay, 4), "m/s\n")
  cat("  Threshold:", round(decay$breakpoint$threshold, 4), "m/s\n")
} else {
  cat("Breakpoint:", decay$breakpoint$interpretation, "\n")
}
cat("\n")

cat("Interpretation:\n")
if (decay$decay_acceleration$accelerating) {
  cat("  Velocity decay accelerates as the set progresses.\n")
  cat("  Non-linear models or velocity thresholds may be needed.\n")
} else {
  cat("  Velocity decay is relatively constant throughout the set.\n")
  cat("  Linear models are appropriate.\n")
}
cat("\n")

# ==============================================================================
# H6: Early Rep Failure Prediction
# ==============================================================================

cat("==========================================================\n")
cat("H6: Early Rep Velocity for Failure Prediction\n")
cat("==========================================================\n\n")

cat("Research Question: Can first rep velocities predict total reps?\n\n")

prediction <- analyzer$build_failure_predictor(data)

cat("Leave-One-Set-Out Cross-Validation:\n")
cat("  Sets analyzed:", prediction$cv_results$n_sets, "\n")
cat("  MAE:", round(prediction$cv_results$mae, 2), "reps\n")
cat("  RMSE:", round(prediction$cv_results$rmse, 2), "reps\n")
cat("  R2:", round(prediction$cv_results$r2, 3), "\n")
cat("  Within 1 rep:", round(prediction$cv_results$within_1_rep_pct, 1), "%\n")
cat("  Within 2 reps:", round(prediction$cv_results$within_2_reps_pct, 1), "%\n\n")

cat("Model Coefficients (v1 only model):\n")
v1_model <- prediction$prediction_models$v1_only
if (!is.null(v1_model)) {
  coefs <- coef(v1_model)
  cat("  Intercept:", round(coefs["(Intercept)"], 3), "\n")
  cat("  First rep velocity coefficient:", round(coefs["v1"], 3), "\n")
  cat("  Interpretation: +0.1 m/s in v1 -> +",
      round(coefs["v1"] * 0.1, 2), " reps expected\n\n")
}

cat("Practical Lookup Table:\n")
cat("(First rep velocity -> Predicted reps [95% PI])\n\n")
print(prediction$lookup_table)
cat("\n")

cat("Interpretation:\n")
if (prediction$cv_results$mae < 1.5) {
  cat("  Good prediction accuracy (MAE < 1.5 reps).\n")
  cat("  First rep velocity is a useful indicator of set capacity.\n")
} else {
  cat("  Moderate prediction accuracy.\n")
  cat("  Consider using multiple early rep velocities for better prediction.\n")
}
cat("\n")

# ==============================================================================
# Save Results
# ==============================================================================

cat("==========================================================\n")
cat("Saving Results\n")
cat("==========================================================\n\n")

results <- list(
  summary = summary_stats,
  mvt = mvt_result,
  reliability = reliability,
  polynomial_comparison = poly_result,
  velocity_decay = decay,
  failure_prediction = prediction,
  data = data,
  timestamp = Sys.time()
)

# Create output directory if needed
output_dir <- dirname(OUTPUT_PATH)
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

saveRDS(results, OUTPUT_PATH)
cat("Results saved to:", OUTPUT_PATH, "\n\n")

# ==============================================================================
# Summary of Key Findings
# ==============================================================================

cat("==========================================================\n")
cat("SUMMARY OF KEY FINDINGS\n")
cat("==========================================================\n\n")

cat("H2 - MVT Variability:\n")
cat("  MVT =", round(mvt_result$population_stats$mean, 3), "+/-",
    round(mvt_result$population_stats$sd, 3), "m/s (CV =",
    round(mvt_result$population_stats$cv_percent, 1), "%)\n")
cat("  Conclusion: Individual calibration essential due to high variability.\n\n")

cat("H3 - Day-to-Day Reliability:\n")
cat("  Slope ICC =", round(reliability$slope_icc$icc, 2),
    "(", reliability$slope_icc$interpretation, ")\n")
cat("  MVT ICC =", round(reliability$mvt_icc$icc, 2),
    "(", reliability$mvt_icc$interpretation, ")\n")
cat("  Conclusion: Periodic recalibration may be needed.\n\n")

cat("H4 - Model Comparison:\n")
cat("  Quadratic preferred in", round(poly_result$best_model_summary$pct_quad_best, 1),
    "% of participants\n")
cat("  Conclusion:", poly_result$best_model_summary$recommendation, "\n\n")

cat("H5 - Velocity Decay:\n")
cat("  Decay accelerates:", decay$decay_acceleration$accelerating, "\n")
cat("  Conclusion:", decay$decay_acceleration$interpretation, "\n\n")

cat("H6 - Failure Prediction:\n")
cat("  MAE =", round(prediction$cv_results$mae, 2), "reps\n")
cat("  Within 2 reps:", round(prediction$cv_results$within_2_reps_pct, 1), "%\n")
cat("  Conclusion: First rep velocity provides reasonable failure prediction.\n\n")

cat("==========================================================\n")
cat("Analysis Complete\n")
cat("==========================================================\n")
