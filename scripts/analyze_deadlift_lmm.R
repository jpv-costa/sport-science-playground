# scripts/analyze_deadlift_lmm.R
# Study 5: Deadlift LMM Analysis with Velocity Stop Tables
#
# Extends Study 4 with:
# - Linear Mixed Effects Models with interaction effects
# - Model assumption testing (normality, homoscedasticity)
# - Model comparison using AIC, BIC, and Bayes Factor approximation
# - Velocity stop tables (global and individual)
# - Conformal prediction intervals
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
  R/calculators/lmm_analyzer[
    LmmAnalyzer,
    LmmModelResult,
    LmmModelComparison,
    LmmDiagnostics
  ],
  R/calculators/velocity_stop_table_generator[
    VelocityStopTableGenerator,
    VelocityStopTable,
    LoadImportanceResult
  ],
  R/calculators/conformal_predictor[
    ConformalPredictor,
    ConformalPredictionResult,
    CoverageComparisonResult
  ]
)

cat("=== Deadlift LMM Analysis with Velocity Stop Tables ===\n")
cat("Study 5: Linear Mixed Effects Modeling\n\n")

# ==============================================================================
# Configuration
# ==============================================================================

DATA_PATH <- "deadlift-study/TESE - Análise estatística (Preliminar).xlsx"
OUTPUT_PATH <- "data/processed/deadlift_lmm_results.rds"

# ==============================================================================
# Data Loading
# ==============================================================================

cat("=== Loading Data ===\n\n")

loader <- DeadliftRirDataLoader$new(DATA_PATH)
data <- loader$load()
summary_stats <- loader$summarize(data)

cat("Data Summary:\n")
cat("  Observations:", summary_stats$n_observations, "\n")
cat("  Participants:", summary_stats$n_participants, "\n")
cat("  Load types:", paste(summary_stats$load_types, collapse = ", "), "\n")
cat("  Days:", paste(summary_stats$days, collapse = ", "), "\n\n")

# ==============================================================================
# Section 1: LMM Model Building and Comparison
# ==============================================================================

cat("=== Section 1: LMM Model Building ===\n\n")

lmm <- LmmAnalyzer$new()

# Model 1: Base model - only RIR (random intercept)
cat("Fitting Model 1: Base (velocity ~ rir)...\n")
model_base <- lmm$fit(
  data = data,
  formula = mean_velocity ~ rir,
  random_formula = ~1 | id,
  model_name = "Base"
)

# Model 2: Random slope model
cat("Fitting Model 2: Random slope (velocity ~ rir + (1 + rir | id))...\n")
model_random_slope <- lmm$fit(
  data = data,
  formula = mean_velocity ~ rir,
  random_formula = ~1 + rir | id,
  model_name = "Random Slope"
)

# Model 3: Add load percentage
cat("Fitting Model 3: With load (velocity ~ rir + load_percentage)...\n")
model_with_load <- lmm$fit(
  data = data,
  formula = mean_velocity ~ rir + load_percentage,
  random_formula = ~1 + rir | id,
  model_name = "With Load"
)

# Model 4: Add day
cat("Fitting Model 4: With day (velocity ~ rir + load_percentage + day)...\n")
model_with_day <- lmm$fit(
  data = data,
  formula = mean_velocity ~ rir + load_percentage + day,
  random_formula = ~1 + rir | id,
  model_name = "With Day"
)

# Model 5: Interaction: RIR x Load
cat("Fitting Model 5: Interaction (velocity ~ rir * load_percentage)...\n")
model_interaction <- lmm$fit_with_interactions(
  data = data,
  base_formula = mean_velocity ~ rir + load_percentage,
  interactions = c("rir:load_percentage"),
  random_formula = ~1 + rir | id,
  model_name = "RIR x Load Interaction"
)

# Model 6: Full interactions
cat("Fitting Model 6: Full interactions (velocity ~ rir * load_percentage * day)...\n")
model_full <- tryCatch({
  lmm$fit_with_interactions(
    data = data,
    base_formula = mean_velocity ~ rir + load_percentage + day,
    interactions = c("rir:load_percentage", "rir:day", "load_percentage:day"),
    random_formula = ~1 + rir | id,
    model_name = "Full Interactions"
  )
}, error = function(e) {
  cat("  Warning: Full model failed to converge, using simpler model\n")
  model_with_day
})

cat("\n")

# ==============================================================================
# Section 2: Model Comparison
# ==============================================================================

cat("=== Section 2: Model Comparison ===\n\n")

models <- list(
  base = model_base,
  random_slope = model_random_slope,
  with_load = model_with_load,
  with_day = model_with_day,
  interaction = model_interaction
)

# Add full model if it's different from with_day
if (!identical(model_full$model_name, model_with_day$model_name)) {
  models$full <- model_full
}

# Compare models
comparison <- lmm$compare_models(models, criterion = "BIC")

cat("Model Comparison Table:\n")
cat("-" |> rep(100) |> paste(collapse = ""), "\n")
cat(sprintf("%-20s %10s %10s %10s %10s %10s %10s\n",
            "Model", "AIC", "BIC", "delta_AIC", "delta_BIC", "BF_approx", "R2_cond"))
cat("-" |> rep(100) |> paste(collapse = ""), "\n")

for (i in seq_len(nrow(comparison$comparison_table))) {
  row <- comparison$comparison_table[i, ]
  cat(sprintf("%-20s %10.1f %10.1f %10.1f %10.1f %10.3f %10.3f\n",
              row$model,
              row$AIC,
              row$BIC,
              row$delta_AIC,
              row$delta_BIC,
              row$bayes_factor_approx,
              row$R2_conditional))
}
cat("\n")

cat("Best model by BIC:", comparison$best_model_name, "\n\n")

# ==============================================================================
# Section 3: Test Load Importance
# ==============================================================================

cat("=== Section 3: Load Percentage Importance Test ===\n\n")

cat("Question: Does load percentage (80% vs 90%) significantly affect\n")
cat("the velocity-RIR relationship?\n\n")

load_test <- lmm$test_variable_importance(
  data = data,
  base_formula = mean_velocity ~ rir,
  full_formula = mean_velocity ~ rir + load_percentage,
  random_formula = ~1 + rir | id
)

cat("Likelihood Ratio Test:\n")
cat(sprintf("  Chi-squared: %.3f (df = %.0f)\n", load_test$lrt_chisq, load_test$lrt_df))
cat(sprintf("  P-value: %.4f\n", load_test$lrt_p_value))
cat(sprintf("  Significant at alpha=0.05: %s\n", ifelse(load_test$variable_significant, "YES", "NO")))
cat("\n")

cat("Information Criteria:\n")
cat(sprintf("  Delta AIC: %.2f (positive = load improves model)\n", load_test$delta_aic))
cat(sprintf("  Delta BIC: %.2f\n", load_test$delta_bic))
cat("\n")

cat("Bayes Factor Analysis:\n")
cat(sprintf("  BF (simpler vs complex): %.3f\n", load_test$bayes_factor))
cat(sprintf("  Interpretation: %s\n", load_test$bf_interpretation))
cat("\n")

# ==============================================================================
# Section 4: Model Diagnostics (Best Model)
# ==============================================================================

cat("=== Section 4: Model Diagnostics ===\n\n")

best_model <- comparison$best_model
diagnostics <- lmm$test_assumptions(best_model)

cat("Normality of Residuals:\n")
cat(sprintf("  Test: %s\n", diagnostics$normality_test$test_name))
cat(sprintf("  Statistic: %.4f\n", diagnostics$normality_test$statistic))
cat(sprintf("  P-value: %.4f\n", diagnostics$normality_test$p_value))
cat(sprintf("  Conclusion: %s\n", diagnostics$normality_test$interpretation))
cat("\n")

cat("Homoscedasticity (Residuals vs Fitted):\n")
cat(sprintf("  Correlation: %.4f\n", diagnostics$homoscedasticity_test$correlation))
cat(sprintf("  P-value: %.4f\n", diagnostics$homoscedasticity_test$p_value))
cat(sprintf("  Conclusion: %s\n", diagnostics$homoscedasticity_test$interpretation))
cat("\n")

# Fixed effects summary
cat("Fixed Effects (Best Model):\n")
cat("-" |> rep(70) |> paste(collapse = ""), "\n")
fixed_effects <- lmm$get_fixed_effects(best_model)
cat(sprintf("%-25s %10s %10s %10s %12s\n", "Term", "Estimate", "SE", "t-value", "p-value"))
cat("-" |> rep(70) |> paste(collapse = ""), "\n")
for (i in seq_len(nrow(fixed_effects))) {
  row <- fixed_effects[i, ]
  cat(sprintf("%-25s %10.4f %10.4f %10.2f %12.4f\n",
              row$term, row$estimate, row$std_error, row$t_value, row$p_value))
}
cat("\n")

# ==============================================================================
# Section 5: Velocity Stop Tables
# ==============================================================================

cat("=== Section 5: Velocity Stop Tables ===\n\n")

vst_generator <- VelocityStopTableGenerator$new()

# Test if load matters for velocity table
cat("Testing if load percentage affects velocity table recommendations...\n")
load_importance <- vst_generator$test_load_importance(data)

cat(sprintf("  LRT P-value: %.4f\n", load_importance$lrt_p_value))
cat(sprintf("  Bayes Factor: %.3f\n", load_importance$bayes_factor))
cat(sprintf("  Recommendation: %s\n\n", load_importance$recommendation))

# Generate appropriate table based on recommendation
if (load_importance$recommendation == "global") {
  cat("Generating GLOBAL velocity stop table (load-independent):\n\n")
  velocity_table <- vst_generator$generate_global_table(data, rir_targets = 0:7)

  cat("Global Velocity Stop Table:\n")
  cat("-" |> rep(50) |> paste(collapse = ""), "\n")
  cat(sprintf("%-5s %12s %12s %12s\n", "RIR", "Velocity", "Lower 95%", "Upper 95%"))
  cat("-" |> rep(50) |> paste(collapse = ""), "\n")
  for (i in seq_len(nrow(velocity_table$table))) {
    row <- velocity_table$table[i, ]
    cat(sprintf("%-5d %12.3f %12.3f %12.3f\n",
                row$rir, row$velocity, row$lower_95, row$upper_95))
  }

} else {
  cat("Generating LOAD-SPECIFIC velocity stop tables:\n\n")
  velocity_table <- vst_generator$generate_load_specific_table(data, rir_targets = 0:7)

  cat("Load-Specific Velocity Stop Table:\n")
  cat("-" |> rep(60) |> paste(collapse = ""), "\n")
  cat(sprintf("%-5s %10s %12s %12s %12s\n", "RIR", "Load", "Velocity", "Lower 95%", "Upper 95%"))
  cat("-" |> rep(60) |> paste(collapse = ""), "\n")
  for (i in seq_len(nrow(velocity_table$table))) {
    row <- velocity_table$table[i, ]
    cat(sprintf("%-5d %10s %12.3f %12.3f %12.3f\n",
                row$rir, row$load_percentage, row$velocity, row$lower_95, row$upper_95))
  }
}
cat("\n")

# Compare general vs individual
cat("Comparison: General vs Individual Tables:\n")
comparison_result <- vst_generator$compare_general_vs_individual(data)
cat(sprintf("  Global MAE: %.4f m/s\n", comparison_result$global_mae))
cat(sprintf("  Individual MAE: %.4f m/s\n", comparison_result$individual_mae))
cat(sprintf("  Improvement: %.1f%%\n", comparison_result$mae_improvement_pct))
cat(sprintf("  %s\n\n", comparison_result$recommendation))

# ==============================================================================
# Section 6: Conformal Prediction Intervals
# ==============================================================================

cat("=== Section 6: Conformal Prediction Intervals ===\n\n")

# Split data: Day 1 for calibration, Day 2 for testing
day1_data <- loader$filter_by_day(data, "Day 1")
day2_data <- loader$filter_by_day(data, "Day 2")

cat(sprintf("Calibration set (Day 1): %d observations\n", nrow(day1_data)))
cat(sprintf("Test set (Day 2): %d observations\n\n", nrow(day2_data)))

# Fit conformal predictor
cp <- ConformalPredictor$new(alpha = 0.05)
cp$fit(best_model$model, day1_data, response_var = "mean_velocity")

# Get prediction intervals on test set
conformal_intervals <- cp$predict_interval(day2_data)

cat("Conformal Prediction (95% intervals):\n")
cat(sprintf("  Calibrated quantile (q_hat): %.4f m/s\n", cp$get_q_hat()))
cat(sprintf("  Interval width: %.4f m/s (constant across predictions)\n\n",
            mean(conformal_intervals$predictions$interval_width)))

# Calculate coverage
conformal_coverage <- cp$calculate_coverage(day2_data, conformal_intervals)
cat(sprintf("Empirical Coverage on Test Set: %.1f%% (target: 95%%)\n\n", conformal_coverage * 100))

# Compare with parametric intervals
cat("Comparison: Conformal vs Parametric Intervals:\n")
interval_comparison <- cp$compare_with_parametric(best_model, day2_data)

cat("-" |> rep(50) |> paste(collapse = ""), "\n")
cat(sprintf("%-20s %15s %15s\n", "Metric", "Parametric", "Conformal"))
cat("-" |> rep(50) |> paste(collapse = ""), "\n")
cat(sprintf("%-20s %14.1f%% %14.1f%%\n", "Coverage",
            interval_comparison$parametric_coverage * 100,
            interval_comparison$conformal_coverage * 100))
cat(sprintf("%-20s %14.4f %14.4f\n", "Avg Width (m/s)",
            interval_comparison$parametric_width,
            interval_comparison$conformal_width))
cat("-" |> rep(50) |> paste(collapse = ""), "\n\n")

if (abs(interval_comparison$conformal_coverage - 0.95) <
    abs(interval_comparison$parametric_coverage - 0.95)) {
  cat("Conformal intervals achieve coverage closer to 95% target.\n")
} else {
  cat("Parametric intervals achieve coverage closer to 95% target.\n")
}
cat("\n")

# ==============================================================================
# Save Results
# ==============================================================================

cat("=== Saving Results ===\n")

output <- list(
  # Data
  data = data,
  summary = summary_stats,
  day1_data = day1_data,
  day2_data = day2_data,

  # Models
  models = list(
    base = model_base$to_list(),
    random_slope = model_random_slope$to_list(),
    with_load = model_with_load$to_list(),
    with_day = model_with_day$to_list(),
    interaction = model_interaction$to_list(),
    best = best_model$to_list()
  ),
  best_model = best_model,  # Keep full object for report

  # Model comparison
  model_comparison = comparison$to_list(),

  # Load importance test
  load_importance_test = load_test,
  load_importance_result = load_importance$to_list(),

  # Diagnostics
  diagnostics = diagnostics$to_list(),
  diagnostics_full = diagnostics,  # Keep full object for plots
  fixed_effects = fixed_effects,

  # Velocity tables
  velocity_table = velocity_table$to_list(),
  velocity_table_full = velocity_table,
  individual_comparison = comparison_result,

  # Conformal prediction
  conformal = list(
    q_hat = cp$get_q_hat(),
    coverage = conformal_coverage,
    interval_width = mean(conformal_intervals$predictions$interval_width),
    comparison = interval_comparison$to_list()
  ),
  conformal_intervals = conformal_intervals$to_list()
)

saveRDS(output, OUTPUT_PATH)
cat("Results saved to:", OUTPUT_PATH, "\n")

cat("\n=== STUDY 5 COMPLETE ===\n")
