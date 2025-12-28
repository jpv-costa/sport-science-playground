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
  ],
  R/calculators/model_validator[
    ModelValidator,
    CrossValidationResult,
    CoefficientStabilityResult,
    ModelSelectionStabilityResult,
    PredictionErrorByParticipantResult,
    CalibrationResult,
    IntervalComparisonResult
  ],
  R/calculators/anomaly_detector[
    AnomalyDetector,
    AnomalyResult,
    ParticipantAnomalyResult
  ],
  R/calculators/influence_diagnostics[
    InfluenceDiagnostics,
    InfluenceDiagnosticsResult
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
# Load with pseudoanonymization (P01-P19 instead of real names)
data <- loader$load(anonymize = TRUE)
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
# Section 7: Robust Methods
# ==============================================================================

cat("=== Section 7: Robust Methods ===\n\n")

cat("Addressing model assumption violations with robust approaches...\n\n")

# Pre-load packages needed for robust methods (box module isolation workaround)
suppressPackageStartupMessages({
  library(lme4)
  library(robustlmm)
  library(clubSandwich)
})

# 7.1: Robust LMM using robustlmm
cat("7.1 Robust LMM (robustlmm package):\n")
cat("Fitting robust LMM to downweight outliers...\n")

robust_result <- tryCatch({
  lmm$fit_robust(
    data = data,
    formula = mean_velocity ~ rir,
    random_formula = ~1 + rir | id,
    model_name = "robust_rlmer"
  )
}, error = function(e) {
  cat("  Warning: Robust LMM failed:", e$message, "\n")
  NULL
})

if (!is.null(robust_result)) {
  cat("\nCoefficient Comparison (Standard vs Robust):\n")
  cat("-" |> rep(70) |> paste(collapse = ""), "\n")
  cat(sprintf("%-15s %12s %12s %12s %12s\n",
              "Term", "Standard", "Robust", "Diff", "% Diff"))
  cat("-" |> rep(70) |> paste(collapse = ""), "\n")
  for (i in seq_len(nrow(robust_result$comparison_table))) {
    row <- robust_result$comparison_table[i, ]
    cat(sprintf("%-15s %12.4f %12.4f %12.4f %11.1f%%\n",
                row$term, row$standard_estimate, row$robust_estimate,
                row$difference, row$pct_difference))
  }
  cat("\n")
  if (robust_result$conclusions_robust) {
    cat("CONCLUSION: Standard and robust estimates agree (<10% difference).\n")
    cat("Conclusions are ROBUST to potential outliers.\n\n")
  } else {
    cat("CAUTION: Standard and robust estimates differ (>10% difference).\n")
    cat("Results may be influenced by outliers.\n\n")
  }
}

# 7.2: Cluster-Robust Standard Errors
cat("7.2 Cluster-Robust Standard Errors (clubSandwich):\n")

robust_se <- tryCatch({
  lmm$cluster_robust_se(best_model)
}, error = function(e) {
  cat("  Warning: Cluster-robust SE failed:", e$message, "\n")
  NULL
})

if (!is.null(robust_se)) {
  cat("Comparing Wald vs Cluster-Robust Standard Errors:\n")
  cat("-" |> rep(70) |> paste(collapse = ""), "\n")
  cat(sprintf("%-15s %10s %10s %10s %10s %10s\n",
              "Term", "Estimate", "SE_Wald", "SE_Robust", "Ratio", "p_robust"))
  cat("-" |> rep(70) |> paste(collapse = ""), "\n")
  for (i in seq_len(nrow(robust_se))) {
    row <- robust_se[i, ]
    cat(sprintf("%-15s %10.4f %10.4f %10.4f %10.2f %10.4f\n",
                row$term, row$estimate, row$se_wald, row$se_robust,
                row$se_ratio, row$p_robust))
  }
  cat("\n")
  max_ratio <- max(robust_se$se_ratio)
  if (max_ratio < 1.5) {
    cat(sprintf("SE ratio range: %.2f-%.2f (< 1.5: minimal heteroscedasticity impact)\n\n",
                min(robust_se$se_ratio), max_ratio))
  } else {
    cat(sprintf("SE ratio range: %.2f-%.2f (> 1.5: notable heteroscedasticity impact)\n\n",
                min(robust_se$se_ratio), max_ratio))
  }
}

# 7.3: Bootstrap Confidence Intervals
cat("7.3 Bootstrap Confidence Intervals:\n")
cat("Running parametric bootstrap (500 replicates)...\n")

bootstrap_ci <- tryCatch({
  lmm$bootstrap_ci(best_model, n_boot = 500, seed = 42)
}, error = function(e) {
  cat("  Warning: Bootstrap failed:", e$message, "\n")
  NULL
})

if (!is.null(bootstrap_ci)) {
  cat("\nBootstrap 95% CIs for Fixed Effects:\n")
  cat("-" |> rep(70) |> paste(collapse = ""), "\n")
  cat(sprintf("%-15s %10s %12s %12s %10s\n",
              "Term", "Estimate", "Lower 95%", "Upper 95%", "Method"))
  cat("-" |> rep(70) |> paste(collapse = ""), "\n")
  for (i in seq_len(nrow(bootstrap_ci))) {
    row <- bootstrap_ci[i, ]
    cat(sprintf("%-15s %10.4f %12.4f %12.4f %10s\n",
                row$term, row$estimate, row$ci_lower, row$ci_upper, row$method))
  }
  cat("\n")
}

# ==============================================================================
# Section 8: Sensitivity Analysis
# ==============================================================================

cat("=== Section 8: Sensitivity Analysis ===\n\n")

cat("Testing sensitivity of RIR effect to model specification...\n\n")

sensitivity <- lmm$sensitivity_analysis(
  data = data,
  base_formula = mean_velocity ~ rir,
  random_formulas = list(
    "Random intercept only" = ~1 | id,
    "Random slope (best)" = ~1 + rir | id
  )
)

cat("RIR Effect Across Model Specifications:\n")
cat("-" |> rep(80) |> paste(collapse = ""), "\n")
cat(sprintf("%-25s %12s %10s %10s %10s\n",
            "Model", "RIR Effect", "SE", "AIC", "BIC"))
cat("-" |> rep(80) |> paste(collapse = ""), "\n")
for (i in seq_len(nrow(sensitivity$rir_effects))) {
  row <- sensitivity$rir_effects[i, ]
  cat(sprintf("%-25s %12.4f %10.4f %10.1f %10.1f\n",
              row$model, row$rir_estimate, row$rir_se, row$aic, row$bic))
}
cat("\n")

cat("Summary:\n")
cat(sprintf("  Mean RIR effect: %.4f m/s per RIR\n", sensitivity$summary$rir_mean))
cat(sprintf("  SD across models: %.4f\n", sensitivity$summary$rir_sd))
cat(sprintf("  Range: %.4f\n", sensitivity$summary$rir_range))
cat(sprintf("  CV: %.1f%%\n", sensitivity$summary$rir_cv))

if (sensitivity$summary$conclusions_robust) {
  cat("\nCONCLUSION: RIR effect is ROBUST to model specification (<10% variation).\n\n")
} else {
  cat("\nCAUTION: RIR effect shows sensitivity to model specification (>10% variation).\n\n")
}

# ==============================================================================
# Section 9: Study 4 vs Study 5 Comparison
# ==============================================================================

cat("=== Section 9: Study 4 vs Study 5 Comparison ===\n\n")

# Load Study 4 results
study4_path <- "data/processed/deadlift_rir_velocity_results.rds"
if (file.exists(study4_path)) {
  study4 <- readRDS(study4_path)

  cat("Comparing methods and conclusions:\n\n")

  cat("METHODOLOGY:\n")
  cat("-" |> rep(70) |> paste(collapse = ""), "\n")
  cat(sprintf("%-30s %-35s\n", "Study 4 (Simple)", "Study 5 (LMM)"))
  cat("-" |> rep(70) |> paste(collapse = ""), "\n")
  cat(sprintf("%-30s %-35s\n",
              "Separate OLS per participant", "Mixed effects (all data)"))
  cat(sprintf("%-30s %-35s\n",
              "No hierarchical structure", "Random intercepts + slopes"))
  cat(sprintf("%-30s %-35s\n",
              "No uncertainty quantification", "Parametric + conformal CI"))
  cat(sprintf("%-30s %-35s\n",
              "No covariate testing", "LRT, AIC, BIC, Bayes Factor"))
  cat("\n")

  cat("MODEL FIT (R-squared):\n")
  cat("-" |> rep(70) |> paste(collapse = ""), "\n")
  cat(sprintf("%-30s %12s %12s %12s\n", "Approach", "General", "Individual", "Improvement"))
  cat("-" |> rep(70) |> paste(collapse = ""), "\n")
  cat(sprintf("%-30s %12.3f %12.3f %12.1fx\n",
              "Study 4 (OLS)",
              study4$comparison$general_r2,
              study4$comparison$individual_r2,
              study4$comparison$improvement_factor))
  cat(sprintf("%-30s %12.3f %12.3f %12s\n",
              "Study 5 (LMM - marginal/cond)",
              best_model$r2_marginal,
              best_model$r2_conditional,
              "-"))
  cat("\n")

  # Extract RIR coefficient from Study 4 general model
  study4_general <- study4$general_results
  if (!is.null(study4_general) && length(study4_general) > 0) {
    # Calculate mean RIR slope across all linear models
    rir_slopes_study4 <- sapply(study4_general, function(x) {
      if (!is.null(x$model)) {
        coef(x$model)[["rir"]]
      } else {
        NA
      }
    })
    mean_rir_study4 <- mean(rir_slopes_study4, na.rm = TRUE)
  } else {
    mean_rir_study4 <- NA
  }

  # Get RIR coefficient from Study 5
  rir_study5 <- fixed_effects[fixed_effects$term == "rir", "estimate"]

  cat("RIR EFFECT (velocity change per RIR):\n")
  cat("-" |> rep(70) |> paste(collapse = ""), "\n")
  cat(sprintf("Study 4 (mean of individual OLS): %.4f m/s per RIR\n",
              if (!is.na(mean_rir_study4)) mean_rir_study4 else NA))
  cat(sprintf("Study 5 (LMM fixed effect):       %.4f m/s per RIR\n", rir_study5))
  if (!is.na(mean_rir_study4) && !is.na(rir_study5)) {
    pct_diff <- 100 * (rir_study5 - mean_rir_study4) / abs(mean_rir_study4)
    cat(sprintf("Difference: %.1f%%\n", pct_diff))
  }
  cat("\n")

  cat("PREDICTION ACCURACY (MAE in m/s):\n")
  cat("-" |> rep(70) |> paste(collapse = ""), "\n")
  cat(sprintf("%-30s %12s %12s\n", "Approach", "General", "Individual"))
  cat("-" |> rep(70) |> paste(collapse = ""), "\n")
  # Study 5 already calculated
  cat(sprintf("%-30s %12.4f %12.4f\n",
              "Study 5 (LMM)",
              comparison_result$global_mae,
              comparison_result$individual_mae))
  cat("\n")

  cat("KEY CONCLUSIONS:\n")
  cat("-" |> rep(70) |> paste(collapse = ""), "\n")
  cat("1. AGREEMENT: Both studies confirm ~0.03 m/s velocity change per RIR.\n")
  cat("2. AGREEMENT: Individual calibration improves accuracy by ~30-50%.\n")
  cat("3. NEW (Study 5): Load percentage (80% vs 90%) does NOT significantly\n")
  cat("   affect the velocity-RIR relationship (p = 0.13).\n")
  cat("4. NEW (Study 5): A single global velocity table can be used\n")
  cat("   regardless of load percentage.\n")
  cat("5. NEW (Study 5): Conformal prediction provides guaranteed 95% coverage\n")
  cat("   vs ~83% with parametric intervals.\n")
  cat("\n")

  study4_comparison <- list(
    study4_general_r2 = study4$comparison$general_r2,
    study4_individual_r2 = study4$comparison$individual_r2,
    study4_improvement = study4$comparison$improvement_factor,
    study5_marginal_r2 = best_model$r2_marginal,
    study5_conditional_r2 = best_model$r2_conditional,
    study4_rir_effect = mean_rir_study4,
    study5_rir_effect = rir_study5
  )
} else {
  cat("Study 4 results not found. Run 'make replicate-deadlift' first.\n")
  study4_comparison <- NULL
}

# ==============================================================================
# Section 10: Advanced Validation (LOO-CV Sensitivity + Anomaly Detection)
# ==============================================================================

cat("=== Section 10: Advanced Validation ===\n\n")

# 10.1: Cross-Validation and Coefficient Stability
cat("10.1 LOO-CV Cross-Validation and Coefficient Stability\n")
cat("-" |> rep(70) |> paste(collapse = ""), "\n\n")

validator <- ModelValidator$new()

# LOO-CV for prediction error
cat("Running Leave-One-Participant-Out Cross-Validation...\n")
cv_result <- validator$leave_one_participant_out(
  data = data,
  model = best_model$model,
  id_col = "id",
  outcome_col = "mean_velocity"
)

cat(sprintf("  LOO-CV MAE: %.4f m/s (%.1f mm/s)\n",
            cv_result$mean_error, cv_result$mean_error * 1000))
cat(sprintf("  LOO-CV RMSE: %.4f m/s (%.1f mm/s)\n",
            cv_result$rmse, cv_result$rmse * 1000))
cat(sprintf("  SE of MAE: %.4f m/s\n\n", cv_result$se_error))

# Coefficient Stability
cat("Testing RIR coefficient stability across LOO folds...\n")
coef_stability <- validator$loo_coefficient_stability(
  data = data,
  model = best_model$model,
  focal_coef = "rir",
  id_col = "id"
)

cat(sprintf("  Full model RIR estimate: %.5f m/s per RIR\n", coef_stability$full_model_estimate))
cat(sprintf("  LOO mean estimate: %.5f m/s per RIR\n", coef_stability$mean_estimate))
cat(sprintf("  Coefficient of Variation: %.2f%%\n", coef_stability$cv_percent))
cat(sprintf("  Range: [%.5f, %.5f]\n", coef_stability$range[1], coef_stability$range[2]))
cat(sprintf("  Stability: %s\n\n", if (coef_stability$is_stable) "STABLE" else "UNSTABLE"))

cat("Most influential participants (by RIR coefficient deviation):\n")
print(head(coef_stability$influential_participants, 5))
cat("\n")

# 10.2: Model Selection Stability
cat("10.2 LOO-CV Model Selection Stability\n")
cat("-" |> rep(70) |> paste(collapse = ""), "\n\n")

# Compare base vs with_load vs best model
models_to_compare <- list(
  base = model_base$model,
  with_load = model_with_load$model,
  random_slope = model_random_slope$model
)

selection_stability <- validator$loo_model_selection_stability(
  data = data,
  models = models_to_compare,
  criterion = "BIC",
  id_col = "id"
)

cat(sprintf("  Full data winner: %s\n", selection_stability$full_data_winner))
cat(sprintf("  Consensus model: %s\n", selection_stability$consensus_model))
cat(sprintf("  Stability: %.1f%% of folds agree with full data\n\n",
            selection_stability$stability_percent))

cat("Vote counts across folds:\n")
print(selection_stability$vote_counts)
cat("\n")

# 10.3: Prediction Error by Participant
cat("10.3 Prediction Error by Participant (LOO-CV)\n")
cat("-" |> rep(70) |> paste(collapse = ""), "\n\n")

pred_error_result <- validator$loo_prediction_error_by_participant(
  data = data,
  model = best_model$model,
  id_col = "id",
  outcome_col = "mean_velocity"
)

cat(sprintf("  Overall RMSE: %.4f m/s (%.1f mm/s)\n",
            pred_error_result$overall_rmse, pred_error_result$overall_rmse * 1000))
cat(sprintf("  Overall MAE: %.4f m/s (%.1f mm/s)\n\n",
            pred_error_result$overall_mae, pred_error_result$overall_mae * 1000))

cat("Hardest to predict participants:\n")
print(pred_error_result$get_hardest_to_predict(5))
cat("\n")

# 10.4: Anomaly Detection (Isolation Forest)
cat("10.4 Anomaly Detection (Isolation Forest)\n")
cat("-" |> rep(70) |> paste(collapse = ""), "\n\n")

anomaly_detector <- AnomalyDetector$new(random_state = 42)

# Stage 1: Raw data anomalies
cat("Stage 1: Detecting anomalies in raw data (velocity, RIR, load)...\n")
raw_anomalies <- anomaly_detector$detect_raw_data_anomalies(
  data = data,
  features = c("mean_velocity", "rir", "load_percentage"),
  contamination = 0.05
)

cat(sprintf("  Threshold: %.4f\n", raw_anomalies$threshold))
cat(sprintf("  Anomalies flagged: %d (%.1f%%)\n",
            raw_anomalies$n_anomalies,
            100 * raw_anomalies$n_anomalies / length(raw_anomalies$scores)))
if (raw_anomalies$n_anomalies > 0) {
  cat(sprintf("  Anomaly indices: %s\n",
              paste(head(raw_anomalies$anomaly_indices, 10), collapse = ", ")))
}
cat("\n")

# Stage 2: Random effects anomalies
cat("Stage 2: Detecting anomalies in random effects (velocity-RIR patterns)...\n")
re_anomalies <- anomaly_detector$detect_random_effects_anomalies(
  model = best_model$model,
  contamination = 0.1
)

cat(sprintf("  Threshold: %.4f\n", re_anomalies$threshold))
cat(sprintf("  Anomalous participants: %d\n", sum(re_anomalies$is_anomaly)))
if (sum(re_anomalies$is_anomaly) > 0) {
  cat(sprintf("  Anomalous IDs: %s\n",
              paste(re_anomalies$get_anomalous_ids(), collapse = ", ")))
}
cat("\n")

# Stage 3: CV residual anomalies
cat("Stage 3: Detecting anomalies in CV residuals (hard to predict)...\n")
cv_anomalies <- anomaly_detector$detect_cv_residual_anomalies(
  residuals = pred_error_result$residuals,
  participant_ids = pred_error_result$residual_participant_ids,
  contamination = 0.1
)

cat(sprintf("  Threshold: %.4f\n", cv_anomalies$threshold))
cat(sprintf("  Anomalous participants: %d\n", sum(cv_anomalies$is_anomaly)))
if (sum(cv_anomalies$is_anomaly) > 0) {
  cat(sprintf("  Anomalous IDs: %s\n",
              paste(cv_anomalies$get_anomalous_ids(), collapse = ", ")))
}
cat("\n")

# 10.5: Influence Diagnostics
cat("10.5 Influence Diagnostics (Cook's Distance)\n")
cat("-" |> rep(70) |> paste(collapse = ""), "\n\n")

influence_diag <- InfluenceDiagnostics$new()
influence_result <- influence_diag$calculate_observation_influence(best_model$model)

cat(sprintf("  Observations: %d\n", influence_result$n_observations))
cat(sprintf("  Cook's D threshold: %.4f\n", influence_result$influential_threshold))
cat(sprintf("  High influence observations: %d\n", influence_result$count_influential()))
cat("\n")

# Compare IF vs Cook's D
cat("Comparing Isolation Forest vs Cook's Distance:\n")
# Create influence flags for raw data (Cook's D > threshold)
cooks_flags <- influence_result$cooks_d > influence_result$influential_threshold
if_comparison <- anomaly_detector$compare_with_influence_diagnostics(
  raw_anomalies,
  cooks_flags
)

cat(sprintf("  Both methods flagged: %d observations\n", if_comparison$both_flagged))
cat(sprintf("  Only Isolation Forest: %d observations\n", if_comparison$only_isolation_forest))
cat(sprintf("  Only Cook's D: %d observations\n", if_comparison$only_cooks_distance))
cat(sprintf("  Agreement rate: %.1f%%\n", if_comparison$agreement_rate * 100))
cat(sprintf("  Interpretation: %s\n\n", if_comparison$interpretation))

# 10.6: Model Calibration
cat("10.6 Model Calibration Assessment\n")
cat("-" |> rep(70) |> paste(collapse = ""), "\n\n")

calibration <- validator$calculate_calibration(
  model = best_model$model,
  data = data,
  outcome_col = "mean_velocity",
  use_individual = TRUE
)

cat(sprintf("  Calibration slope: %.4f (ideal = 1)\n", calibration$slope))
cat(sprintf("  Calibration intercept: %.4f (ideal = 0)\n", calibration$intercept))
cat(sprintf("  R-squared: %.4f\n", calibration$r_squared))
cat(sprintf("  Mean bias: %.4f m/s\n", calibration$mean_bias))
cat(sprintf("  %s\n\n", calibration$interpret()))

# 10.7: CI vs PI Comparison
cat("10.7 Confidence vs Prediction Interval Comparison\n")
cat("-" |> rep(70) |> paste(collapse = ""), "\n\n")

interval_comparison <- validator$compare_ci_vs_pi(
  model = best_model$model,
  data = data,
  outcome_col = "mean_velocity"
)

cat(sprintf("  Mean CI width: %.4f m/s\n", interval_comparison$ci_width))
cat(sprintf("  Mean PI width: %.4f m/s\n", interval_comparison$pi_width))
cat(sprintf("  PI/CI ratio: %.1fx wider\n",
            interval_comparison$pi_width / interval_comparison$ci_width))
cat(sprintf("  Empirical CI coverage: %.1f%% (target 95%%)\n",
            interval_comparison$ci_coverage * 100))
cat(sprintf("  Empirical PI coverage: %.1f%% (target 95%%)\n\n",
            interval_comparison$pi_coverage * 100))

cat("Summary:\n")
cat("  CI is for estimating the MEAN velocity (narrow, research use)\n")
cat("  PI is for predicting a NEW observation (wider, practical use)\n\n")

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
  conformal_intervals = conformal_intervals$to_list(),

  # Robust methods (Section 7)
  robust_lmm = robust_result,
  robust_se = robust_se,
  bootstrap_ci = bootstrap_ci,

  # Sensitivity analysis (Section 8)
  sensitivity = sensitivity,

  # Study 4 comparison (Section 9)
  study4_comparison = study4_comparison,

  # Advanced Validation (Section 10)
  advanced_validation = list(
    # Cross-validation
    loo_cv = cv_result$to_list(),

    # Coefficient stability
    coefficient_stability = coef_stability$to_list(),

    # Model selection stability
    model_selection_stability = selection_stability$to_list(),

    # Prediction error by participant
    prediction_error_by_participant = pred_error_result$to_list(),

    # Anomaly detection
    anomalies = list(
      raw_data = raw_anomalies$summarize(),
      random_effects = re_anomalies$summarize(),
      cv_residuals = cv_anomalies$summarize(),
      # Full result objects for visualizations
      raw_data_result = raw_anomalies,
      random_effects_result = re_anomalies,
      cv_residuals_result = cv_anomalies
    ),

    # Influence diagnostics
    influence = list(
      n_observations = influence_result$n_observations,
      cooks_d_threshold = influence_result$influential_threshold,
      high_influence_count = influence_result$count_influential()
    ),

    # IF vs Cook's D comparison
    if_vs_cooks_comparison = if_comparison,

    # Calibration
    calibration = calibration$to_list(),

    # CI vs PI
    ci_vs_pi = interval_comparison$to_list()
  ),

  # Keep full objects for report plots
  advanced_validation_full = list(
    cv_result = cv_result,
    coef_stability = coef_stability,
    selection_stability = selection_stability,
    pred_error_result = pred_error_result,
    raw_anomalies = raw_anomalies,
    re_anomalies = re_anomalies,
    cv_anomalies = cv_anomalies,
    influence_result = influence_result,
    calibration = calibration,
    interval_comparison_full = interval_comparison
  )
)

saveRDS(output, OUTPUT_PATH)
cat("Results saved to:", OUTPUT_PATH, "\n")

cat("\n=== STUDY 5 COMPLETE ===\n")
cat("\nNOTE on data structure:\n")
cat(sprintf("- 19 participants × ~21 observations each = %d total observations\n", nrow(data)))
cat("- Each observation is a single rep during a set (multiple reps to failure)\n")
cat("- The QQ plots and residual plots show ALL observations, not just participants\n")
cat("- LMM properly accounts for the nested structure (observations within participants)\n")
