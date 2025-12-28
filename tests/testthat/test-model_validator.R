# tests/testthat/test-model_validator.R
# Tests for ModelValidator, InfluenceDiagnostics, MarginalEffectsAnalyzer, PowerAnalyzer

box::use(
  testthat[describe, it, expect_true, expect_s3_class, expect_type,
           expect_length, expect_gt, expect_equal, expect_lt, skip_if_not_installed]
)

# Set box path relative to project root
options(box.path = c("../../R", getOption("box.path")))

# Load modules under test
box::use(
  calculators/model_validator[
    ModelValidator, CrossValidationResult, CalibrationResult,
    IntervalComparisonResult, CoefficientStabilityResult,
    ModelSelectionStabilityResult, PredictionErrorByParticipantResult,
    leave_one_out_cv, calculate_calibration
  ],
  calculators/influence_diagnostics[
    InfluenceDiagnostics, InfluenceDiagnosticsResult,
    ParticipantInfluenceResult, calculate_influence
  ],
  calculators/marginal_effects_analyzer[
    MarginalEffectsAnalyzer, MarginalEffectResult,
    calculate_marginal_effect
  ],
  calculators/power_analyzer[
    PowerAnalyzer, PowerAnalysisResult, SampleSizeResult,
    calculate_power
  ]
)

# =============================================================================
# TEST FIXTURES
# =============================================================================

#' Create test data mimicking deadlift study structure
create_test_data <- function(n_participants = 15, obs_per_participant = 20) {
  set.seed(42)

  # Individual intercepts and slopes
  intercepts <- rnorm(n_participants, mean = 0.4, sd = 0.05)
  slopes <- rnorm(n_participants, mean = 0.025, sd = 0.008)

  data_list <- lapply(seq_len(n_participants), function(i) {
    rir <- sample(0:7, obs_per_participant, replace = TRUE)
    velocity <- intercepts[i] + slopes[i] * rir + rnorm(obs_per_participant, sd = 0.03)
    data.frame(
      id = factor(rep(i, obs_per_participant)),
      rir = rir,
      mean_velocity = velocity,
      load_pct = rep(runif(1, 70, 90), obs_per_participant)
    )
  })

  do.call(rbind, data_list)
}

#' Create fitted model for testing
create_test_model <- function(data) {
  skip_if_not_installed("lme4")
  lme4::lmer(mean_velocity ~ rir + (1 + rir | id), data = data, REML = FALSE)
}

# =============================================================================
# MODEL VALIDATOR TESTS
# =============================================================================

describe("ModelValidator", {

  describe("initialization", {
    it("creates instance", {
      validator <- ModelValidator$new()
      expect_s3_class(validator, "ModelValidator")
    })
  })

  describe("leave_one_participant_out", {
    skip_if_not_installed("lme4")

    it("performs LOO-CV", {
      data <- create_test_data(n_participants = 10, obs_per_participant = 15)
      model <- create_test_model(data)
      validator <- ModelValidator$new()

      result <- validator$leave_one_participant_out(data, model)

      expect_s3_class(result, "CrossValidationResult")
      expect_equal(result$cv_type, "Leave-One-Participant-Out")
      expect_equal(result$n_folds, 10)
      expect_gt(length(result$fold_errors), 0)
      expect_gt(result$mean_error, 0)
    })
  })

  describe("k_fold_cv", {
    skip_if_not_installed("lme4")

    it("performs k-fold CV", {
      data <- create_test_data(n_participants = 15, obs_per_participant = 15)
      model <- create_test_model(data)
      validator <- ModelValidator$new()

      result <- validator$k_fold_cv(data, model, k = 5)

      expect_s3_class(result, "CrossValidationResult")
      expect_true(grepl("5-Fold", result$cv_type))
      expect_length(result$fold_errors, 5)
    })
  })

  describe("calculate_calibration", {
    skip_if_not_installed("lme4")

    it("calculates calibration metrics", {
      data <- create_test_data()
      model <- create_test_model(data)
      validator <- ModelValidator$new()

      result <- validator$calculate_calibration(model, data)

      expect_s3_class(result, "CalibrationResult")
      expect_true(is.numeric(result$slope))
      expect_true(is.numeric(result$intercept))
      expect_true(is.numeric(result$r_squared))
      # Slope should be close to 1 for well-calibrated model
      expect_lt(abs(result$slope - 1), 0.3)
    })

    it("interprets calibration quality", {
      data <- create_test_data()
      model <- create_test_model(data)
      validator <- ModelValidator$new()

      result <- validator$calculate_calibration(model, data)
      interpretation <- result$interpret()

      expect_type(interpretation, "character")
      expect_true(nchar(interpretation) > 0)
    })
  })

  describe("compare_ci_vs_pi", {
    skip_if_not_installed("lme4")

    it("compares confidence and prediction intervals", {
      data <- create_test_data()
      model <- create_test_model(data)
      validator <- ModelValidator$new()

      result <- validator$compare_ci_vs_pi(model, data)

      expect_s3_class(result, "IntervalComparisonResult")
      expect_gt(result$pi_width, result$ci_width)  # PI should be wider
      expect_gt(result$pi_coverage, 0)
      expect_lt(result$pi_coverage, 1)
    })
  })

  # =========================================================================
  # LOO-CV SENSITIVITY ANALYSIS TESTS
  # =========================================================================

  describe("loo_coefficient_stability", {
    skip_if_not_installed("lme4")

    it("calculates coefficient stability", {
      data <- create_test_data(n_participants = 10, obs_per_participant = 15)
      model <- create_test_model(data)
      validator <- ModelValidator$new()

      result <- validator$loo_coefficient_stability(data, model, "rir", "id")

      expect_s3_class(result, "CoefficientStabilityResult")
      expect_equal(result$focal_coefficient, "rir")
      expect_true(is.numeric(result$full_model_estimate))
      expect_length(result$fold_estimates, 10)
      expect_true(is.numeric(result$cv_percent))
    })

    it("identifies influential participants", {
      data <- create_test_data(n_participants = 10, obs_per_participant = 15)
      model <- create_test_model(data)
      validator <- ModelValidator$new()

      result <- validator$loo_coefficient_stability(data, model, "rir", "id")

      expect_s3_class(result$influential_participants, "data.frame")
      expect_true("participant" %in% names(result$influential_participants))
      expect_true("deviation" %in% names(result$influential_participants))
    })

    it("determines stability status", {
      data <- create_test_data(n_participants = 10, obs_per_participant = 15)
      model <- create_test_model(data)
      validator <- ModelValidator$new()

      result <- validator$loo_coefficient_stability(data, model, "rir", "id")

      expect_type(result$is_stable, "logical")
    })

    it("interprets stability", {
      data <- create_test_data(n_participants = 10, obs_per_participant = 15)
      model <- create_test_model(data)
      validator <- ModelValidator$new()

      result <- validator$loo_coefficient_stability(data, model, "rir", "id")
      interpretation <- result$interpret()

      expect_type(interpretation, "character")
      expect_true(nchar(interpretation) > 0)
      expect_true(grepl("rir", interpretation))
    })
  })

  describe("loo_model_selection_stability", {
    skip_if_not_installed("lme4")

    it("compares model selection across folds", {
      data <- create_test_data(n_participants = 10, obs_per_participant = 15)
      model1 <- lme4::lmer(mean_velocity ~ rir + (1 | id), data = data, REML = FALSE)
      model2 <- create_test_model(data)

      models <- list(
        simple = model1,
        full = model2
      )

      validator <- ModelValidator$new()
      result <- validator$loo_model_selection_stability(data, models, "BIC", "id")

      expect_s3_class(result, "ModelSelectionStabilityResult")
      expect_true(result$full_data_winner %in% c("simple", "full"))
      expect_length(result$fold_winners, 10)
      expect_true(is.numeric(result$stability_percent))
    })

    it("calculates vote counts", {
      data <- create_test_data(n_participants = 10, obs_per_participant = 15)
      model1 <- lme4::lmer(mean_velocity ~ rir + (1 | id), data = data, REML = FALSE)
      model2 <- create_test_model(data)

      models <- list(simple = model1, full = model2)

      validator <- ModelValidator$new()
      result <- validator$loo_model_selection_stability(data, models, "BIC", "id")

      expect_true(inherits(result$vote_counts, "table"))
      expect_equal(sum(result$vote_counts), 10)
    })

    it("interprets selection stability", {
      data <- create_test_data(n_participants = 10, obs_per_participant = 15)
      model1 <- lme4::lmer(mean_velocity ~ rir + (1 | id), data = data, REML = FALSE)
      model2 <- create_test_model(data)

      models <- list(simple = model1, full = model2)

      validator <- ModelValidator$new()
      result <- validator$loo_model_selection_stability(data, models, "BIC", "id")
      interpretation <- result$interpret()

      expect_type(interpretation, "character")
      expect_true(nchar(interpretation) > 0)
    })
  })

  describe("loo_prediction_error_by_participant", {
    skip_if_not_installed("lme4")

    it("calculates per-participant prediction error", {
      data <- create_test_data(n_participants = 10, obs_per_participant = 15)
      model <- create_test_model(data)
      validator <- ModelValidator$new()

      result <- validator$loo_prediction_error_by_participant(
        data, model, "id", "mean_velocity"
      )

      expect_s3_class(result, "PredictionErrorByParticipantResult")
      expect_length(result$participant_id, 10)
      expect_length(result$rmse, 10)
      expect_length(result$mae, 10)
      expect_true(is.numeric(result$overall_rmse))
    })

    it("stores residuals for anomaly detection", {
      data <- create_test_data(n_participants = 10, obs_per_participant = 15)
      model <- create_test_model(data)
      validator <- ModelValidator$new()

      result <- validator$loo_prediction_error_by_participant(
        data, model, "id", "mean_velocity"
      )

      expect_true(length(result$residuals) > 0)
      expect_equal(length(result$residuals), length(result$residual_participant_ids))
    })

    it("identifies hardest to predict participants", {
      data <- create_test_data(n_participants = 10, obs_per_participant = 15)
      model <- create_test_model(data)
      validator <- ModelValidator$new()

      result <- validator$loo_prediction_error_by_participant(
        data, model, "id", "mean_velocity"
      )

      hardest <- result$get_hardest_to_predict(3)

      expect_s3_class(hardest, "data.frame")
      expect_equal(nrow(hardest), 3)
      expect_true("rmse" %in% names(hardest))
      # Should be sorted by RMSE descending
      expect_true(hardest$rmse[1] >= hardest$rmse[2])
    })
  })

  describe("plot_coefficient_stability", {
    skip_if_not_installed("lme4")
    skip_if_not_installed("ggplot2")

    it("creates coefficient stability plot", {
      data <- create_test_data(n_participants = 10, obs_per_participant = 15)
      model <- create_test_model(data)
      validator <- ModelValidator$new()

      stability <- validator$loo_coefficient_stability(data, model, "rir", "id")
      plot <- validator$plot_coefficient_stability(stability)

      expect_s3_class(plot, "ggplot")
    })
  })
})

# =============================================================================
# INFLUENCE DIAGNOSTICS TESTS
# =============================================================================

describe("InfluenceDiagnostics", {

  describe("initialization", {
    it("creates instance", {
      diagnostics <- InfluenceDiagnostics$new()
      expect_s3_class(diagnostics, "InfluenceDiagnostics")
    })
  })

  describe("calculate_observation_influence", {
    skip_if_not_installed("lme4")

    it("calculates influence measures", {
      data <- create_test_data()
      model <- create_test_model(data)
      diagnostics <- InfluenceDiagnostics$new()

      result <- diagnostics$calculate_observation_influence(model)

      expect_s3_class(result, "InfluenceDiagnosticsResult")
      expect_equal(length(result$leverage), nrow(data))
      expect_equal(length(result$cooks_d), nrow(data))
      expect_gt(result$influential_threshold, 0)
    })

    it("identifies influential observations", {
      data <- create_test_data()
      model <- create_test_model(data)
      diagnostics <- InfluenceDiagnostics$new()

      result <- diagnostics$calculate_observation_influence(model)
      n_influential <- result$count_influential()

      expect_true(is.numeric(n_influential))
      expect_true(n_influential >= 0)
    })
  })

  describe("calculate_participant_influence", {
    skip_if_not_installed("lme4")

    it("calculates participant-level influence", {
      data <- create_test_data(n_participants = 10, obs_per_participant = 15)
      model <- create_test_model(data)
      diagnostics <- InfluenceDiagnostics$new()

      result <- diagnostics$calculate_participant_influence(model, data)

      expect_s3_class(result, "ParticipantInfluenceResult")
      expect_true(is.numeric(result$original_estimate))
      expect_true(is.character(result$most_influential_id))
      expect_equal(nrow(result$participant_effects), 10)
    })

    it("interprets influence magnitude", {
      data <- create_test_data(n_participants = 10, obs_per_participant = 15)
      model <- create_test_model(data)
      diagnostics <- InfluenceDiagnostics$new()

      result <- diagnostics$calculate_participant_influence(model, data)
      interpretation <- result$interpret()

      expect_type(interpretation, "character")
      expect_true(interpretation %in% c("minimal", "small", "moderate", "substantial"))
    })
  })
})

# =============================================================================
# MARGINAL EFFECTS TESTS
# =============================================================================

describe("MarginalEffectsAnalyzer", {

  describe("initialization", {
    it("creates instance", {
      analyzer <- MarginalEffectsAnalyzer$new()
      expect_s3_class(analyzer, "MarginalEffectsAnalyzer")
    })
  })

  describe("calculate_marginal_effect", {
    skip_if_not_installed("lme4")

    it("calculates marginal effect", {
      data <- create_test_data()
      model <- create_test_model(data)
      analyzer <- MarginalEffectsAnalyzer$new()

      result <- analyzer$calculate_marginal_effect(model, "rir", data)

      expect_s3_class(result, "MarginalEffectResult")
      expect_equal(result$focal_variable, "rir")
      expect_true(is.numeric(result$effect_estimate))
      expect_gt(result$effect_estimate, 0)  # RIR should have positive effect
      expect_lt(result$ci_lower, result$ci_upper)
    })

    it("generates prediction data", {
      data <- create_test_data()
      model <- create_test_model(data)
      analyzer <- MarginalEffectsAnalyzer$new()

      result <- analyzer$calculate_marginal_effect(model, "rir", data, n_points = 20)

      expect_equal(nrow(result$prediction_data), 20)
      expect_true("predicted" %in% names(result$prediction_data))
    })

    it("interprets effect", {
      data <- create_test_data()
      model <- create_test_model(data)
      analyzer <- MarginalEffectsAnalyzer$new()

      result <- analyzer$calculate_marginal_effect(model, "rir", data)
      interpretation <- result$interpret()

      expect_type(interpretation, "character")
      expect_true(grepl("rir", interpretation))
    })
  })

  describe("create_effect_table", {
    skip_if_not_installed("lme4")

    it("creates effect summary table", {
      data <- create_test_data()
      model <- create_test_model(data)
      analyzer <- MarginalEffectsAnalyzer$new()

      table <- analyzer$create_effect_table(model, data)

      expect_s3_class(table, "data.frame")
      expect_true("predictor" %in% names(table))
      expect_true("estimate" %in% names(table))
      expect_true("ci_lower" %in% names(table))
      expect_true("p_value" %in% names(table))
    })
  })
})

# =============================================================================
# POWER ANALYZER TESTS
# =============================================================================

describe("PowerAnalyzer", {

  describe("initialization", {
    it("creates instance", {
      analyzer <- PowerAnalyzer$new()
      expect_s3_class(analyzer, "PowerAnalyzer")
    })
  })

  describe("calculate_power", {
    skip_if_not_installed("lme4")

    it("calculates power analysis", {
      data <- create_test_data()
      model <- create_test_model(data)
      analyzer <- PowerAnalyzer$new()

      result <- analyzer$calculate_power(model, data)

      expect_s3_class(result, "PowerAnalysisResult")
      expect_equal(result$n_participants, 15)
      expect_gt(result$n_observations, 0)
      expect_gt(result$effective_n, 0)
      expect_lt(result$effective_n, result$n_observations)  # Should be < total due to clustering
    })

    it("calculates ICC", {
      data <- create_test_data()
      model <- create_test_model(data)
      analyzer <- PowerAnalyzer$new()

      result <- analyzer$calculate_power(model, data)

      expect_true(is.numeric(result$icc))
      expect_gt(result$icc, 0)
      expect_lt(result$icc, 1)
    })

    it("calculates minimum detectable effect", {
      data <- create_test_data()
      model <- create_test_model(data)
      analyzer <- PowerAnalyzer$new()

      result <- analyzer$calculate_power(model, data)

      expect_true(is.numeric(result$detectable_effect))
      expect_gt(result$detectable_effect, 0)
    })

    it("interprets power", {
      data <- create_test_data()
      model <- create_test_model(data)
      analyzer <- PowerAnalyzer$new()

      result <- analyzer$calculate_power(model, data)
      interpretation <- result$interpret()

      expect_type(interpretation, "character")
      expect_true(grepl("power", interpretation, ignore.case = TRUE))
    })
  })

  describe("calculate_required_n", {
    it("calculates required sample size", {
      analyzer <- PowerAnalyzer$new()

      result <- analyzer$calculate_required_n(
        target_effect = 0.02,
        residual_sd = 0.03,
        icc = 0.3,
        obs_per_participant = 20
      )

      expect_s3_class(result, "SampleSizeResult")
      expect_gt(result$required_n_participants, 0)
      expect_equal(result$target_power, 0.80)
    })
  })

  describe("create_power_summary", {
    skip_if_not_installed("lme4")

    it("creates comprehensive summary", {
      data <- create_test_data()
      model <- create_test_model(data)
      analyzer <- PowerAnalyzer$new()

      summary <- analyzer$create_power_summary(model, data)

      expect_type(summary, "list")
      expect_true("current_study" %in% names(summary))
      expect_true("future_study" %in% names(summary))
      expect_true("interpretation" %in% names(summary))
      expect_true("caveats" %in% names(summary))
    })
  })
})

# =============================================================================
# CONVENIENCE FUNCTION TESTS
# =============================================================================

describe("Convenience functions", {
  skip_if_not_installed("lme4")

  it("leave_one_out_cv works", {
    data <- create_test_data(n_participants = 8, obs_per_participant = 10)
    model <- create_test_model(data)

    result <- leave_one_out_cv(data, model)
    expect_s3_class(result, "CrossValidationResult")
  })

  it("calculate_calibration works", {
    data <- create_test_data()
    model <- create_test_model(data)

    result <- calculate_calibration(model, data)
    expect_s3_class(result, "CalibrationResult")
  })

  it("calculate_influence works", {
    data <- create_test_data()
    model <- create_test_model(data)

    result <- calculate_influence(model)
    expect_s3_class(result, "InfluenceDiagnosticsResult")
  })

  it("calculate_marginal_effect works", {
    data <- create_test_data()
    model <- create_test_model(data)

    result <- calculate_marginal_effect(model, "rir", data)
    expect_s3_class(result, "MarginalEffectResult")
  })

  it("calculate_power works", {
    data <- create_test_data()
    model <- create_test_model(data)

    result <- calculate_power(model, data)
    expect_s3_class(result, "PowerAnalysisResult")
  })
})
