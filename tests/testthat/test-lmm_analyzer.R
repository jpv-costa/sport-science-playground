# tests/testthat/test-lmm_analyzer.R
# Unit tests for LmmAnalyzer

box::use(
  testthat[...],
  ../../R/calculators/lmm_analyzer[
    LmmAnalyzer,
    LmmModelResult,
    LmmModelComparison,
    LmmDiagnostics
  ]
)

# ==============================================================================
# Test Fixtures
# ==============================================================================

create_nested_data <- function(n_participants = 10, n_obs_per = 20, seed = 42) {
  set.seed(seed)

  # Create nested data structure similar to deadlift study
  data <- expand.grid(
    participant_num = 1:n_participants,
    obs = 1:n_obs_per
  )

  data$id <- paste0("P", data$participant_num)

  # Random intercepts and slopes per participant
  participant_effects <- data.frame(
    id = paste0("P", 1:n_participants),
    intercept_effect = rnorm(n_participants, 0, 0.05),
    slope_effect = rnorm(n_participants, 0, 0.01)
  )

  data <- merge(data, participant_effects, by = "id")

  # Generate RIR values (0-7)
  data$rir <- sample(0:7, nrow(data), replace = TRUE)

  # Generate load percentage (80% or 90%)
  data$load_percentage <- sample(c("80%", "90%"), nrow(data), replace = TRUE)
  data$load_numeric <- ifelse(data$load_percentage == "80%", 80, 90)

  # Generate day (Day 1 or Day 2)
  data$day <- sample(c("Day 1", "Day 2"), nrow(data), replace = TRUE)

  # Generate velocity: higher RIR = higher velocity
  # Base relationship: velocity = 0.2 + 0.04 * rir - 0.002 * load_numeric
  # Add participant random effects and noise
  data$mean_velocity <- 0.2 +
    0.04 * data$rir -
    0.002 * data$load_numeric +
    data$intercept_effect +
    data$slope_effect * data$rir +
    rnorm(nrow(data), 0, 0.02)

  # Clean up
  data <- data[, c("id", "rir", "load_percentage", "day", "mean_velocity")]
  rownames(data) <- NULL

  data
}

# ==============================================================================
# LmmModelResult Tests
# ==============================================================================

describe("LmmModelResult", {

  describe("initialization", {
    it("stores model metadata correctly", {
      # Create a mock model (NULL is fine for testing the result class)
      result <- LmmModelResult$new(
        model = NULL,
        model_name = "test_model",
        formula_string = "mean_velocity ~ rir + (1 | id)",
        aic = 100.5,
        bic = 110.2,
        log_likelihood = -45.3,
        r2_marginal = 0.25,
        r2_conditional = 0.65,
        n_observations = 200,
        n_groups = 10
      )

      expect_equal(result$model_name, "test_model")
      expect_equal(result$aic, 100.5)
      expect_equal(result$bic, 110.2)
      expect_equal(result$r2_marginal, 0.25)
      expect_equal(result$r2_conditional, 0.65)
      expect_equal(result$n_observations, 200)
      expect_equal(result$n_groups, 10)
    })
  })

  describe("to_list", {
    it("returns all fields as list", {
      result <- LmmModelResult$new(
        model = NULL,
        model_name = "test",
        formula_string = "y ~ x",
        aic = 50,
        bic = 55,
        log_likelihood = -20,
        r2_marginal = 0.3,
        r2_conditional = 0.5,
        n_observations = 100,
        n_groups = 5
      )

      as_list <- result$to_list()

      expect_true("model_name" %in% names(as_list))
      expect_true("aic" %in% names(as_list))
      expect_true("r2_marginal" %in% names(as_list))
    })
  })
})

# ==============================================================================
# LmmAnalyzer Tests
# ==============================================================================

describe("LmmAnalyzer", {

  describe("initialization", {
    it("creates analyzer without error", {
      analyzer <- LmmAnalyzer$new()
      expect_s3_class(analyzer, "LmmAnalyzer")
    })
  })

  describe("fit", {
    it("fits LMM with random intercepts", {
      skip_if_not_installed("lme4")

      analyzer <- LmmAnalyzer$new()
      data <- create_nested_data()

      result <- analyzer$fit(
        data = data,
        formula = mean_velocity ~ rir,
        random_formula = ~1 | id,
        model_name = "random_intercept"
      )

      expect_s3_class(result, "LmmModelResult")
      expect_equal(result$model_name, "random_intercept")
      expect_true(result$n_observations > 0)
      expect_true(result$n_groups > 0)
      expect_true(!is.na(result$aic))
    })

    it("fits LMM with random intercepts and slopes", {
      skip_if_not_installed("lme4")

      analyzer <- LmmAnalyzer$new()
      data <- create_nested_data()

      result <- analyzer$fit(
        data = data,
        formula = mean_velocity ~ rir,
        random_formula = ~1 + rir | id,
        model_name = "random_slope"
      )

      expect_s3_class(result, "LmmModelResult")
      expect_true(!is.na(result$r2_marginal))
      expect_true(!is.na(result$r2_conditional))
    })

    it("handles multiple fixed effects", {
      skip_if_not_installed("lme4")

      analyzer <- LmmAnalyzer$new()
      data <- create_nested_data()

      result <- analyzer$fit(
        data = data,
        formula = mean_velocity ~ rir + load_percentage + day,
        random_formula = ~1 | id
      )

      expect_s3_class(result, "LmmModelResult")
      expect_true(grepl("load_percentage", result$formula_string))
    })
  })

  describe("fit_with_interactions", {
    it("adds interaction terms to formula", {
      skip_if_not_installed("lme4")

      analyzer <- LmmAnalyzer$new()
      data <- create_nested_data()

      result <- analyzer$fit_with_interactions(
        data = data,
        base_formula = mean_velocity ~ rir + load_percentage,
        interactions = c("rir:load_percentage"),
        random_formula = ~1 | id
      )

      expect_s3_class(result, "LmmModelResult")
      expect_true(grepl("rir:load_percentage", result$formula_string))
    })
  })

  describe("compare_models", {
    it("compares multiple models and returns comparison", {
      skip_if_not_installed("lme4")

      analyzer <- LmmAnalyzer$new()
      data <- create_nested_data()

      model1 <- analyzer$fit(data, mean_velocity ~ rir, ~1 | id, "base")
      model2 <- analyzer$fit(data, mean_velocity ~ rir + load_percentage, ~1 | id, "with_load")

      comparison <- analyzer$compare_models(
        list(base = model1, with_load = model2),
        criterion = "AIC"
      )

      expect_s3_class(comparison, "LmmModelComparison")
      expect_true(nrow(comparison$comparison_table) == 2)
      expect_true("delta_AIC" %in% names(comparison$comparison_table))
      expect_true("bayes_factor_approx" %in% names(comparison$comparison_table))
      expect_true(comparison$best_model_name %in% c("base", "with_load"))
    })

    it("calculates delta values correctly", {
      skip_if_not_installed("lme4")

      analyzer <- LmmAnalyzer$new()
      data <- create_nested_data()

      model1 <- analyzer$fit(data, mean_velocity ~ rir, ~1 | id, "m1")
      model2 <- analyzer$fit(data, mean_velocity ~ rir + load_percentage, ~1 | id, "m2")

      comparison <- analyzer$compare_models(list(m1 = model1, m2 = model2))

      # Best model should have delta_AIC = 0
      best_row <- comparison$comparison_table[
        comparison$comparison_table$model == comparison$best_model_name,
      ]
      expect_equal(best_row$delta_AIC, 0)
    })
  })

  describe("test_variable_importance", {
    it("tests if variable improves model fit", {
      skip_if_not_installed("lme4")

      analyzer <- LmmAnalyzer$new()
      data <- create_nested_data()

      result <- analyzer$test_variable_importance(
        data = data,
        base_formula = mean_velocity ~ rir,
        full_formula = mean_velocity ~ rir + load_percentage,
        random_formula = ~1 | id
      )

      expect_true("lrt_p_value" %in% names(result))
      expect_true("delta_aic" %in% names(result))
      expect_true("bayes_factor" %in% names(result))
      expect_true("bf_interpretation" %in% names(result))
      expect_true("variable_significant" %in% names(result))
      expect_true(is.logical(result$variable_significant))
    })
  })

  describe("test_assumptions", {
    it("returns diagnostics object", {
      skip_if_not_installed("lme4")

      analyzer <- LmmAnalyzer$new()
      data <- create_nested_data()

      model_result <- analyzer$fit(data, mean_velocity ~ rir, ~1 | id)
      diagnostics <- analyzer$test_assumptions(model_result)

      expect_s3_class(diagnostics, "LmmDiagnostics")
      expect_true(length(diagnostics$residuals) > 0)
      expect_true(length(diagnostics$fitted) > 0)
    })

    it("includes normality test results", {
      skip_if_not_installed("lme4")

      analyzer <- LmmAnalyzer$new()
      data <- create_nested_data()

      model_result <- analyzer$fit(data, mean_velocity ~ rir, ~1 | id)
      diagnostics <- analyzer$test_assumptions(model_result)

      normality <- diagnostics$normality_test
      expect_true("p_value" %in% names(normality))
      expect_true("is_normal" %in% names(normality))
      expect_true(is.logical(normality$is_normal))
    })

    it("includes homoscedasticity assessment", {
      skip_if_not_installed("lme4")

      analyzer <- LmmAnalyzer$new()
      data <- create_nested_data()

      model_result <- analyzer$fit(data, mean_velocity ~ rir, ~1 | id)
      diagnostics <- analyzer$test_assumptions(model_result)

      homoscedasticity <- diagnostics$homoscedasticity_test
      expect_true("correlation" %in% names(homoscedasticity))
      expect_true("is_homoscedastic" %in% names(homoscedasticity))
    })
  })

  describe("qq_plot_data", {
    it("returns data for Q-Q plot", {
      skip_if_not_installed("lme4")

      analyzer <- LmmAnalyzer$new()
      data <- create_nested_data()

      model_result <- analyzer$fit(data, mean_velocity ~ rir, ~1 | id)
      diagnostics <- analyzer$test_assumptions(model_result)
      qq_data <- analyzer$qq_plot_data(diagnostics)

      expect_true("theoretical" %in% names(qq_data))
      expect_true("sample" %in% names(qq_data))
      expect_equal(nrow(qq_data), length(diagnostics$residuals))
    })
  })

  describe("residual_plot_data", {
    it("returns data for residual plot", {
      skip_if_not_installed("lme4")

      analyzer <- LmmAnalyzer$new()
      data <- create_nested_data()

      model_result <- analyzer$fit(data, mean_velocity ~ rir, ~1 | id)
      diagnostics <- analyzer$test_assumptions(model_result)
      resid_data <- analyzer$residual_plot_data(diagnostics)

      expect_true("fitted" %in% names(resid_data))
      expect_true("residuals" %in% names(resid_data))
    })
  })

  describe("get_fixed_effects", {
    it("returns fixed effects table", {
      skip_if_not_installed("lme4")
      skip_if_not_installed("lmerTest")

      analyzer <- LmmAnalyzer$new()
      data <- create_nested_data()

      model_result <- analyzer$fit(data, mean_velocity ~ rir + load_percentage, ~1 | id)
      fixed_effects <- analyzer$get_fixed_effects(model_result)

      expect_true("term" %in% names(fixed_effects))
      expect_true("estimate" %in% names(fixed_effects))
      expect_true("std_error" %in% names(fixed_effects))
      expect_true(nrow(fixed_effects) >= 2)  # Intercept + predictors
    })
  })

  describe("get_random_effects", {
    it("returns random effects", {
      skip_if_not_installed("lme4")

      analyzer <- LmmAnalyzer$new()
      data <- create_nested_data()

      model_result <- analyzer$fit(data, mean_velocity ~ rir, ~1 | id)
      random_effects <- analyzer$get_random_effects(model_result)

      expect_true("estimates" %in% names(random_effects))
      expect_true("variances" %in% names(random_effects))
    })
  })

  describe("predict_values", {
    it("predicts with random effects", {
      skip_if_not_installed("lme4")

      analyzer <- LmmAnalyzer$new()
      data <- create_nested_data()

      model_result <- analyzer$fit(data, mean_velocity ~ rir, ~1 | id)
      predictions <- analyzer$predict_values(model_result, data, include_random = TRUE)

      expect_equal(length(predictions), nrow(data))
      expect_true(all(is.finite(predictions)))
    })

    it("predicts without random effects (population average)", {
      skip_if_not_installed("lme4")

      analyzer <- LmmAnalyzer$new()
      data <- create_nested_data()

      model_result <- analyzer$fit(data, mean_velocity ~ rir, ~1 | id)
      pred_with <- analyzer$predict_values(model_result, data, include_random = TRUE)
      pred_without <- analyzer$predict_values(model_result, data, include_random = FALSE)

      # Predictions should differ (random effects matter)
      expect_false(all(pred_with == pred_without))
    })
  })
})
