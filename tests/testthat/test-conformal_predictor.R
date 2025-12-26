# tests/testthat/test-conformal_predictor.R
# Unit tests for ConformalPredictor

box::use(
  testthat[...],
  ../../R/calculators/conformal_predictor[
    ConformalPredictor,
    ConformalPredictionResult,
    CoverageComparisonResult
  ]
)

# ==============================================================================
# Test Fixtures
# ==============================================================================

create_calibration_data <- function(n = 100, seed = 42) {
  set.seed(seed)

  data <- data.frame(
    id = rep(paste0("P", 1:10), each = n / 10),
    rir = sample(0:7, n, replace = TRUE)
  )

  # Create velocity with relationship to RIR
  data$mean_velocity <- 0.25 + 0.04 * data$rir + rnorm(n, 0, 0.03)

  data
}

create_test_data_for_conformal <- function(n = 50, seed = 123) {
  set.seed(seed)

  data <- data.frame(
    id = rep(paste0("P", 1:10), length.out = n),
    rir = sample(0:7, n, replace = TRUE)
  )

  data$mean_velocity <- 0.25 + 0.04 * data$rir + rnorm(n, 0, 0.03)

  data
}

fit_simple_model <- function(data) {
  stats::lm(mean_velocity ~ rir, data = data)
}

# ==============================================================================
# ConformalPredictionResult Tests
# ==============================================================================

describe("ConformalPredictionResult", {

  describe("initialization", {
    it("stores prediction data", {
      predictions <- data.frame(
        prediction = c(0.25, 0.30, 0.35),
        lower = c(0.20, 0.25, 0.30),
        upper = c(0.30, 0.35, 0.40),
        interval_width = c(0.10, 0.10, 0.10)
      )

      result <- ConformalPredictionResult$new(
        predictions = predictions,
        coverage_target = 0.95,
        q_hat = 0.05,
        n_calibration = 100
      )

      expect_equal(result$coverage_target, 0.95)
      expect_equal(result$q_hat, 0.05)
      expect_equal(result$n_calibration, 100)
      expect_equal(nrow(result$predictions), 3)
    })
  })

  describe("to_list", {
    it("returns all fields as list", {
      predictions <- data.frame(
        prediction = 0.3,
        lower = 0.25,
        upper = 0.35,
        interval_width = 0.1
      )

      result <- ConformalPredictionResult$new(predictions, 0.95, 0.05, 50)
      as_list <- result$to_list()

      expect_true("predictions" %in% names(as_list))
      expect_true("coverage_target" %in% names(as_list))
      expect_true("q_hat" %in% names(as_list))
    })
  })
})

# ==============================================================================
# CoverageComparisonResult Tests
# ==============================================================================

describe("CoverageComparisonResult", {

  describe("initialization", {
    it("stores comparison data", {
      result <- CoverageComparisonResult$new(
        parametric_coverage = 0.90,
        conformal_coverage = 0.95,
        parametric_width = 0.08,
        conformal_width = 0.10,
        target_coverage = 0.95
      )

      expect_equal(result$parametric_coverage, 0.90)
      expect_equal(result$conformal_coverage, 0.95)
      expect_equal(result$target_coverage, 0.95)
    })
  })

  describe("to_list", {
    it("calculates coverage errors", {
      result <- CoverageComparisonResult$new(0.90, 0.94, 0.08, 0.10, 0.95)
      as_list <- result$to_list()

      expect_true("parametric_coverage_error" %in% names(as_list))
      expect_true("conformal_coverage_error" %in% names(as_list))
      expect_equal(as_list$parametric_coverage_error, 0.05)  # |0.90 - 0.95|
      expect_equal(as_list$conformal_coverage_error, 0.01)   # |0.94 - 0.95|
    })
  })
})

# ==============================================================================
# ConformalPredictor Tests
# ==============================================================================

describe("ConformalPredictor", {

  describe("initialization", {
    it("creates predictor with default alpha", {
      cp <- ConformalPredictor$new()
      expect_s3_class(cp, "ConformalPredictor")
    })

    it("creates predictor with custom alpha", {
      cp <- ConformalPredictor$new(alpha = 0.10)
      expect_s3_class(cp, "ConformalPredictor")
    })
  })

  describe("fit", {
    it("fits on calibration data", {
      cp <- ConformalPredictor$new(alpha = 0.05)
      cal_data <- create_calibration_data()
      model <- fit_simple_model(cal_data)

      cp$fit(model, cal_data, response_var = "mean_velocity")

      # Should be able to get q_hat after fitting
      q_hat <- cp$get_q_hat()
      expect_true(is.numeric(q_hat))
      expect_true(q_hat > 0)
    })

    it("calculates nonconformity scores", {
      cp <- ConformalPredictor$new()
      cal_data <- create_calibration_data()
      model <- fit_simple_model(cal_data)

      cp$fit(model, cal_data, response_var = "mean_velocity")

      scores <- cp$get_nonconformity_scores()
      expect_equal(length(scores), nrow(cal_data))
      expect_true(all(scores >= 0))  # Absolute residuals are non-negative
    })

    it("errors before fitting", {
      cp <- ConformalPredictor$new()

      expect_error(cp$get_q_hat(), "must be fitted")
      expect_error(cp$get_nonconformity_scores(), "must be fitted")
    })
  })

  describe("predict_interval", {
    it("generates prediction intervals", {
      cp <- ConformalPredictor$new(alpha = 0.05)
      cal_data <- create_calibration_data()
      test_data <- create_test_data_for_conformal()
      model <- fit_simple_model(cal_data)

      cp$fit(model, cal_data, "mean_velocity")
      result <- cp$predict_interval(test_data)

      expect_s3_class(result, "ConformalPredictionResult")
      expect_equal(nrow(result$predictions), nrow(test_data))
      expect_true("lower" %in% names(result$predictions))
      expect_true("upper" %in% names(result$predictions))
    })

    it("intervals have correct structure", {
      cp <- ConformalPredictor$new()
      cal_data <- create_calibration_data()
      test_data <- create_test_data_for_conformal()
      model <- fit_simple_model(cal_data)

      cp$fit(model, cal_data, "mean_velocity")
      result <- cp$predict_interval(test_data)

      predictions <- result$predictions
      expect_true(all(predictions$lower < predictions$prediction))
      expect_true(all(predictions$upper > predictions$prediction))
    })

    it("interval width equals 2 * q_hat", {
      cp <- ConformalPredictor$new()
      cal_data <- create_calibration_data()
      test_data <- create_test_data_for_conformal()
      model <- fit_simple_model(cal_data)

      cp$fit(model, cal_data, "mean_velocity")
      result <- cp$predict_interval(test_data)

      q_hat <- cp$get_q_hat()
      expected_width <- 2 * q_hat

      actual_widths <- result$predictions$interval_width
      expect_true(all(abs(actual_widths - expected_width) < 1e-10))
    })

    it("errors if not fitted", {
      cp <- ConformalPredictor$new()
      test_data <- create_test_data_for_conformal()

      expect_error(cp$predict_interval(test_data), "must be fitted")
    })
  })

  describe("calculate_coverage", {
    it("calculates empirical coverage on test set", {
      cp <- ConformalPredictor$new(alpha = 0.05)
      cal_data <- create_calibration_data()
      test_data <- create_test_data_for_conformal()
      model <- fit_simple_model(cal_data)

      cp$fit(model, cal_data, "mean_velocity")
      conformal_result <- cp$predict_interval(test_data)

      coverage <- cp$calculate_coverage(test_data, conformal_result)

      expect_true(is.numeric(coverage))
      expect_true(coverage >= 0 && coverage <= 1)
    })

    it("coverage is approximately at target level", {
      # Use larger samples for more stable coverage estimate
      set.seed(42)
      cp <- ConformalPredictor$new(alpha = 0.10)  # 90% coverage

      cal_data <- create_calibration_data(n = 200)
      test_data <- create_test_data_for_conformal(n = 100)
      model <- fit_simple_model(cal_data)

      cp$fit(model, cal_data, "mean_velocity")
      conformal_result <- cp$predict_interval(test_data)
      coverage <- cp$calculate_coverage(test_data, conformal_result)

      # Coverage should be roughly around target (allow ±15% for finite sample)
      expect_true(coverage >= 0.75)  # At least 75%
      expect_true(coverage <= 1.0)   # At most 100%
    })
  })

  describe("compare_with_parametric", {
    it("compares conformal and parametric intervals", {
      skip_if_not_installed("lme4")

      box::use(../../R/calculators/lmm_analyzer[LmmAnalyzer])

      cp <- ConformalPredictor$new(alpha = 0.05)
      cal_data <- create_calibration_data()
      test_data <- create_test_data_for_conformal()

      # Fit LMM
      analyzer <- LmmAnalyzer$new()
      model_result <- analyzer$fit(
        cal_data,
        mean_velocity ~ rir,
        ~1 | id,
        "test_model"
      )

      cp$fit(model_result$model, cal_data, "mean_velocity")
      comparison <- cp$compare_with_parametric(model_result, test_data)

      expect_s3_class(comparison, "CoverageComparisonResult")
      expect_true(comparison$parametric_coverage >= 0)
      expect_true(comparison$conformal_coverage >= 0)
    })

    it("returns width comparison", {
      skip_if_not_installed("lme4")

      box::use(../../R/calculators/lmm_analyzer[LmmAnalyzer])

      cp <- ConformalPredictor$new()
      cal_data <- create_calibration_data()
      test_data <- create_test_data_for_conformal()

      analyzer <- LmmAnalyzer$new()
      model_result <- analyzer$fit(cal_data, mean_velocity ~ rir, ~1 | id)

      cp$fit(model_result$model, cal_data, "mean_velocity")
      comparison <- cp$compare_with_parametric(model_result, test_data)

      expect_true(comparison$parametric_width > 0)
      expect_true(comparison$conformal_width > 0)
    })
  })

  describe("get_interval_width_distribution", {
    it("returns distribution statistics", {
      cp <- ConformalPredictor$new()
      cal_data <- create_calibration_data()
      test_data <- create_test_data_for_conformal()
      model <- fit_simple_model(cal_data)

      cp$fit(model, cal_data, "mean_velocity")
      result <- cp$predict_interval(test_data)

      dist <- cp$get_interval_width_distribution(result)

      expect_true("mean" %in% names(dist))
      expect_true("median" %in% names(dist))
      expect_true("sd" %in% names(dist))
      expect_true("min" %in% names(dist))
      expect_true("max" %in% names(dist))
    })

    it("all intervals have same width (split conformal)", {
      cp <- ConformalPredictor$new()
      cal_data <- create_calibration_data()
      test_data <- create_test_data_for_conformal()
      model <- fit_simple_model(cal_data)

      cp$fit(model, cal_data, "mean_velocity")
      result <- cp$predict_interval(test_data)
      dist <- cp$get_interval_width_distribution(result)

      # In split conformal with absolute residual score, all widths are equal
      expect_equal(dist$sd, 0)
      expect_equal(dist$min, dist$max)
    })
  })

  describe("works with LMM models", {
    it("handles lmerMod objects", {
      skip_if_not_installed("lme4")

      box::use(../../R/calculators/lmm_analyzer[LmmAnalyzer])

      cp <- ConformalPredictor$new()
      cal_data <- create_calibration_data()
      test_data <- create_test_data_for_conformal()

      analyzer <- LmmAnalyzer$new()
      model_result <- analyzer$fit(cal_data, mean_velocity ~ rir, ~1 | id)

      cp$fit(model_result$model, cal_data, "mean_velocity")
      result <- cp$predict_interval(test_data)

      expect_equal(nrow(result$predictions), nrow(test_data))
    })

    it("handles LmmModelResult wrapper", {
      skip_if_not_installed("lme4")

      box::use(../../R/calculators/lmm_analyzer[LmmAnalyzer])

      cp <- ConformalPredictor$new()
      cal_data <- create_calibration_data()
      test_data <- create_test_data_for_conformal()

      analyzer <- LmmAnalyzer$new()
      model_result <- analyzer$fit(cal_data, mean_velocity ~ rir, ~1 | id)

      # Pass the wrapper object
      cp$fit(model_result, cal_data, "mean_velocity")
      result <- cp$predict_interval(test_data)

      expect_equal(nrow(result$predictions), nrow(test_data))
    })
  })
})
