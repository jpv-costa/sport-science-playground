# tests/testthat/test-meta_regression.R
# Unit tests for MetaRegressionModel

box::use(
  testthat[...],
  ../../R/models/meta_regression[
    MetaRegressionModel,
    MetaRegressionResult
  ]
)

# ==============================================================================
# Test Fixtures
# ==============================================================================

create_mock_model <- function() {
  # Create a minimal mock object that mimics rma.mv structure
  mock_model <- list(
    b = matrix(c(0.5, -0.02, 0.1), ncol = 1, dimnames = list(
      c("intrcpt", "avg.rir", "weeks"), "estimate"
    )),
    se = c(0.15, 0.01, 0.03),
    zval = c(3.33, -2.0, 3.33),
    pval = c(0.001, 0.045, 0.001),
    ci.lb = c(0.2, -0.04, 0.04),
    ci.ub = c(0.8, -0.001, 0.16),
    sigma2 = c(0.1, 0.05, 0.02),
    fit.stats = list(
      ML = c(-50, 100, 110, 112, 111)
    )
  )
  class(mock_model) <- "rma.mv"
  mock_model
}

create_synthetic_meta_data <- function(n_studies = 5, effects_per_study = 4) {
  set.seed(42)
  n <- n_studies * effects_per_study

  data.frame(
    study = rep(1:n_studies, each = effects_per_study),
    group = rep(1:2, times = n / 2),
    obs = 1:n,
    yi = rnorm(n, mean = 0.5, sd = 0.3),
    vi = runif(n, 0.01, 0.1),
    avg.rir = runif(n, 0, 4),
    load.set = runif(n, 60, 90),
    set.rep.equated = sample(0:1, n, replace = TRUE),
    weeks = sample(4:16, n, replace = TRUE),
    train.status = sample(c("trained", "untrained"), n, replace = TRUE)
  )
}

# ==============================================================================
# MetaRegressionResult Tests
# ==============================================================================

describe("MetaRegressionResult", {

  describe("initialization", {
    it("stores model and metadata", {
      mock_model <- create_mock_model()

      result <- MetaRegressionResult$new(
        model = mock_model,
        outcome = "Strength",
        n_effects = 50,
        n_studies = 10
      )

      expect_equal(result$outcome, "Strength")
      expect_equal(result$n_effects, 50)
      expect_equal(result$n_studies, 10)
    })
  })

  describe("get_coefficient", {
    it("extracts coefficient by exact name", {
      mock_model <- create_mock_model()
      result <- MetaRegressionResult$new(mock_model, "Test", 20, 5)

      coef <- result$get_coefficient("avg.rir")

      expect_equal(coef$estimate, -0.02)
      expect_equal(coef$standard_error, 0.01)
      expect_equal(coef$p_value, 0.045)
      expect_equal(coef$ci_lower, -0.04)
      expect_equal(coef$ci_upper, -0.001)
    })

    it("uses partial matching when exact fails", {
      mock_model <- create_mock_model()
      result <- MetaRegressionResult$new(mock_model, "Test", 20, 5)

      coef <- result$get_coefficient("avg")

      expect_equal(coef$estimate, -0.02)  # Matches avg.rir
    })

    it("returns NAs for nonexistent predictor", {
      mock_model <- create_mock_model()
      result <- MetaRegressionResult$new(mock_model, "Test", 20, 5)

      coef <- result$get_coefficient("nonexistent")

      expect_true(is.na(coef$estimate))
      expect_true(is.na(coef$p_value))
    })
  })

  describe("is_rir_significant", {
    it("returns TRUE when p < alpha", {
      mock_model <- create_mock_model()
      result <- MetaRegressionResult$new(mock_model, "Test", 20, 5)

      expect_true(result$is_rir_significant(alpha = 0.05))
    })

    it("returns FALSE when p > alpha", {
      mock_model <- create_mock_model()
      result <- MetaRegressionResult$new(mock_model, "Test", 20, 5)

      expect_false(result$is_rir_significant(alpha = 0.01))
    })
  })

  describe("get_coefficients_table", {
    it("returns all coefficients as data frame", {
      mock_model <- create_mock_model()
      result <- MetaRegressionResult$new(mock_model, "Test", 20, 5)

      coef_table <- result$get_coefficients_table()

      expect_s3_class(coef_table, "data.frame")
      expect_equal(nrow(coef_table), 3)
      expect_true("predictor" %in% names(coef_table))
      expect_true("estimate" %in% names(coef_table))
      expect_true("p_value" %in% names(coef_table))
    })
  })

  describe("get_variance_components", {
    it("extracts sigma2 values", {
      mock_model <- create_mock_model()
      result <- MetaRegressionResult$new(mock_model, "Test", 20, 5)

      var_comp <- result$get_variance_components()

      expect_equal(var_comp$sigma2_study, 0.1)
      expect_equal(var_comp$sigma2_group, 0.05)
      expect_equal(var_comp$sigma2_obs, 0.02)
      expect_equal(var_comp$total_heterogeneity, 0.17)
    })
  })

  describe("summarize", {
    it("creates comprehensive summary", {
      mock_model <- create_mock_model()
      result <- MetaRegressionResult$new(mock_model, "Hypertrophy", 100, 20)

      summary <- result$summarize()

      expect_equal(summary$outcome, "Hypertrophy")
      expect_equal(summary$n_effects, 100)
      expect_equal(summary$n_studies, 20)
      expect_equal(summary$rir_effect, -0.02)
      expect_equal(summary$rir_p_value, 0.045)
      expect_true(summary$rir_significant)
    })
  })
})

# ==============================================================================
# MetaRegressionModel Tests
# ==============================================================================

describe("MetaRegressionModel", {

  describe("initialization", {
    it("accepts custom moderators", {
      model <- MetaRegressionModel$new(
        moderators = c("avg.rir", "weeks")
      )

      expect_s3_class(model, "MetaRegressionModel")
    })

    it("uses default moderators when not specified", {
      model <- MetaRegressionModel$new()

      expect_s3_class(model, "MetaRegressionModel")
    })
  })

  describe("fit", {
    it("validates required columns exist", {
      model <- MetaRegressionModel$new(
        moderators = c("avg.rir")
      )

      incomplete_data <- data.frame(
        yi = c(0.5, 0.6),
        vi = c(0.1, 0.1)
        # Missing study, group, obs, avg.rir
      )

      expect_error(
        model$fit(incomplete_data, "Test"),
        "Missing required columns"
      )
    })

    it("validates no NA in effect sizes", {
      model <- MetaRegressionModel$new(
        moderators = c("avg.rir")
      )

      data_with_na <- data.frame(
        study = 1:2,
        group = 1:2,
        obs = 1:2,
        yi = c(0.5, NA),
        vi = c(0.1, 0.1),
        avg.rir = c(2, 3)
      )

      expect_error(
        model$fit(data_with_na, "Test"),
        "cannot contain NA"
      )
    })

    it("fits model with valid data", {
      skip_if_not_installed("metafor")

      model <- MetaRegressionModel$new(
        moderators = c("avg.rir", "weeks"),
        random_effects = list(~1|study, ~1|group, ~1|obs)
      )

      data <- create_synthetic_meta_data(n_studies = 5, effects_per_study = 4)

      result <- model$fit(data, outcome = "Strength")

      expect_s3_class(result, "MetaRegressionResult")
      expect_equal(result$outcome, "Strength")
      expect_equal(result$n_studies, 5)
      expect_equal(result$n_effects, 20)
    })
  })

  describe("fit_both_outcomes", {
    it("returns results for both outcomes", {
      skip_if_not_installed("metafor")

      model <- MetaRegressionModel$new(
        moderators = c("avg.rir"),
        random_effects = list(~1|study, ~1|group, ~1|obs)
      )

      strength_data <- create_synthetic_meta_data(n_studies = 3, effects_per_study = 4)
      hypertrophy_data <- create_synthetic_meta_data(n_studies = 4, effects_per_study = 3)

      results <- model$fit_both_outcomes(strength_data, hypertrophy_data)

      expect_true("strength" %in% names(results))
      expect_true("hypertrophy" %in% names(results))
      expect_equal(results$strength$outcome, "Strength")
      expect_equal(results$hypertrophy$outcome, "Hypertrophy")
    })
  })
})
