# tests/testthat/test-correlation_imputer.R
# Unit tests for CorrelationImputer

box::use(
  testthat[...],
  ../../R/calculators/correlation_imputer[
    CorrelationImputer,
    CorrelationImputationResult
  ]
)

describe("CorrelationImputationResult", {

  describe("initialization", {
    it("stores all imputation metadata", {
      result <- CorrelationImputationResult$new(
        imputed_correlation = 0.75,
        fisher_z = 0.973,
        standard_error = 0.05,
        n_studies = 10,
        method = "meta-analytic"
      )

      expect_equal(result$imputed_correlation, 0.75)
      expect_equal(result$fisher_z, 0.973)
      expect_equal(result$standard_error, 0.05)
      expect_equal(result$n_studies, 10)
      expect_equal(result$method, "meta-analytic")
    })
  })

  describe("to_list", {
    it("returns all fields as list", {
      result <- CorrelationImputationResult$new(
        imputed_correlation = 0.8,
        fisher_z = 1.099,
        standard_error = 0.03,
        n_studies = 15,
        method = "robust-meta-analysis"
      )
      as_list <- result$to_list()

      expect_equal(as_list$imputed_correlation, 0.8)
      expect_equal(as_list$fisher_z, 1.099)
      expect_equal(as_list$standard_error, 0.03)
      expect_equal(as_list$n_studies, 15)
      expect_equal(as_list$method, "robust-meta-analysis")
    })
  })
})

describe("CorrelationImputer", {

  describe("calculate_correlation", {
    it("computes correlation from SDs using standard formula", {
      imputer <- CorrelationImputer$new()

      # Known case: ri = (sd_pre^2 + sd_post^2 - sd_delta^2) / (2 * sd_pre * sd_post)
      # If pre_sd = 10, post_sd = 10, delta_sd = 5
      # ri = (100 + 100 - 25) / (2 * 10 * 10) = 175 / 200 = 0.875
      result <- imputer$calculate_correlation(pre_sd = 10, post_sd = 10, delta_sd = 5)

      expect_equal(result, 0.875, tolerance = 0.001)
    })

    it("returns NA for zero denominator", {
      imputer <- CorrelationImputer$new()

      result <- imputer$calculate_correlation(pre_sd = 0, post_sd = 10, delta_sd = 5)

      expect_true(is.na(result))
    })

    it("returns NA when calculated r exceeds bounds", {
      imputer <- CorrelationImputer$new()

      # If delta_sd is very small, r > 1, which is invalid
      # ri = (100 + 100 - 0) / (2 * 10 * 10) = 200 / 200 = 1.0 (edge case, valid)
      result_valid <- imputer$calculate_correlation(pre_sd = 10, post_sd = 10, delta_sd = 0)
      expect_equal(result_valid, 1.0, tolerance = 0.001)

      # Large delta_sd gives r < -1
      # ri = (100 + 100 - 900) / (2 * 10 * 10) = -700 / 200 = -3.5 (invalid)
      result_invalid <- imputer$calculate_correlation(pre_sd = 10, post_sd = 10, delta_sd = 30)
      expect_true(is.na(result_invalid))
    })

    it("handles realistic strength training data", {
      imputer <- CorrelationImputer$new()

      # Typical values from strength studies: high correlation expected
      # Pre: mean=100kg, SD=15kg; Post: mean=110kg, SD=16kg; Delta: SD=5kg
      result <- imputer$calculate_correlation(pre_sd = 15, post_sd = 16, delta_sd = 5)

      # ri = (225 + 256 - 25) / (2 * 15 * 16) = 456 / 480 = 0.95
      expect_equal(result, 0.95, tolerance = 0.001)
    })
  })

  describe("apply_imputation", {
    it("fills NA values with imputed correlation", {
      imputer <- CorrelationImputer$new()

      data <- data.frame(
        study = c("A", "A", "B", "B"),
        ri = c(0.8, NA, 0.7, NA)
      )

      result <- imputer$apply_imputation(data, imputed_value = 0.75)

      expect_equal(result$ri[1], 0.8)   # Original preserved
      expect_equal(result$ri[2], 0.75)  # Imputed
      expect_equal(result$ri[3], 0.7)   # Original preserved
      expect_equal(result$ri[4], 0.75)  # Imputed
    })

    it("leaves all values unchanged when no NAs present", {
      imputer <- CorrelationImputer$new()

      data <- data.frame(
        study = c("A", "B"),
        ri = c(0.8, 0.9)
      )

      result <- imputer$apply_imputation(data, imputed_value = 0.5)

      expect_equal(result$ri, c(0.8, 0.9))
    })
  })

  describe("impute_via_meta_analysis", {
    it("performs meta-analytic imputation on valid data", {
      skip_if_not_installed("metafor")
      skip_if_not_installed("psych")

      imputer <- CorrelationImputer$new()

      # Create synthetic data with known correlations
      set.seed(42)
      data <- data.frame(
        study = rep(1:5, each = 4),
        group = rep(1:2, times = 10),
        obs = 1:20,
        ri = c(0.7, 0.75, 0.8, 0.72,   # Study 1
               0.65, 0.70, 0.68, 0.73,  # Study 2
               0.78, 0.82, 0.79, 0.81,  # Study 3
               0.71, 0.74, 0.69, 0.76,  # Study 4
               0.80, 0.77, 0.83, 0.79), # Study 5
        n = rep(c(20, 25, 30, 22), 5)
      )

      result <- imputer$impute_via_meta_analysis(data, use_robust = FALSE)

      expect_s3_class(result, "CorrelationImputationResult")
      expect_true(result$imputed_correlation > 0.5)
      expect_true(result$imputed_correlation < 1.0)
      expect_equal(result$n_studies, 5)
      expect_equal(result$method, "meta-analysis")
    })

    it("uses robust variance estimation when requested", {
      skip_if_not_installed("metafor")
      skip_if_not_installed("psych")

      imputer <- CorrelationImputer$new()

      # Use stable correlations (less variance) for convergence
      data <- data.frame(
        study = rep(1:5, each = 4),
        group = rep(1:2, times = 10),
        obs = 1:20,
        ri = c(0.70, 0.72, 0.71, 0.73,   # Study 1: tight cluster
               0.68, 0.70, 0.69, 0.71,   # Study 2: tight cluster
               0.75, 0.77, 0.76, 0.74,   # Study 3: tight cluster
               0.69, 0.71, 0.70, 0.72,   # Study 4: tight cluster
               0.73, 0.75, 0.74, 0.76),  # Study 5: tight cluster
        n = rep(c(30, 35, 40, 32), 5)    # Larger sample sizes
      )

      result <- imputer$impute_via_meta_analysis(data, use_robust = TRUE)

      expect_equal(result$method, "robust-meta-analysis")
    })

    it("errors when insufficient valid correlations", {
      imputer <- CorrelationImputer$new()

      data <- data.frame(
        study = 1,
        group = 1,
        obs = 1,
        ri = 0.7,
        n = 20
      )

      expect_error(
        imputer$impute_via_meta_analysis(data),
        "At least 2 valid correlations"
      )
    })
  })
})
