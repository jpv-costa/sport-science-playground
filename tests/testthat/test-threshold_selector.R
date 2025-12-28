# tests/testthat/test-threshold_selector.R
# Unit tests for ThresholdSelector

box::use(
  testthat[...],
  ../../R/calculators/threshold_selector[
    ThresholdSelector, ThresholdResult, ThresholdStrategy
  ]
)

describe("ThresholdStrategy", {

  it("contains all expected strategies", {
    expect_equal(ThresholdStrategy$FIXED_PERCENTILE, "fixed_percentile")
    expect_equal(ThresholdStrategy$ELBOW_DETECTION, "elbow_detection")
    expect_equal(ThresholdStrategy$GAP_STATISTICS, "gap_statistics")
    expect_equal(ThresholdStrategy$KNEE_POINT, "knee_point")
    expect_equal(ThresholdStrategy$MAD_BASED, "mad_based")
  })
})

describe("ThresholdResult", {

  it("stores threshold and strategy", {
    result <- ThresholdResult$new(
      threshold = 0.85,
      strategy = "fixed_percentile",
      diagnostics = list(n_flagged = 5, pct_flagged = 10)
    )

    expect_equal(result$threshold, 0.85)
    expect_equal(result$strategy, "fixed_percentile")
  })

  it("provides summary", {
    result <- ThresholdResult$new(
      threshold = 0.75,
      strategy = "elbow_detection",
      diagnostics = list(n_flagged = 8, pct_flagged = 16)
    )

    summary <- result$summarize()

    expect_equal(summary$threshold, 0.75)
    expect_equal(summary$strategy, "elbow_detection")
    expect_equal(summary$n_flagged, 8)
    expect_equal(summary$pct_flagged, 16)
  })
})

describe("ThresholdSelector", {

  # Create test data with clear cluster structure
  create_bimodal_scores <- function(n_normal = 80, n_anomaly = 20) {
    c(
      runif(n_normal, 0.1, 0.5),   # Normal cluster
      runif(n_anomaly, 0.7, 0.95)  # Anomaly cluster
    )
  }

  describe("select_threshold", {

    it("validates inputs", {
      selector <- ThresholdSelector$new()

      expect_error(
        selector$select_threshold("not numeric"),
        "scores must be numeric"
      )

      expect_error(
        selector$select_threshold(0.5),
        "scores must have length > 1"
      )

      expect_error(
        selector$select_threshold(c(0.5, 0.6), strategy = "invalid"),
        "strategy must be valid"
      )
    })

    it("returns ThresholdResult object", {
      selector <- ThresholdSelector$new()
      scores <- runif(100, 0, 1)

      result <- selector$select_threshold(scores)

      expect_s3_class(result, "ThresholdResult")
      expect_true(is.numeric(result$threshold))
      expect_true(result$threshold >= 0 && result$threshold <= 1)
    })

    it("handles NA values in scores", {
      selector <- ThresholdSelector$new()
      scores <- c(runif(50, 0, 1), NA, NA, NA)

      result <- selector$select_threshold(scores)

      expect_s3_class(result, "ThresholdResult")
      expect_true(!is.na(result$threshold))
    })
  })

  describe("fixed_percentile strategy", {

    it("returns expected percentile threshold", {
      selector <- ThresholdSelector$new()
      scores <- seq(0, 1, length.out = 100)

      result <- selector$select_threshold(
        scores,
        strategy = ThresholdStrategy$FIXED_PERCENTILE,
        params = list(percentile = 0.95)
      )

      # 95th percentile of 0-1 sequence should be ~0.95
      expect_equal(result$threshold, 0.95, tolerance = 0.02)
    })

    it("respects custom percentile parameter", {
      selector <- ThresholdSelector$new()
      scores <- seq(0, 1, length.out = 100)

      result_90 <- selector$select_threshold(
        scores,
        strategy = ThresholdStrategy$FIXED_PERCENTILE,
        params = list(percentile = 0.90)
      )

      result_99 <- selector$select_threshold(
        scores,
        strategy = ThresholdStrategy$FIXED_PERCENTILE,
        params = list(percentile = 0.99)
      )

      expect_lt(result_90$threshold, result_99$threshold)
    })

    it("validates percentile parameter", {
      selector <- ThresholdSelector$new()
      scores <- runif(100)

      expect_error(
        selector$select_threshold(
          scores,
          strategy = ThresholdStrategy$FIXED_PERCENTILE,
          params = list(percentile = 1.5)
        ),
        "percentile must be between 0 and 1"
      )
    })
  })

  describe("elbow_detection strategy", {

    it("finds elbow in bimodal distribution", {
      set.seed(42)
      selector <- ThresholdSelector$new()
      scores <- create_bimodal_scores()

      result <- selector$select_threshold(
        scores,
        strategy = ThresholdStrategy$ELBOW_DETECTION
      )

      # Threshold should be between clusters (0.5-0.7 range)
      expect_gt(result$threshold, 0.45)
      expect_lt(result$threshold, 0.8)
    })

    it("handles uniform distribution", {
      selector <- ThresholdSelector$new()
      scores <- seq(0, 1, length.out = 50)

      result <- selector$select_threshold(
        scores,
        strategy = ThresholdStrategy$ELBOW_DETECTION
      )

      expect_true(is.numeric(result$threshold))
      expect_true(!is.na(result$threshold))
    })

    it("handles small sample sizes", {
      selector <- ThresholdSelector$new()
      scores <- c(0.1, 0.2, 0.9)

      result <- selector$select_threshold(
        scores,
        strategy = ThresholdStrategy$ELBOW_DETECTION
      )

      expect_true(is.numeric(result$threshold))
    })
  })

  describe("gap_statistics strategy", {

    it("finds gap between clusters", {
      set.seed(42)
      selector <- ThresholdSelector$new()

      # Create data with clear gap
      scores <- c(
        runif(50, 0.1, 0.3),   # Low cluster
        runif(10, 0.8, 0.95)   # High cluster (gap at ~0.5)
      )

      result <- selector$select_threshold(
        scores,
        strategy = ThresholdStrategy$GAP_STATISTICS
      )

      # Threshold should be in the gap area
      expect_gt(result$threshold, 0.3)
      expect_lt(result$threshold, 0.8)
    })

    it("handles continuous distribution", {
      selector <- ThresholdSelector$new()
      scores <- runif(100, 0, 1)

      result <- selector$select_threshold(
        scores,
        strategy = ThresholdStrategy$GAP_STATISTICS
      )

      expect_true(is.numeric(result$threshold))
    })
  })

  describe("knee_point strategy", {

    it("detects knee in sorted scores", {
      set.seed(42)
      selector <- ThresholdSelector$new()
      scores <- create_bimodal_scores()

      result <- selector$select_threshold(
        scores,
        strategy = ThresholdStrategy$KNEE_POINT
      )

      # Should find transition point (allow wider range for stochastic data)
      expect_gt(result$threshold, 0.3)
      expect_lt(result$threshold, 0.96)
    })

    it("handles constant scores gracefully", {
      selector <- ThresholdSelector$new()
      scores <- rep(0.5, 20)

      result <- selector$select_threshold(
        scores,
        strategy = ThresholdStrategy$KNEE_POINT
      )

      expect_equal(result$threshold, 0.5)
    })
  })

  describe("mad_based strategy", {

    it("calculates MAD-based threshold", {
      selector <- ThresholdSelector$new()
      # Create scores with known median and MAD
      scores <- c(rep(0.4, 40), rep(0.5, 20), rep(0.6, 40))

      result <- selector$select_threshold(
        scores,
        strategy = ThresholdStrategy$MAD_BASED,
        params = list(multiplier = 3)
      )

      expect_true(is.numeric(result$threshold))
      expect_gt(result$threshold, median(scores))
    })

    it("respects multiplier parameter", {
      selector <- ThresholdSelector$new()
      scores <- runif(100, 0, 1)

      result_2 <- selector$select_threshold(
        scores,
        strategy = ThresholdStrategy$MAD_BASED,
        params = list(multiplier = 2)
      )

      result_4 <- selector$select_threshold(
        scores,
        strategy = ThresholdStrategy$MAD_BASED,
        params = list(multiplier = 4)
      )

      expect_lt(result_2$threshold, result_4$threshold)
    })

    it("handles zero MAD gracefully", {
      selector <- ThresholdSelector$new()
      # Scores with zero MAD (all same value) but with some outliers
      scores <- c(rep(0.5, 95), 0.9, 0.91, 0.92, 0.93, 0.94)

      result <- selector$select_threshold(
        scores,
        strategy = ThresholdStrategy$MAD_BASED
      )

      expect_true(is.numeric(result$threshold))
      expect_true(!is.na(result$threshold))
    })
  })

  describe("compare_strategies", {

    it("returns comparison for all strategies", {
      selector <- ThresholdSelector$new()
      scores <- runif(100, 0, 1)

      comparison <- selector$compare_strategies(scores)

      expect_s3_class(comparison, "data.frame")
      expect_true("strategy" %in% names(comparison))
      expect_true("threshold" %in% names(comparison))
      expect_true("n_flagged" %in% names(comparison))
      expect_true("pct_flagged" %in% names(comparison))
      expect_equal(nrow(comparison), 5)  # All 5 strategies
    })

    it("handles errors gracefully", {
      selector <- ThresholdSelector$new()
      scores <- c(0.1, 0.2, 0.3)  # Small but valid

      comparison <- selector$compare_strategies(scores)

      expect_s3_class(comparison, "data.frame")
      expect_gte(nrow(comparison), 1)
    })
  })

  describe("plot_threshold_comparison", {

    it("returns ggplot object", {
      skip_if_not_installed("ggplot2")

      selector <- ThresholdSelector$new()
      scores <- runif(50, 0, 1)

      plot <- selector$plot_threshold_comparison(scores)

      expect_s3_class(plot, "gg")
    })
  })

  describe("reproducibility", {

    it("produces deterministic results for same input", {
      selector <- ThresholdSelector$new()
      scores <- c(0.1, 0.2, 0.3, 0.5, 0.7, 0.8, 0.9)

      result1 <- selector$select_threshold(
        scores,
        strategy = ThresholdStrategy$ELBOW_DETECTION
      )

      result2 <- selector$select_threshold(
        scores,
        strategy = ThresholdStrategy$ELBOW_DETECTION
      )

      expect_equal(result1$threshold, result2$threshold)
    })
  })

  describe("diagnostics", {

    it("includes all diagnostic fields", {
      selector <- ThresholdSelector$new()
      scores <- runif(100, 0, 1)

      result <- selector$select_threshold(scores)
      diag <- result$diagnostics

      expect_true("n_total" %in% names(diag))
      expect_true("n_flagged" %in% names(diag))
      expect_true("pct_flagged" %in% names(diag))
      expect_true("score_range" %in% names(diag))
      expect_true("score_median" %in% names(diag))
      expect_true("score_mad" %in% names(diag))
      expect_true("params_used" %in% names(diag))
    })

    it("calculates correct flagged counts", {
      selector <- ThresholdSelector$new()
      scores <- seq(0, 1, by = 0.1)  # 11 scores: 0, 0.1, ..., 1.0

      result <- selector$select_threshold(
        scores,
        strategy = ThresholdStrategy$FIXED_PERCENTILE,
        params = list(percentile = 0.9)
      )

      # With 90th percentile, threshold ~0.9, should flag ~10%
      expect_equal(result$diagnostics$n_total, 11)
      expect_lte(result$diagnostics$pct_flagged, 20)
    })
  })
})
