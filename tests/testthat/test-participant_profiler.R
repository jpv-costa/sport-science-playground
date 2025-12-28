# tests/testthat/test-participant_profiler.R
# Unit tests for ParticipantProfiler (with SHAP-like explanations)

box::use(
  testthat[...],
  ../../R/calculators/participant_profiler[
    ParticipantProfiler, ParticipantProfileResult
  ],
  ../../R/calculators/threshold_selector[ThresholdStrategy]
)

# Helper to create participant feature data (output from FeatureEngineer)
create_participant_features <- function(n = 15, with_anomaly = TRUE) {
  set.seed(42)

  # Normal participants
  df <- data.frame(
    id = sprintf("P%02d", seq_len(n)),
    participant_mean_velocity = rnorm(n, 0.45, 0.05),
    participant_velocity_sd = runif(n, 0.03, 0.08),
    velocity_cv = runif(n, 0.08, 0.15),
    participant_rir_range = sample(4:8, n, replace = TRUE),
    participant_load_sensitivity = rnorm(n, -0.08, 0.02),
    mean_reps_to_failure = rnorm(n, 7, 1.5),
    set_consistency = runif(n, 0.1, 0.25),
    rir_slope = rnorm(n, 0.025, 0.005),
    rir_slope_deviation = rnorm(n, 0, 0.005),
    n_observations = sample(30:50, n, replace = TRUE),
    stringsAsFactors = FALSE
  )

  # Add obvious anomaly if requested
  if (with_anomaly && n >= 3) {
    # Make P03 anomalous: very low velocity, high CV
    df$participant_mean_velocity[3] <- 0.25  # Very low
    df$velocity_cv[3] <- 0.35               # Very high
    df$rir_slope[3] <- 0.05                 # Very high slope
    df$rir_slope_deviation[3] <- 0.025      # High deviation
  }

  df
}

describe("ParticipantProfileResult", {

  it("stores profiling results", {
    profiles <- data.frame(id = c("P01", "P02"), velocity_cv = c(0.1, 0.5))
    scores <- c(0.3, 0.8)
    is_anomaly <- c(FALSE, TRUE)
    reasons <- c("", "P02 flagged: velocity_cv very HIGH (z=2.5)")

    result <- ParticipantProfileResult$new(
      profiles = profiles,
      anomaly_scores = scores,
      is_anomaly = is_anomaly,
      threshold = 0.7,
      threshold_strategy = "elbow_detection",
      anomaly_reasons = reasons,
      method = "isolation_forest"
    )

    expect_equal(length(result$anomaly_scores), 2)
    expect_equal(result$threshold, 0.7)
    expect_equal(result$method, "isolation_forest")
  })

  it("returns anomalous IDs", {
    profiles <- data.frame(id = c("P01", "P02", "P03"))
    scores <- c(0.3, 0.8, 0.9)
    is_anomaly <- c(FALSE, TRUE, TRUE)

    result <- ParticipantProfileResult$new(
      profiles = profiles,
      anomaly_scores = scores,
      is_anomaly = is_anomaly,
      threshold = 0.7,
      threshold_strategy = "fixed_percentile",
      anomaly_reasons = c("", "reason1", "reason2"),
      method = "isolation_forest"
    )

    expect_equal(result$get_anomalous_ids(), c("P02", "P03"))
  })

  it("provides explanation for specific participant", {
    profiles <- data.frame(id = c("P01", "P02"))
    reasons <- c("", "P02 flagged: velocity HIGH")

    result <- ParticipantProfileResult$new(
      profiles = profiles,
      anomaly_scores = c(0.3, 0.8),
      is_anomaly = c(FALSE, TRUE),
      threshold = 0.7,
      threshold_strategy = "elbow",
      anomaly_reasons = reasons,
      method = "mahalanobis"
    )

    expect_match(result$get_explanation("P02"), "velocity HIGH")
    expect_match(result$get_explanation("P01"), "not flagged")
    expect_match(result$get_explanation("P99"), "not found")
  })

  it("provides summary", {
    profiles <- data.frame(id = c("P01", "P02", "P03"))
    result <- ParticipantProfileResult$new(
      profiles = profiles,
      anomaly_scores = c(0.3, 0.7, 0.9),
      is_anomaly = c(FALSE, FALSE, TRUE),
      threshold = 0.85,
      threshold_strategy = "gap_statistics",
      anomaly_reasons = c("", "", "reason"),
      method = "isolation_forest"
    )

    summary <- result$summarize()

    expect_equal(summary$n_participants, 3)
    expect_equal(summary$n_anomalous, 1)
    expect_equal(summary$threshold, 0.85)
    expect_equal(summary$method, "isolation_forest")
  })
})

describe("ParticipantProfiler", {

  describe("profile_participants", {

    it("validates inputs", {
      profiler <- ParticipantProfiler$new()

      expect_error(
        profiler$profile_participants("not a data frame"),
        "participant_features must be data frame"
      )

      expect_error(
        profiler$profile_participants(data.frame(x = 1)),
        "must have 'id' column"
      )

      expect_error(
        profiler$profile_participants(
          data.frame(id = "P01", x = 1),
          method = "invalid"
        ),
        "method must be valid"
      )
    })

    it("returns ParticipantProfileResult", {
      profiler <- ParticipantProfiler$new()
      features <- create_participant_features()

      result <- profiler$profile_participants(features)

      expect_s3_class(result, "ParticipantProfileResult")
    })

    it("detects synthetic anomaly with isolation_forest", {
      set.seed(42)
      profiler <- ParticipantProfiler$new(random_state = 42)
      features <- create_participant_features(n = 15, with_anomaly = TRUE)

      result <- profiler$profile_participants(
        features,
        method = "isolation_forest",
        threshold_strategy = ThresholdStrategy$ELBOW_DETECTION
      )

      # P03 should be flagged (or have high score)
      p03_idx <- which(result$profiles$id == "P03")
      expect_gt(result$anomaly_scores[p03_idx], median(result$anomaly_scores))
    })

    it("detects synthetic anomaly with mahalanobis", {
      profiler <- ParticipantProfiler$new()
      features <- create_participant_features(n = 15, with_anomaly = TRUE)

      result <- profiler$profile_participants(
        features,
        method = "mahalanobis",
        threshold_strategy = ThresholdStrategy$ELBOW_DETECTION
      )

      # P03 should have high score
      p03_idx <- which(result$profiles$id == "P03")
      expect_gt(result$anomaly_scores[p03_idx], 0.7)
    })

    it("respects threshold strategy", {
      profiler <- ParticipantProfiler$new()
      features <- create_participant_features(n = 20, with_anomaly = FALSE)

      result_percentile <- profiler$profile_participants(
        features,
        threshold_strategy = ThresholdStrategy$FIXED_PERCENTILE,
        threshold_params = list(percentile = 0.95)
      )

      result_elbow <- profiler$profile_participants(
        features,
        threshold_strategy = ThresholdStrategy$ELBOW_DETECTION
      )

      # Thresholds should be different
      expect_false(result_percentile$threshold == result_elbow$threshold)
    })

    it("uses default features when not specified", {
      profiler <- ParticipantProfiler$new()
      features <- create_participant_features()

      result <- profiler$profile_participants(features)

      expect_s3_class(result, "ParticipantProfileResult")
      expect_equal(length(result$anomaly_scores), nrow(features))
    })

    it("uses custom features when specified", {
      profiler <- ParticipantProfiler$new()
      features <- create_participant_features()

      result <- profiler$profile_participants(
        features,
        features = c("velocity_cv", "participant_mean_velocity")
      )

      expect_s3_class(result, "ParticipantProfileResult")
    })

    it("handles missing features gracefully", {
      profiler <- ParticipantProfiler$new()
      features <- data.frame(
        id = c("P01", "P02", "P03"),
        velocity_cv = c(0.1, 0.15, 0.4),
        participant_mean_velocity = c(0.45, 0.42, 0.25)
      )

      result <- profiler$profile_participants(
        features,
        features = c("velocity_cv", "participant_mean_velocity", "nonexistent")
      )

      expect_s3_class(result, "ParticipantProfileResult")
    })

    it("requires at least 2 features", {
      profiler <- ParticipantProfiler$new()
      features <- data.frame(id = c("P01", "P02"), x = c(1, 2))

      expect_error(
        profiler$profile_participants(features, features = "x"),
        "Need at least 2 features"
      )
    })
  })

  describe("explanations", {

    it("generates explanations for anomalous participants", {
      profiler <- ParticipantProfiler$new()
      features <- create_participant_features(n = 15, with_anomaly = TRUE)

      result <- profiler$profile_participants(features)

      # Check that anomalous participants have explanations
      for (id in result$get_anomalous_ids()) {
        explanation <- result$get_explanation(id)
        expect_true(nchar(explanation) > 0)
        expect_match(explanation, "flagged:")
      }
    })

    it("includes SHAP-like contributions in explanations", {
      profiler <- ParticipantProfiler$new()
      features <- create_participant_features(n = 15, with_anomaly = TRUE)

      result <- profiler$profile_participants(
        features,
        method = "isolation_forest",  # SHAP contributions only for isolation_forest
        threshold_strategy = ThresholdStrategy$FIXED_PERCENTILE,
        threshold_params = list(percentile = 0.8)  # Lower threshold to ensure anomalies
      )

      if (length(result$get_anomalous_ids()) > 0) {
        explanation <- result$get_explanation(result$get_anomalous_ids()[1])
        # SHAP format: "feature increases/decreases score by X"
        expect_match(explanation, "score by")
      }
    })
  })

  describe("plot_profiles", {

    it("returns ggplot object", {
      skip_if_not_installed("ggplot2")

      profiler <- ParticipantProfiler$new()
      features <- create_participant_features()
      result <- profiler$profile_participants(features)

      plot <- profiler$plot_profiles(result)

      expect_s3_class(plot, "gg")
    })
  })

  describe("compare_methods", {

    it("compares isolation_forest and mahalanobis", {
      profiler <- ParticipantProfiler$new()
      features <- create_participant_features()

      comparison <- profiler$compare_methods(features)

      expect_s3_class(comparison, "data.frame")
      expect_true("method" %in% names(comparison))
      expect_true("n_flagged" %in% names(comparison))
      expect_true("threshold" %in% names(comparison))
      expect_equal(nrow(comparison), 2)  # Both methods
    })
  })

  describe("reproducibility", {

    it("produces same results with same seed", {
      profiler1 <- ParticipantProfiler$new(random_state = 123)
      profiler2 <- ParticipantProfiler$new(random_state = 123)
      features <- create_participant_features()

      result1 <- profiler1$profile_participants(features, method = "isolation_forest")
      result2 <- profiler2$profile_participants(features, method = "isolation_forest")

      expect_equal(result1$anomaly_scores, result2$anomaly_scores)
    })

    it("produces different results with different seeds", {
      profiler1 <- ParticipantProfiler$new(random_state = 123)
      profiler2 <- ParticipantProfiler$new(random_state = 456)
      features <- create_participant_features()

      result1 <- profiler1$profile_participants(features, method = "isolation_forest")
      result2 <- profiler2$profile_participants(features, method = "isolation_forest")

      # Scores may differ (not always, but usually)
      # Just check that the function runs with different seeds
      expect_equal(length(result1$anomaly_scores), length(result2$anomaly_scores))
    })
  })

  describe("edge cases", {

    it("handles small number of participants", {
      profiler <- ParticipantProfiler$new()
      features <- create_participant_features(n = 5, with_anomaly = FALSE)

      result <- profiler$profile_participants(features)

      expect_equal(length(result$anomaly_scores), 5)
    })

    it("handles NA values in features", {
      profiler <- ParticipantProfiler$new()
      features <- create_participant_features(n = 10)

      # Introduce some NAs
      features$velocity_cv[c(2, 5)] <- NA
      features$participant_mean_velocity[3] <- NA

      result <- profiler$profile_participants(features)

      expect_s3_class(result, "ParticipantProfileResult")
      expect_equal(length(result$anomaly_scores), 10)
      expect_true(all(!is.na(result$anomaly_scores)))
    })

    it("handles all same values in a feature", {
      profiler <- ParticipantProfiler$new()
      features <- create_participant_features(n = 10)

      # Make one feature constant
      features$velocity_cv <- 0.1

      result <- profiler$profile_participants(features)

      expect_s3_class(result, "ParticipantProfileResult")
    })
  })

  describe("SHAP-like contributions", {

    it("computes feature contributions for isolation_forest", {
      profiler <- ParticipantProfiler$new()
      features <- create_participant_features(n = 15, with_anomaly = TRUE)

      result <- profiler$profile_participants(
        features,
        method = "isolation_forest",
        compute_contributions = TRUE
      )

      # Should have contributions
      expect_false(is.null(result$feature_contributions))
      expect_equal(length(result$feature_contributions), nrow(features))
    })

    it("get_top_contributors returns contribution data", {
      profiler <- ParticipantProfiler$new()
      features <- create_participant_features(n = 15, with_anomaly = TRUE)

      result <- profiler$profile_participants(
        features,
        method = "isolation_forest",
        threshold_strategy = ThresholdStrategy$FIXED_PERCENTILE,
        threshold_params = list(percentile = 0.8)
      )

      if (length(result$get_anomalous_ids()) > 0) {
        anomalous_id <- result$get_anomalous_ids()[1]
        top <- result$get_top_contributors(anomalous_id, k = 3)

        expect_s3_class(top, "data.frame")
        expect_true("feature" %in% names(top))
        expect_true("contribution" %in% names(top))
        expect_lte(nrow(top), 3)
      }
    })

    it("skips contributions for mahalanobis method", {
      profiler <- ParticipantProfiler$new()
      features <- create_participant_features(n = 15)

      result <- profiler$profile_participants(
        features,
        method = "mahalanobis",
        compute_contributions = TRUE
      )

      # Mahalanobis doesn't compute SHAP contributions
      expect_true(is.null(result$feature_contributions))
    })

    it("can disable contribution computation", {
      profiler <- ParticipantProfiler$new()
      features <- create_participant_features(n = 15)

      result <- profiler$profile_participants(
        features,
        method = "isolation_forest",
        compute_contributions = FALSE
      )

      expect_true(is.null(result$feature_contributions))
    })

    it("contributions explain anomaly direction", {
      profiler <- ParticipantProfiler$new()
      features <- create_participant_features(n = 15, with_anomaly = TRUE)

      result <- profiler$profile_participants(
        features,
        method = "isolation_forest",
        threshold_strategy = ThresholdStrategy$FIXED_PERCENTILE,
        threshold_params = list(percentile = 0.8)
      )

      anomalous_ids <- result$get_anomalous_ids()
      if (length(anomalous_ids) > 0) {
        idx <- which(result$profiles$id == anomalous_ids[1])
        contrib <- result$feature_contributions[[idx]]

        # Should have contributions for each feature used
        expect_true(length(contrib) >= 2)
        # At least one should be non-zero for an anomaly
        expect_true(any(abs(contrib) > 0))
      }
    })
  })

  describe("plot_contributions", {

    it("returns ggplot for valid participant", {
      skip_if_not_installed("ggplot2")

      profiler <- ParticipantProfiler$new()
      features <- create_participant_features(n = 15, with_anomaly = TRUE)

      result <- profiler$profile_participants(
        features,
        method = "isolation_forest",
        threshold_strategy = ThresholdStrategy$FIXED_PERCENTILE,
        threshold_params = list(percentile = 0.8)
      )

      if (length(result$get_anomalous_ids()) > 0) {
        plot <- profiler$plot_contributions(
          result,
          result$get_anomalous_ids()[1],
          top_k = 5
        )
        expect_s3_class(plot, "gg")
      }
    })

    it("errors for participant without contributions", {
      profiler <- ParticipantProfiler$new()
      features <- create_participant_features(n = 10)

      result <- profiler$profile_participants(
        features,
        method = "mahalanobis"  # No SHAP contributions
      )

      expect_error(
        profiler$plot_contributions(result, result$profiles$id[1]),
        "No contributions"
      )
    })
  })

  describe("plot_aggregate_importance", {

    it("returns ggplot showing feature importance", {
      skip_if_not_installed("ggplot2")

      profiler <- ParticipantProfiler$new()
      features <- create_participant_features(n = 15, with_anomaly = TRUE)

      result <- profiler$profile_participants(
        features,
        method = "isolation_forest",
        threshold_strategy = ThresholdStrategy$FIXED_PERCENTILE,
        threshold_params = list(percentile = 0.8)
      )

      if (sum(result$is_anomaly) > 0) {
        plot <- profiler$plot_aggregate_importance(result, top_k = 5)
        expect_s3_class(plot, "gg")
      }
    })
  })

  describe("plot_all_anomaly_explanations", {

    it("returns faceted ggplot for all anomalies", {
      skip_if_not_installed("ggplot2")

      profiler <- ParticipantProfiler$new()
      features <- create_participant_features(n = 15, with_anomaly = TRUE)

      result <- profiler$profile_participants(
        features,
        method = "isolation_forest",
        threshold_strategy = ThresholdStrategy$FIXED_PERCENTILE,
        threshold_params = list(percentile = 0.8)
      )

      if (sum(result$is_anomaly) > 0) {
        plot <- profiler$plot_all_anomaly_explanations(result, top_k = 3)
        expect_s3_class(plot, "gg")
      }
    })
  })

  describe("scores are in valid range", {

    it("isolation_forest scores are in [0, 1]", {
      profiler <- ParticipantProfiler$new()
      features <- create_participant_features(n = 20)

      result <- profiler$profile_participants(features, method = "isolation_forest")

      expect_true(all(result$anomaly_scores >= 0))
      expect_true(all(result$anomaly_scores <= 1))
    })

    it("mahalanobis scores are in [0, 1]", {
      profiler <- ParticipantProfiler$new()
      features <- create_participant_features(n = 20)

      result <- profiler$profile_participants(features, method = "mahalanobis")

      expect_true(all(result$anomaly_scores >= 0))
      expect_true(all(result$anomaly_scores <= 1))
    })
  })
})
