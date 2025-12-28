# tests/testthat/test-anomaly_detector.R
# Tests for AnomalyDetector using Isolation Forest

box::use(
  testthat[describe, it, expect_true, expect_s3_class, expect_type,
           expect_length, expect_gt, expect_equal, expect_lt, skip_if_not_installed]
)

# Set box path relative to project root
options(box.path = c("../../R", getOption("box.path")))

# Load module under test
box::use(
  calculators/anomaly_detector[
    AnomalyDetector, AnomalyResult, ParticipantAnomalyResult,
    detect_raw_anomalies, detect_re_anomalies
  ]
)

# =============================================================================
# TEST FIXTURES
# =============================================================================

#' Create test data with known outliers
create_test_data_with_outliers <- function(n_normal = 95, n_outliers = 5) {
  set.seed(42)

  # Normal data
  normal_data <- data.frame(
    x = rnorm(n_normal, mean = 0, sd = 1),
    y = rnorm(n_normal, mean = 0, sd = 1),
    z = rnorm(n_normal, mean = 0, sd = 1)
  )

  # Outlier data (shifted mean)
  outlier_data <- data.frame(
    x = rnorm(n_outliers, mean = 5, sd = 0.5),
    y = rnorm(n_outliers, mean = 5, sd = 0.5),
    z = rnorm(n_outliers, mean = 5, sd = 0.5)
  )

  rbind(normal_data, outlier_data)
}

#' Create test data mimicking deadlift study structure
create_deadlift_test_data <- function(n_participants = 10, obs_per_participant = 15) {
  set.seed(42)

  # Individual intercepts and slopes
  intercepts <- rnorm(n_participants, mean = 0.4, sd = 0.05)
  slopes <- rnorm(n_participants, mean = 0.025, sd = 0.008)

  data_list <- lapply(seq_len(n_participants), function(i) {
    rir <- sample(0:7, obs_per_participant, replace = TRUE)
    velocity <- intercepts[i] + slopes[i] * rir + rnorm(obs_per_participant, sd = 0.03)
    data.frame(
      id = paste0("P", i),
      rir = rir,
      mean_velocity = velocity,
      load_percentage = rep(runif(1, 70, 90), obs_per_participant)
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
# ANOMALY RESULT TESTS
# =============================================================================

describe("AnomalyResult", {

  it("stores anomaly detection results", {
    scores <- runif(100, 0, 1)
    threshold <- 0.8
    contamination <- 0.05

    result <- AnomalyResult$new(scores, threshold, contamination)

    expect_s3_class(result, "AnomalyResult")
    expect_length(result$scores, 100)
    expect_equal(result$threshold, 0.8)
    expect_equal(result$contamination_rate, 0.05)
  })

  it("calculates anomaly counts", {
    scores <- c(rep(0.3, 90), rep(0.9, 10))  # 10 outliers
    threshold <- 0.8

    result <- AnomalyResult$new(scores, threshold, 0.1)

    expect_equal(result$n_anomalies, 10)
    expect_equal(length(result$anomaly_indices), 10)
    expect_true(all(result$anomaly_indices >= 91))
  })

  it("summarizes results", {
    scores <- runif(100, 0, 1)
    result <- AnomalyResult$new(scores, 0.8, 0.05)

    summary <- result$summarize()

    expect_type(summary, "list")
    expect_true("n_total" %in% names(summary))
    expect_true("n_anomalies" %in% names(summary))
    expect_true("anomaly_rate" %in% names(summary))
  })
})

# =============================================================================
# PARTICIPANT ANOMALY RESULT TESTS
# =============================================================================

describe("ParticipantAnomalyResult", {

  it("stores participant-level results", {
    participant_summary <- data.frame(
      participant_id = paste0("P", 1:10),
      mean_score = seq(0.1, 1.0, by = 0.1),
      max_score = seq(0.15, 1.05, by = 0.1),
      n_anomalies = c(rep(0, 8), 1, 2),
      n_observations = rep(10, 10)
    )

    result <- ParticipantAnomalyResult$new(participant_summary, threshold = 0.8)

    expect_s3_class(result, "ParticipantAnomalyResult")
    expect_length(result$participant_id, 10)
  })

  it("classifies participants by interpretation", {
    participant_summary <- data.frame(
      participant_id = paste0("P", 1:5),
      mean_score = c(0.3, 0.5, 0.6, 0.85, 0.95),
      max_score = c(0.4, 0.6, 0.7, 0.9, 1.0),
      n_anomalies = c(0, 0, 1, 1, 2),
      n_observations = rep(10, 5)
    )

    result <- ParticipantAnomalyResult$new(participant_summary, threshold = 0.8)

    expect_true("normal" %in% result$interpretation)
    expect_true("anomalous" %in% result$interpretation)
  })

  it("gets anomalous participant IDs", {
    participant_summary <- data.frame(
      participant_id = paste0("P", 1:5),
      mean_score = c(0.3, 0.5, 0.6, 0.85, 0.95),
      max_score = c(0.4, 0.6, 0.7, 0.9, 1.0),
      n_anomalies = c(0, 0, 1, 1, 2),
      n_observations = rep(10, 5)
    )

    result <- ParticipantAnomalyResult$new(participant_summary, threshold = 0.8)
    anomalous_ids <- result$get_anomalous_ids()

    expect_type(anomalous_ids, "character")
    expect_true("P4" %in% anomalous_ids || "P5" %in% anomalous_ids)
  })
})

# =============================================================================
# ANOMALY DETECTOR TESTS
# =============================================================================

describe("AnomalyDetector", {

  describe("initialization", {
    it("creates instance with default parameters", {
      detector <- AnomalyDetector$new()
      expect_s3_class(detector, "AnomalyDetector")
    })

    it("accepts custom random state", {
      detector <- AnomalyDetector$new(random_state = 123)
      expect_s3_class(detector, "AnomalyDetector")
    })
  })

  describe("detect_raw_data_anomalies", {
    skip_if_not_installed("isotree")

    it("detects anomalies in raw data", {
      data <- create_test_data_with_outliers(n_normal = 95, n_outliers = 5)
      detector <- AnomalyDetector$new(random_state = 42)

      result <- detector$detect_raw_data_anomalies(
        data, c("x", "y", "z"), contamination = 0.05
      )

      expect_s3_class(result, "AnomalyResult")
      expect_length(result$scores, 100)
      expect_gt(result$n_anomalies, 0)
      expect_lt(result$n_anomalies, 20)  # Reasonable upper bound
    })

    it("respects contamination parameter", {
      data <- create_test_data_with_outliers()
      detector <- AnomalyDetector$new(random_state = 42)

      result_low <- detector$detect_raw_data_anomalies(
        data, c("x", "y"), contamination = 0.01
      )
      result_high <- detector$detect_raw_data_anomalies(
        data, c("x", "y"), contamination = 0.10
      )

      # Higher contamination = lower threshold = more observations flagged
      # (threshold is the (1-contamination) quantile, so more contamination = lower quantile)
      expect_gt(result_low$threshold, result_high$threshold)
    })

    it("validates input data", {
      detector <- AnomalyDetector$new()
      data <- data.frame(x = 1:10, y = 1:10)

      expect_error(
        detector$detect_raw_data_anomalies(data, c("x", "nonexistent")),
        "features must exist"
      )
    })

    it("produces reproducible results", {
      skip_if_not_installed("isotree")
      data <- create_test_data_with_outliers()

      result1 <- AnomalyDetector$new(random_state = 42)$detect_raw_data_anomalies(
        data, c("x", "y", "z"), 0.05
      )
      result2 <- AnomalyDetector$new(random_state = 42)$detect_raw_data_anomalies(
        data, c("x", "y", "z"), 0.05
      )

      expect_equal(result1$scores, result2$scores)
      expect_equal(result1$anomaly_indices, result2$anomaly_indices)
    })
  })

  describe("detect_random_effects_anomalies", {
    skip_if_not_installed("lme4")
    skip_if_not_installed("isotree")

    it("detects anomalies in random effects", {
      data <- create_deadlift_test_data()
      model <- create_test_model(data)
      detector <- AnomalyDetector$new(random_state = 42)

      result <- detector$detect_random_effects_anomalies(model, contamination = 0.1)

      expect_s3_class(result, "ParticipantAnomalyResult")
      expect_length(result$participant_id, 10)
      expect_true(is.logical(result$is_anomaly))
    })

    it("validates model type", {
      detector <- AnomalyDetector$new()

      expect_error(
        detector$detect_random_effects_anomalies("not a model"),
        "model must be an lmerMod"
      )
    })
  })

  describe("detect_cv_residual_anomalies", {
    skip_if_not_installed("isotree")

    it("detects anomalies in CV residuals", {
      # Simulate CV residuals
      set.seed(42)
      residuals <- c(
        rnorm(80, 0, 0.02),   # Normal residuals
        rnorm(10, 0.05, 0.01),  # Participant with larger residuals
        rnorm(10, 0, 0.02)
      )
      participant_ids <- rep(paste0("P", 1:10), each = 10)

      detector <- AnomalyDetector$new(random_state = 42)
      result <- detector$detect_cv_residual_anomalies(
        residuals, participant_ids, contamination = 0.1
      )

      expect_s3_class(result, "ParticipantAnomalyResult")
      expect_length(result$participant_id, 10)
    })
  })

  describe("compare_with_influence_diagnostics", {
    skip_if_not_installed("isotree")

    it("compares IF with Cook's D", {
      data <- create_test_data_with_outliers()
      detector <- AnomalyDetector$new(random_state = 42)

      if_result <- detector$detect_raw_data_anomalies(
        data, c("x", "y", "z"), contamination = 0.05
      )

      # Simulate Cook's D flags (some overlap, some different)
      cooks_flags <- rep(FALSE, 100)
      cooks_flags[c(96, 97, 98, 50, 51)] <- TRUE  # 3 shared, 2 unique

      comparison <- detector$compare_with_influence_diagnostics(if_result, cooks_flags)

      expect_type(comparison, "list")
      expect_true("both_flagged" %in% names(comparison))
      expect_true("only_isolation_forest" %in% names(comparison))
      expect_true("only_cooks_distance" %in% names(comparison))
      expect_true("agreement_rate" %in% names(comparison))
      expect_true("interpretation" %in% names(comparison))
    })

    it("handles perfect agreement", {
      set.seed(42)
      scores <- c(rep(0.3, 90), rep(0.9, 10))
      result <- AnomalyResult$new(scores, 0.8, 0.1)

      # Same flags as IF
      cooks_flags <- c(rep(FALSE, 90), rep(TRUE, 10))

      detector <- AnomalyDetector$new()
      comparison <- detector$compare_with_influence_diagnostics(result, cooks_flags)

      expect_equal(comparison$agreement_rate, 1.0)
      expect_equal(comparison$only_isolation_forest, 0)
      expect_equal(comparison$only_cooks_distance, 0)
    })
  })

  describe("plot_anomaly_scores", {
    skip_if_not_installed("isotree")
    skip_if_not_installed("ggplot2")

    it("creates anomaly score distribution plot", {
      data <- create_test_data_with_outliers()
      detector <- AnomalyDetector$new(random_state = 42)

      result <- detector$detect_raw_data_anomalies(
        data, c("x", "y", "z"), 0.05
      )
      plot <- detector$plot_anomaly_scores(result)

      expect_s3_class(plot, "ggplot")
    })
  })

  describe("plot_participant_anomalies", {
    skip_if_not_installed("lme4")
    skip_if_not_installed("isotree")
    skip_if_not_installed("ggplot2")

    it("creates participant anomaly plot", {
      data <- create_deadlift_test_data()
      model <- create_test_model(data)
      detector <- AnomalyDetector$new(random_state = 42)

      result <- detector$detect_random_effects_anomalies(model, 0.1)
      plot <- detector$plot_participant_anomalies(result)

      expect_s3_class(plot, "ggplot")
    })
  })
})

# =============================================================================
# CONVENIENCE FUNCTION TESTS
# =============================================================================

describe("Convenience functions", {
  skip_if_not_installed("isotree")

  it("detect_raw_anomalies works", {
    data <- create_test_data_with_outliers()

    result <- detect_raw_anomalies(data, c("x", "y", "z"))
    expect_s3_class(result, "AnomalyResult")
  })

  it("detect_re_anomalies works", {
    skip_if_not_installed("lme4")

    data <- create_deadlift_test_data()
    model <- create_test_model(data)

    result <- detect_re_anomalies(model)
    expect_s3_class(result, "ParticipantAnomalyResult")
  })
})
