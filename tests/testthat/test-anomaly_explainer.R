# tests/testthat/test-anomaly_explainer.R
# Unit tests for AnomalyExplainer (SHAP-like explanations)

box::use(
  testthat[...],
  ../../R/calculators/anomaly_explainer[
    AnomalyExplainer, AnomalyExplanation
  ]
)

# Helper to create test data with anomalies
create_test_data <- function(n = 50, n_anomalies = 3) {
  set.seed(42)

  # Normal observations
  df <- data.frame(
    id = seq_len(n),
    feature1 = rnorm(n, 0.5, 0.1),
    feature2 = rnorm(n, 0.3, 0.05),
    feature3 = rnorm(n, 0.8, 0.1)
  )

  # Create obvious anomalies
  if (n_anomalies > 0 && n >= n_anomalies) {
    anomaly_idx <- seq_len(n_anomalies)
    df$feature1[anomaly_idx] <- 1.5  # Very high
    df$feature2[anomaly_idx] <- 0.8  # Very high
  }

  df
}

describe("AnomalyExplanation", {

  it("stores explanation data", {
    contributions <- c(feature1 = 0.2, feature2 = -0.1, feature3 = 0.05)
    feature_values <- c(feature1 = 0.8, feature2 = 0.2, feature3 = 0.5)

    explanation <- AnomalyExplanation$new(
      observation_id = "obs_1",
      anomaly_score = 0.85,
      threshold = 0.7,
      contributions = contributions,
      feature_values = feature_values
    )

    expect_equal(explanation$observation_id, "obs_1")
    expect_equal(explanation$anomaly_score, 0.85)
    expect_equal(explanation$threshold, 0.7)
    expect_equal(length(explanation$contributions), 3)
  })

  it("returns top K contributors", {
    contributions <- c(feature1 = 0.25, feature2 = -0.15, feature3 = 0.05)
    feature_values <- c(feature1 = 0.8, feature2 = 0.2, feature3 = 0.5)

    explanation <- AnomalyExplanation$new(
      observation_id = "obs_1",
      anomaly_score = 0.85,
      threshold = 0.7,
      contributions = contributions,
      feature_values = feature_values
    )

    top <- explanation$get_top_contributors(k = 2)

    expect_equal(nrow(top), 2)
    expect_true("feature" %in% names(top))
    expect_true("contribution" %in% names(top))
    # Should be sorted by absolute contribution
    expect_true(abs(top$contribution[1]) >= abs(top$contribution[2]))
  })

  it("generates human-readable explanation", {
    contributions <- c(feature1 = 0.25, feature2 = -0.15, feature3 = 0.05)
    feature_values <- c(feature1 = 0.8, feature2 = 0.2, feature3 = 0.5)

    explanation <- AnomalyExplanation$new(
      observation_id = "obs_1",
      anomaly_score = 0.85,
      threshold = 0.7,
      contributions = contributions,
      feature_values = feature_values
    )

    text <- explanation$explain(k = 2)

    expect_true(is.character(text))
    expect_true(grepl("obs_1", text))
    expect_true(grepl("score", text, ignore.case = TRUE))
    expect_true(grepl("increases|decreases", text, ignore.case = TRUE))
  })
})

describe("AnomalyExplainer", {

  describe("explain_anomalies", {

    it("returns list of AnomalyExplanation objects", {
      explainer <- AnomalyExplainer$new()
      data <- create_test_data(n = 30, n_anomalies = 2)
      feature_names <- c("feature1", "feature2", "feature3")

      explanations <- explainer$explain_anomalies(
        data,
        feature_names = feature_names,
        threshold = 0.6,  # Ensure some anomalies
        id_column = "id"
      )

      expect_true(length(explanations) > 0)
      expect_s3_class(explanations[[1]], "AnomalyExplanation")
    })

    it("identifies correct anomaly indices", {
      explainer <- AnomalyExplainer$new()
      data <- create_test_data(n = 30, n_anomalies = 3)
      feature_names <- c("feature1", "feature2", "feature3")

      # Use specific indices
      explanations <- explainer$explain_anomalies(
        data,
        feature_names = feature_names,
        anomaly_indices = c(1, 2, 3),  # Our created anomalies
        threshold = 0.7,
        id_column = "id"
      )

      expect_equal(length(explanations), 3)
    })

    it("computes feature contributions", {
      explainer <- AnomalyExplainer$new()
      data <- create_test_data(n = 30, n_anomalies = 2)
      feature_names <- c("feature1", "feature2", "feature3")

      explanations <- explainer$explain_anomalies(
        data,
        feature_names = feature_names,
        anomaly_indices = c(1),
        threshold = 0.7,
        id_column = "id"
      )

      exp <- explanations[[1]]

      # Should have contributions for all features
      expect_equal(length(exp$contributions), 3)
      expect_true(all(names(exp$contributions) %in% feature_names))
    })

    it("uses row numbers as IDs when id_column is NULL", {
      explainer <- AnomalyExplainer$new()
      data <- create_test_data(n = 20, n_anomalies = 1)
      data$id <- NULL  # Remove id column
      feature_names <- c("feature1", "feature2", "feature3")

      explanations <- explainer$explain_anomalies(
        data,
        feature_names = feature_names,
        anomaly_indices = c(1),
        threshold = 0.7
      )

      # Should use row number as ID
      expect_equal(explanations[[1]]$observation_id, "1")
    })

    it("returns empty list when no anomalies", {
      explainer <- AnomalyExplainer$new()
      # Create data with no anomalies
      set.seed(123)
      data <- data.frame(
        id = 1:20,
        feature1 = rnorm(20, 0, 0.01),
        feature2 = rnorm(20, 0, 0.01),
        feature3 = rnorm(20, 0, 0.01)
      )
      feature_names <- c("feature1", "feature2", "feature3")

      explanations <- explainer$explain_anomalies(
        data,
        feature_names = feature_names,
        threshold = 0.99,  # Very high threshold
        id_column = "id"
      )

      expect_equal(length(explanations), 0)
    })

    it("requires at least 2 features", {
      explainer <- AnomalyExplainer$new()
      data <- data.frame(id = 1:10, feature1 = rnorm(10))

      expect_error(
        explainer$explain_anomalies(data, feature_names = "feature1"),
        "at least 2"
      )
    })
  })

  describe("summarize_explanations", {

    it("returns data frame with all explanations", {
      explainer <- AnomalyExplainer$new()
      data <- create_test_data(n = 30, n_anomalies = 3)
      feature_names <- c("feature1", "feature2", "feature3")

      explanations <- explainer$explain_anomalies(
        data,
        feature_names = feature_names,
        anomaly_indices = c(1, 2),
        threshold = 0.7,
        id_column = "id"
      )

      summary <- explainer$summarize_explanations(explanations, top_k = 2)

      expect_s3_class(summary, "data.frame")
      expect_true("observation_id" %in% names(summary))
      expect_true("feature" %in% names(summary))
      expect_true("contribution" %in% names(summary))
    })

    it("returns empty data frame for no explanations", {
      explainer <- AnomalyExplainer$new()

      summary <- explainer$summarize_explanations(list(), top_k = 3)

      expect_equal(nrow(summary), 0)
    })
  })

  describe("plot_explanation", {

    it("returns ggplot object", {
      skip_if_not_installed("ggplot2")

      explainer <- AnomalyExplainer$new()
      data <- create_test_data(n = 30, n_anomalies = 2)
      feature_names <- c("feature1", "feature2", "feature3")

      explanations <- explainer$explain_anomalies(
        data,
        feature_names = feature_names,
        anomaly_indices = c(1),
        threshold = 0.7,
        id_column = "id"
      )

      plot <- explainer$plot_explanation(explanations[[1]], top_k = 3)

      expect_s3_class(plot, "gg")
    })

    it("validates input is AnomalyExplanation", {
      explainer <- AnomalyExplainer$new()

      expect_error(
        explainer$plot_explanation(list()),
        "AnomalyExplanation"
      )
    })
  })

  describe("plot_aggregate_importance", {

    it("returns ggplot object", {
      skip_if_not_installed("ggplot2")

      explainer <- AnomalyExplainer$new()
      data <- create_test_data(n = 30, n_anomalies = 3)
      feature_names <- c("feature1", "feature2", "feature3")

      explanations <- explainer$explain_anomalies(
        data,
        feature_names = feature_names,
        anomaly_indices = c(1, 2, 3),
        threshold = 0.7,
        id_column = "id"
      )

      plot <- explainer$plot_aggregate_importance(explanations)

      expect_s3_class(plot, "gg")
    })

    it("errors when no explanations provided", {
      explainer <- AnomalyExplainer$new()

      expect_error(
        explainer$plot_aggregate_importance(list()),
        "No explanations"
      )
    })
  })

  describe("reproducibility", {

    it("produces same results with same seed", {
      explainer1 <- AnomalyExplainer$new(random_state = 123)
      explainer2 <- AnomalyExplainer$new(random_state = 123)
      data <- create_test_data(n = 30, n_anomalies = 2)
      feature_names <- c("feature1", "feature2", "feature3")

      explanations1 <- explainer1$explain_anomalies(
        data,
        feature_names = feature_names,
        anomaly_indices = c(1),
        threshold = 0.7,
        id_column = "id"
      )

      explanations2 <- explainer2$explain_anomalies(
        data,
        feature_names = feature_names,
        anomaly_indices = c(1),
        threshold = 0.7,
        id_column = "id"
      )

      expect_equal(
        explanations1[[1]]$contributions,
        explanations2[[1]]$contributions
      )
    })
  })

  describe("contribution properties", {

    it("positive contribution means feature increases anomaly score", {
      explainer <- AnomalyExplainer$new(random_state = 42)

      # Create data where feature1 is clearly anomalous
      set.seed(42)
      data <- data.frame(
        id = 1:20,
        feature1 = c(2.0, rnorm(19, 0, 0.1)),  # First obs is clearly anomalous
        feature2 = rnorm(20, 0, 0.1),
        feature3 = rnorm(20, 0, 0.1)
      )
      feature_names <- c("feature1", "feature2", "feature3")

      explanations <- explainer$explain_anomalies(
        data,
        feature_names = feature_names,
        anomaly_indices = c(1),
        threshold = 0.5,
        id_column = "id"
      )

      exp <- explanations[[1]]

      # feature1 should have highest absolute contribution
      sorted_contributions <- sort(abs(exp$contributions), decreasing = TRUE)
      expect_equal(names(sorted_contributions)[1], "feature1")
    })
  })

  describe("threshold computation", {

    it("auto-computes threshold when not provided", {
      explainer <- AnomalyExplainer$new()
      data <- create_test_data(n = 50, n_anomalies = 5)
      feature_names <- c("feature1", "feature2", "feature3")

      # Don't provide threshold - should auto-compute
      explanations <- explainer$explain_anomalies(
        data,
        feature_names = feature_names,
        id_column = "id"
      )

      # Should have found at least one anomaly
      expect_true(length(explanations) >= 1)

      # Each explanation should have a threshold
      if (length(explanations) > 0) {
        expect_true(!is.na(explanations[[1]]$threshold))
        expect_true(explanations[[1]]$threshold > 0)
      }
    })
  })
})
