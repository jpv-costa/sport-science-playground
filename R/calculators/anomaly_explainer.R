# R/calculators/anomaly_explainer.R
# Service: SHAP-like Explanations for Isolation Forest Anomalies
#
# =============================================================================
# APPROACH:
# For Isolation Forest, we compute feature contributions using permutation:
# 1. Get baseline anomaly score for each observation
# 2. For each feature, replace with population median
# 3. Compute perturbed anomaly score
# 4. Contribution = baseline_score - perturbed_score
#
# Positive contribution = feature increases anomaly score (makes it more anomalous)
# Negative contribution = feature decreases anomaly score (makes it less anomalous)
#
# This is similar to SHAP's marginal contribution concept.
# =============================================================================

box::use(
  R6[R6Class],
  stats[median, quantile, predict],
  utils[head],
  isotree[isolation.forest, predict.isolation_forest],
  dplyr[mutate, arrange, filter, select, .data, slice_head, bind_rows],
  ggplot2[
    ggplot, aes, geom_col, geom_hline, geom_text, coord_flip,
    labs, theme_minimal, theme, element_text, scale_fill_gradient2, .data
  ]
)

#' Anomaly Explanation Result
#'
#' Contains feature contributions for anomalous observations.
#'
#' @export
AnomalyExplanation <- R6Class(
  classname = "AnomalyExplanation",
  cloneable = FALSE,

  public = list(
    #' @field observation_id ID of the observation (row index or actual ID)
    observation_id = NULL,
    #' @field anomaly_score Original anomaly score
    anomaly_score = NULL,
    #' @field threshold Threshold used for flagging
    threshold = NULL,
    #' @field contributions Named vector of feature contributions
    contributions = NULL,
    #' @field feature_values Named vector of original feature values
    feature_values = NULL,

    #' @description Create explanation
    initialize = function(observation_id, anomaly_score, threshold,
                          contributions, feature_values) {
      self$observation_id <- observation_id
      self$anomaly_score <- anomaly_score
      self$threshold <- threshold
      self$contributions <- contributions
      self$feature_values <- feature_values
    },

    #' @description Get top K contributing features
    #' @param k Number of features to return
    #' @return Data frame with feature, contribution, value
    get_top_contributors = function(k = 5) {
      sorted_idx <- order(abs(self$contributions), decreasing = TRUE)
      top_idx <- head(sorted_idx, k)

      data.frame(
        feature = names(self$contributions)[top_idx],
        contribution = self$contributions[top_idx],
        value = self$feature_values[top_idx],
        stringsAsFactors = FALSE
      )
    },

    #' @description Generate human-readable explanation
    #' @param k Number of features to include
    #' @return Character string
    explain = function(k = 3) {
      top <- self$get_top_contributors(k)

      parts <- sapply(seq_len(nrow(top)), function(i) {
        direction <- if (top$contribution[i] > 0) "increases" else "decreases"
        sprintf("%s (%.3f) %s anomaly score by %.3f",
                top$feature[i], top$value[i], direction, abs(top$contribution[i]))
      })

      paste0(
        sprintf("Observation %s (score: %.3f, threshold: %.3f)\n",
                self$observation_id, self$anomaly_score, self$threshold),
        "Top contributors:\n",
        paste("  -", parts, collapse = "\n")
      )
    }
  )
)

#' Anomaly Explainer
#'
#' Computes SHAP-like feature contributions for Isolation Forest anomalies.
#' Uses permutation-based importance to estimate each feature's contribution.
#'
#' @export
AnomalyExplainer <- R6Class(
  classname = "AnomalyExplainer",
  cloneable = FALSE,

  public = list(
    #' @description Create explainer
    #' @param random_state Seed for reproducibility
    #' @param n_trees Number of trees in Isolation Forest
    initialize = function(random_state = 42, n_trees = 100) {
      private$.random_state <- random_state
      private$.n_trees <- n_trees
    },

    #' @description Explain anomalies in data
    #' @param data Data frame with features
    #' @param feature_names Features to use for detection (or NULL for all numeric)
    #' @param anomaly_indices Indices of observations to explain (or NULL for all flagged)
    #' @param threshold Anomaly threshold (or NULL to compute via elbow)
    #' @param id_column Name of ID column (or NULL to use row numbers)
    #' @return List of AnomalyExplanation objects
    explain_anomalies = function(data, feature_names = NULL,
                                  anomaly_indices = NULL,
                                  threshold = NULL,
                                  id_column = NULL) {
      # Determine features
      if (is.null(feature_names)) {
        feature_names <- names(data)[sapply(data, is.numeric)]
      }
      feature_names <- intersect(feature_names, names(data))

      if (length(feature_names) < 2) {
        stop("Need at least 2 numeric features for explanation")
      }

      # Extract feature matrix
      feature_matrix <- private$.prepare_matrix(data, feature_names)

      # Fit Isolation Forest
      iso_model <- private$.fit_model(feature_matrix)

      # Get baseline scores
      baseline_scores <- predict(iso_model, feature_matrix)

      # Determine threshold
      if (is.null(threshold)) {
        threshold <- private$.compute_elbow_threshold(baseline_scores)
      }

      # Determine which observations to explain
      if (is.null(anomaly_indices)) {
        anomaly_indices <- which(baseline_scores >= threshold)
      }

      if (length(anomaly_indices) == 0) {
        message("No anomalies to explain")
        return(list())
      }

      # Compute explanations for each anomaly
      explanations <- lapply(anomaly_indices, function(idx) {
        private$.explain_single(
          idx = idx,
          feature_matrix = feature_matrix,
          feature_names = feature_names,
          iso_model = iso_model,
          baseline_score = baseline_scores[idx],
          threshold = threshold,
          id_column = id_column,
          data = data
        )
      })

      names(explanations) <- as.character(anomaly_indices)
      explanations
    },

    #' @description Get summary of all explanations
    #' @param explanations List of AnomalyExplanation objects
    #' @param top_k Top K features per explanation
    #' @return Data frame with all explanations
    summarize_explanations = function(explanations, top_k = 3) {
      if (length(explanations) == 0) return(data.frame())

      result <- lapply(names(explanations), function(idx) {
        exp <- explanations[[idx]]
        top <- exp$get_top_contributors(top_k)
        top$observation_id <- exp$observation_id
        top$anomaly_score <- exp$anomaly_score
        top
      })

      do.call(bind_rows, result)
    },

    #' @description Plot feature contributions for single observation
    #' @param explanation AnomalyExplanation object
    #' @param top_k Number of features to show
    #' @return ggplot object
    plot_explanation = function(explanation, top_k = 10) {
      stopifnot(inherits(explanation, "AnomalyExplanation"))

      df <- explanation$get_top_contributors(top_k)
      df$feature <- factor(df$feature, levels = rev(df$feature))

      ggplot(df, aes(x = .data$feature, y = .data$contribution,
                     fill = .data$contribution)) +
        geom_col(width = 0.7) +
        geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
        coord_flip() +
        scale_fill_gradient2(
          low = "#2E86AB", mid = "white", high = "#E63946",
          midpoint = 0, guide = "none"
        ) +
        labs(
          title = sprintf("Feature Contributions - %s", explanation$observation_id),
          subtitle = sprintf("Anomaly Score: %.3f | Threshold: %.3f",
                            explanation$anomaly_score, explanation$threshold),
          x = NULL,
          y = "Contribution to Anomaly Score"
        ) +
        theme_minimal(base_size = 11) +
        theme(
          plot.title = element_text(face = "bold"),
          axis.text.y = element_text(size = 9)
        )
    },

    #' @description Plot aggregate feature importance across all explanations
    #' @param explanations List of AnomalyExplanation objects
    #' @return ggplot object
    plot_aggregate_importance = function(explanations) {
      if (length(explanations) == 0) {
        stop("No explanations to plot")
      }

      # Aggregate contributions
      all_contrib <- lapply(explanations, function(exp) {
        data.frame(
          feature = names(exp$contributions),
          contribution = abs(exp$contributions),
          stringsAsFactors = FALSE
        )
      })

      agg <- do.call(bind_rows, all_contrib) |>
        dplyr::group_by(.data$feature) |>
        dplyr::summarize(
          mean_abs_contribution = mean(.data$contribution, na.rm = TRUE),
          .groups = "drop"
        ) |>
        arrange(dplyr::desc(.data$mean_abs_contribution)) |>
        slice_head(n = 15)

      agg$feature <- factor(agg$feature, levels = rev(agg$feature))

      ggplot(agg, aes(x = .data$feature, y = .data$mean_abs_contribution)) +
        geom_col(fill = "#2E86AB", width = 0.7) +
        coord_flip() +
        labs(
          title = "Aggregate Feature Importance",
          subtitle = sprintf("Mean |contribution| across %d anomalies", length(explanations)),
          x = NULL,
          y = "Mean Absolute Contribution"
        ) +
        theme_minimal(base_size = 11) +
        theme(plot.title = element_text(face = "bold"))
    }
  ),

  private = list(
    .random_state = NULL,
    .n_trees = NULL,

    .prepare_matrix = function(data, feature_names) {
      mat <- as.matrix(data[, feature_names, drop = FALSE])
      # Impute NAs with column medians
      for (j in seq_len(ncol(mat))) {
        na_idx <- is.na(mat[, j])
        if (any(na_idx)) {
          mat[na_idx, j] <- median(mat[!na_idx, j], na.rm = TRUE)
        }
      }
      mat
    },

    .fit_model = function(feature_matrix) {
      isolation.forest(
        feature_matrix,
        ntrees = private$.n_trees,
        seed = private$.random_state,
        nthreads = 1
      )
    },

    .compute_elbow_threshold = function(scores) {
      sorted <- sort(scores)
      n <- length(sorted)
      if (n < 10) return(quantile(scores, 0.95))

      # Simple elbow detection
      x <- seq_len(n)
      y <- sorted

      # Line from first to last point
      p1 <- c(1, y[1])
      p2 <- c(n, y[n])

      # Distance from each point to line
      distances <- abs((p2[2] - p1[2]) * x - (p2[1] - p1[1]) * y +
                       p2[1] * p1[2] - p2[2] * p1[1]) /
                   sqrt((p2[2] - p1[2])^2 + (p2[1] - p1[1])^2)

      elbow_idx <- which.max(distances)
      sorted[elbow_idx]
    },

    .explain_single = function(idx, feature_matrix, feature_names, iso_model,
                                baseline_score, threshold, id_column, data) {
      n_features <- length(feature_names)
      contributions <- numeric(n_features)
      names(contributions) <- feature_names

      # Compute median for each feature (baseline for perturbation)
      medians <- apply(feature_matrix, 2, median, na.rm = TRUE)

      # For each feature, compute contribution via perturbation
      for (j in seq_len(n_features)) {
        # Create perturbed data (replace feature j with median)
        perturbed <- feature_matrix[idx, , drop = FALSE]
        perturbed[1, j] <- medians[j]

        # Get perturbed score
        perturbed_score <- predict(iso_model, perturbed)

        # Contribution = how much removing this feature's influence changes the score
        contributions[j] <- baseline_score - perturbed_score
      }

      # Get observation ID
      if (!is.null(id_column) && id_column %in% names(data)) {
        obs_id <- data[[id_column]][idx]
      } else {
        obs_id <- as.character(idx)
      }

      # Get original feature values
      feature_values <- feature_matrix[idx, ]
      names(feature_values) <- feature_names

      AnomalyExplanation$new(
        observation_id = obs_id,
        anomaly_score = baseline_score,
        threshold = threshold,
        contributions = contributions,
        feature_values = feature_values
      )
    }
  )
)
