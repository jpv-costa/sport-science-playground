# R/calculators/participant_profiler.R
# Service: Participant-Level Anomaly Detection via Aggregate Behavior Patterns
#
# =============================================================================
# PRINCIPLES APPLIED (from CLAUDE.md)
# =============================================================================
#
# NAMING PRINCIPLES:
# - Understandability: Profiler describes the action (building profiles)
# - Consistency: Methods follow verb-based naming (profile_, get_, plot_)
# - Distinguishability: Clear separation from observation-level detection
# - Conciseness: Short but meaningful names
#
# SOLID PRINCIPLES:
# - SRP: Single responsibility - only profiles participants
# - OCP: Open for extension with new profiling methods
# - DIP: Depends on abstractions (data frames, ThresholdSelector interface)
#
# CUPID PRINCIPLES:
# - Composable: Works with FeatureEngineer output
# - Unix: Profiles participants based on aggregate behavior
# - Predictable: Same inputs + seed -> same outputs
# - Idiomatic: Follows R6/R conventions
# - Domain-based: Names reflect athlete profiling concepts
#
# SCIENTIFIC VALIDITY:
# - Isolation Forest appropriate for multivariate anomaly detection
# - Mahalanobis distance provides interpretable statistical distance
# - Z-score explanations grounded in statistical theory
# =============================================================================

box::use(
  R6[R6Class],
  stats[mahalanobis, cov, median, sd, pchisq, predict],
  utils[head],
  grid[unit],
  isotree[isolation.forest, predict.isolation_forest],
  ggplot2[
    ggplot, aes, geom_point, geom_hline, geom_text, geom_col, geom_segment,
    labs, coord_flip, facet_wrap, scale_fill_gradient2,
    theme_minimal, theme, element_text, element_blank, scale_color_manual, .data
  ],
  dplyr[bind_rows, group_by, summarize, arrange, desc, slice_head, .data],
  ./threshold_selector[ThresholdSelector, ThresholdStrategy]
)

#' Participant Profile Result Value Object
#'
#' Stores participant-level anomaly detection results with SHAP-like explanations.
#'
#' @export
ParticipantProfileResult <- R6Class(
  classname = "ParticipantProfileResult",
  cloneable = FALSE,

  public = list(
    #' @field profiles Data frame with all participant feature profiles
    profiles = NULL,
    #' @field anomaly_scores Numeric vector of scores per participant
    anomaly_scores = NULL,
    #' @field is_anomaly Logical vector of anomaly flags
    is_anomaly = NULL,
    #' @field threshold Threshold used for flagging
    threshold = NULL,
    #' @field threshold_strategy Strategy used to select threshold
    threshold_strategy = NULL,
    #' @field anomaly_reasons Character vector with explanations
    anomaly_reasons = NULL,
    #' @field method Method used (isolation_forest or mahalanobis)
    method = NULL,
    #' @field feature_contributions List of named vectors (SHAP-like contributions per participant)
    feature_contributions = NULL,

    #' @description Create ParticipantProfileResult
    #' @param profiles Data frame with participant features
    #' @param anomaly_scores Numeric scores per participant
    #' @param is_anomaly Logical flags
    #' @param threshold Threshold value
    #' @param threshold_strategy Strategy name
    #' @param anomaly_reasons Explanation strings
    #' @param method Detection method used
    #' @param feature_contributions List of SHAP-like contributions per participant
    initialize = function(profiles, anomaly_scores, is_anomaly,
                          threshold, threshold_strategy, anomaly_reasons, method,
                          feature_contributions = NULL) {
      self$profiles <- profiles
      self$anomaly_scores <- anomaly_scores
      self$is_anomaly <- is_anomaly
      self$threshold <- threshold
      self$threshold_strategy <- threshold_strategy
      self$anomaly_reasons <- anomaly_reasons
      self$method <- method
      self$feature_contributions <- feature_contributions
    },

    #' @description Get IDs of anomalous participants
    #' @return Character vector of participant IDs
    get_anomalous_ids = function() {
      self$profiles$id[self$is_anomaly]
    },

    #' @description Get detailed explanation for a specific participant
    #' @param participant_id ID of participant to explain
    #' @return Character string with explanation
    get_explanation = function(participant_id) {
      idx <- which(self$profiles$id == participant_id)
      if (length(idx) == 0) {
        return(paste("Participant", participant_id, "not found"))
      }
      if (!self$is_anomaly[idx]) {
        return(paste(participant_id, "is not flagged as anomalous"))
      }
      self$anomaly_reasons[idx]
    },

    #' @description Get top K contributing features for a participant
    #' @param participant_id ID of participant
    #' @param k Number of top features to return
    #' @return Data frame with feature, contribution, value
    get_top_contributors = function(participant_id, k = 5) {
      idx <- which(self$profiles$id == participant_id)
      if (length(idx) == 0 || is.null(self$feature_contributions)) {
        return(data.frame(feature = character(), contribution = numeric(),
                          value = numeric(), stringsAsFactors = FALSE))
      }

      contributions <- self$feature_contributions[[idx]]
      if (is.null(contributions)) {
        return(data.frame(feature = character(), contribution = numeric(),
                          value = numeric(), stringsAsFactors = FALSE))
      }

      sorted_idx <- order(abs(contributions), decreasing = TRUE)
      top_idx <- head(sorted_idx, k)

      data.frame(
        feature = names(contributions)[top_idx],
        contribution = contributions[top_idx],
        stringsAsFactors = FALSE
      )
    },

    #' @description Summary of profiling results
    #' @return Named list with key metrics
    summarize = function() {
      list(
        n_participants = nrow(self$profiles),
        n_anomalous = sum(self$is_anomaly),
        anomalous_ids = self$get_anomalous_ids(),
        threshold = self$threshold,
        threshold_strategy = self$threshold_strategy,
        method = self$method,
        score_range = range(self$anomaly_scores, na.rm = TRUE)
      )
    }
  )
)

#' Participant Profiler
#'
#' Detects participants with unusual aggregate behavior patterns.
#' Uses Isolation Forest or Mahalanobis distance on participant-level features.
#'
#' @section Detection Methods:
#' - isolation_forest: Multivariate Isolation Forest on participant features
#' - mahalanobis: Robust Mahalanobis distance (normalized via chi-squared CDF)
#'
#' @export
ParticipantProfiler <- R6Class(

  classname = "ParticipantProfiler",
  cloneable = FALSE,

  public = list(

    #' @description Create profiler
    #' @param random_state Seed for reproducibility (default: 42)
    #' @param n_trees Number of trees in Isolation Forest (default: 100)
    initialize = function(random_state = 42, n_trees = 100) {
      private$.random_state <- random_state
      private$.n_trees <- n_trees
    },

    #' @description Profile participants and detect anomalies
    #' @param participant_features Data frame from FeatureEngineer with id column
    #' @param features Character vector of feature names to use (or NULL for defaults)
    #' @param method "isolation_forest" or "mahalanobis"
    #' @param threshold_strategy ThresholdStrategy value or string
    #' @param threshold_params List of threshold parameters
    #' @param compute_contributions Whether to compute SHAP-like contributions
    #' @return ParticipantProfileResult object
    profile_participants = function(participant_features,
                                     features = NULL,
                                     method = "isolation_forest",
                                     threshold_strategy = ThresholdStrategy$ELBOW_DETECTION,
                                     threshold_params = list(),
                                     compute_contributions = TRUE) {
      stopifnot(
        "participant_features must be data frame" = is.data.frame(participant_features),
        "participant_features must have 'id' column" = "id" %in% names(participant_features),
        "method must be valid" = method %in% c("isolation_forest", "mahalanobis")
      )

      # Default features if not specified
      if (is.null(features)) {
        features <- c(
          "p_mean_velocity", "p_velocity_cv", "p_velocity_range",
          "p_rir_slope", "p_rir_r2", "p_load_sensitivity",
          "p_mean_reps_per_set", "p_reps_cv"
        )
      }

      # Filter to available features
      available <- intersect(features, names(participant_features))
      if (length(available) < 2) {
        stop("Need at least 2 features for profiling. Available: ",
             paste(intersect(features, names(participant_features)), collapse = ", "))
      }

      # Extract feature matrix
      feature_matrix <- as.matrix(participant_features[, available, drop = FALSE])

      # Impute NAs with column medians
      feature_matrix <- private$.impute_nas(feature_matrix)

      # Store model for contribution computation
      iso_model <- NULL

      # Calculate anomaly scores based on method
      if (method == "isolation_forest") {
        iso_result <- private$.isolation_forest_with_model(feature_matrix)
        scores <- iso_result$scores
        iso_model <- iso_result$model
      } else {
        scores <- private$.mahalanobis_scores(feature_matrix)
      }

      # Select threshold
      threshold_selector <- ThresholdSelector$new()
      threshold_result <- threshold_selector$select_threshold(
        scores, threshold_strategy, threshold_params
      )
      threshold <- threshold_result$threshold

      # Flag anomalies
      is_anomaly <- scores >= threshold

      # Compute SHAP-like contributions for all participants
      contributions <- NULL
      if (compute_contributions && method == "isolation_forest" && !is.null(iso_model)) {
        contributions <- private$.compute_contributions(
          feature_matrix, available, iso_model, scores
        )
      }

      # Generate explanations using contributions
      reasons <- private$.generate_shap_explanations(
        participant_features, available, is_anomaly, contributions
      )

      ParticipantProfileResult$new(
        profiles = participant_features,
        anomaly_scores = scores,
        is_anomaly = is_anomaly,
        threshold = threshold,
        threshold_strategy = threshold_strategy,
        anomaly_reasons = reasons,
        method = method,
        feature_contributions = contributions
      )
    },

    #' @description Plot participant profiles with anomalies highlighted
    #' @param result ParticipantProfileResult object
    #' @param title Plot title
    #' @return ggplot object
    plot_profiles = function(result, title = "Participant Anomaly Profiles") {
      stopifnot(
        "result must be ParticipantProfileResult" =
          inherits(result, "ParticipantProfileResult")
      )

      df <- data.frame(
        id = result$profiles$id,
        score = result$anomaly_scores,
        is_anomaly = factor(
          ifelse(result$is_anomaly, "Anomalous", "Normal"),
          levels = c("Normal", "Anomalous")
        )
      )

      # Order by score
      df <- df[order(df$score), ]
      df$rank <- seq_len(nrow(df))

      ggplot(df, aes(x = .data$rank, y = .data$score, color = .data$is_anomaly)) +
        geom_hline(
          yintercept = result$threshold,
          linetype = "dashed",
          color = "#E63946",
          linewidth = 1
        ) +
        geom_point(size = 4, alpha = 0.8) +
        geom_text(
          aes(label = .data$id),
          vjust = -1, size = 3, color = "gray30"
        ) +
        scale_color_manual(values = c("Normal" = "#2E86AB", "Anomalous" = "#E63946")) +
        labs(
          title = title,
          subtitle = sprintf(
            "Method: %s | Threshold: %.3f (%s) | %d of %d flagged",
            result$method,
            result$threshold,
            result$threshold_strategy,
            sum(result$is_anomaly),
            length(result$is_anomaly)
          ),
          x = "Participant (ranked by score)",
          y = "Anomaly Score",
          color = "Classification"
        ) +
        theme_minimal(base_size = 12) +
        theme(
          plot.title = element_text(face = "bold"),
          legend.position = "bottom"
        )
    },

    #' @description Compare detection methods
    #' @param participant_features Data frame from FeatureEngineer
    #' @param features Feature names to use
    #' @return Data frame comparing methods
    compare_methods = function(participant_features, features = NULL) {
      methods <- c("isolation_forest", "mahalanobis")

      results <- lapply(methods, function(m) {
        result <- tryCatch(
          self$profile_participants(participant_features, features, method = m,
                                    compute_contributions = FALSE),
          error = function(e) NULL
        )
        if (!is.null(result)) {
          data.frame(
            method = m,
            n_flagged = sum(result$is_anomaly),
            threshold = result$threshold,
            score_mean = mean(result$anomaly_scores, na.rm = TRUE),
            score_sd = sd(result$anomaly_scores, na.rm = TRUE),
            stringsAsFactors = FALSE
          )
        }
      })

      do.call(rbind, results[!sapply(results, is.null)])
    },

    #' @description Plot SHAP-like feature contributions for a participant
    #' @param result ParticipantProfileResult object
    #' @param participant_id Participant ID to plot
    #' @param top_k Number of top features to show
    #' @return ggplot object
    plot_contributions = function(result, participant_id, top_k = 10) {
      stopifnot(inherits(result, "ParticipantProfileResult"))

      contrib_df <- result$get_top_contributors(participant_id, top_k)
      if (nrow(contrib_df) == 0) {
        stop("No contributions available for participant: ", participant_id)
      }

      contrib_df$feature <- factor(contrib_df$feature, levels = rev(contrib_df$feature))
      idx <- which(result$profiles$id == participant_id)
      score <- result$anomaly_scores[idx]

      ggplot(contrib_df, aes(x = .data$feature, y = .data$contribution,
                             fill = .data$contribution)) +
        geom_col(width = 0.7) +
        geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
        coord_flip() +
        scale_fill_gradient2(
          low = "#2E86AB", mid = "white", high = "#E63946",
          midpoint = 0, guide = "none"
        ) +
        labs(
          title = sprintf("Feature Contributions - %s", participant_id),
          subtitle = sprintf("Anomaly Score: %.3f | Threshold: %.3f",
                             score, result$threshold),
          x = NULL,
          y = "Contribution to Anomaly Score"
        ) +
        theme_minimal(base_size = 11) +
        theme(
          plot.title = element_text(face = "bold"),
          axis.text.y = element_text(size = 9)
        )
    },

    #' @description Plot aggregate feature importance across all anomalous participants
    #' @param result ParticipantProfileResult object
    #' @param top_k Number of top features to show
    #' @return ggplot object
    plot_aggregate_importance = function(result, top_k = 15) {
      stopifnot(inherits(result, "ParticipantProfileResult"))

      if (is.null(result$feature_contributions)) {
        stop("No feature contributions available")
      }

      # Aggregate across anomalous participants
      anomaly_idx <- which(result$is_anomaly)
      if (length(anomaly_idx) == 0) {
        stop("No anomalies to aggregate")
      }

      all_contrib <- lapply(anomaly_idx, function(i) {
        contrib <- result$feature_contributions[[i]]
        if (!is.null(contrib)) {
          data.frame(
            feature = names(contrib),
            abs_contribution = abs(contrib),
            stringsAsFactors = FALSE
          )
        }
      })

      agg <- bind_rows(all_contrib) |>
        group_by(.data$feature) |>
        summarize(
          mean_abs_contribution = mean(.data$abs_contribution, na.rm = TRUE),
          .groups = "drop"
        ) |>
        arrange(desc(.data$mean_abs_contribution)) |>
        slice_head(n = top_k)

      agg$feature <- factor(agg$feature, levels = rev(agg$feature))

      ggplot(agg, aes(x = .data$feature, y = .data$mean_abs_contribution)) +
        geom_col(fill = "#2E86AB", width = 0.7) +
        coord_flip() +
        labs(
          title = "Aggregate Feature Importance",
          subtitle = sprintf("Mean |contribution| across %d anomalous participants",
                             length(anomaly_idx)),
          x = NULL,
          y = "Mean Absolute Contribution"
        ) +
        theme_minimal(base_size = 11) +
        theme(plot.title = element_text(face = "bold"))
    },

    #' @description Plot all anomalous participants with their top contributors
    #' @param result ParticipantProfileResult object
    #' @param top_k Number of top features per participant
    #' @return ggplot object
    plot_all_anomaly_explanations = function(result, top_k = 3) {
      stopifnot(inherits(result, "ParticipantProfileResult"))

      anomalous_ids <- result$get_anomalous_ids()
      if (length(anomalous_ids) == 0) {
        stop("No anomalies to plot")
      }

      all_contrib <- lapply(anomalous_ids, function(id) {
        df <- result$get_top_contributors(id, top_k)
        if (nrow(df) > 0) {
          df$participant_id <- id
          idx <- which(result$profiles$id == id)
          df$anomaly_score <- result$anomaly_scores[idx]
          df
        }
      })

      plot_data <- bind_rows(all_contrib)
      if (nrow(plot_data) == 0) {
        stop("No contribution data available")
      }

      # Order participants by anomaly score
      participant_order <- plot_data |>
        group_by(.data$participant_id) |>
        summarize(score = unique(.data$anomaly_score), .groups = "drop") |>
        arrange(desc(.data$score))

      plot_data$participant_id <- factor(
        plot_data$participant_id,
        levels = participant_order$participant_id
      )

      ggplot(plot_data, aes(x = .data$feature, y = .data$contribution,
                            fill = .data$contribution)) +
        geom_col(width = 0.7) +
        geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.3) +
        coord_flip() +
        facet_wrap(~.data$participant_id, scales = "free_y", ncol = 2) +
        scale_fill_gradient2(
          low = "#2E86AB", mid = "white", high = "#E63946",
          midpoint = 0, guide = "none"
        ) +
        labs(
          title = "Feature Contributions for Anomalous Participants",
          subtitle = sprintf("Top %d features per participant | Threshold: %.3f",
                             top_k, result$threshold),
          x = NULL,
          y = "Contribution"
        ) +
        theme_minimal(base_size = 10) +
        theme(
          plot.title = element_text(face = "bold"),
          strip.text = element_text(face = "bold"),
          axis.text.y = element_text(size = 8),
          panel.spacing = unit(1, "lines")
        )
    }
  ),

  private = list(
    .random_state = NULL,
    .n_trees = NULL,

    #' Impute NA values with column medians
    .impute_nas = function(feature_matrix) {
      for (i in seq_len(ncol(feature_matrix))) {
        nas <- is.na(feature_matrix[, i])
        if (any(nas)) {
          med <- median(feature_matrix[!nas, i], na.rm = TRUE)
          if (is.na(med)) med <- 0  # Fallback if all NA
          feature_matrix[nas, i] <- med
        }
      }
      feature_matrix
    },

    #' Fit Isolation Forest and return both model and scores
    .isolation_forest_with_model = function(feature_matrix) {
      # Ensure at least 2 columns (isotree requirement)
      if (ncol(feature_matrix) == 1) {
        feature_matrix <- cbind(feature_matrix, feature_matrix)
      }

      model <- isolation.forest(
        feature_matrix,
        ntrees = private$.n_trees,
        seed = private$.random_state,
        nthreads = 1  # Reproducibility
      )

      scores <- predict(model, feature_matrix)
      scores <- pmin(pmax(scores, 0), 1)

      list(model = model, scores = scores)
    },

    #' Calculate Mahalanobis distance-based scores
    .mahalanobis_scores = function(feature_matrix) {
      n <- nrow(feature_matrix)
      p <- ncol(feature_matrix)

      # Robust center (median)
      center <- apply(feature_matrix, 2, median, na.rm = TRUE)

      # Covariance matrix
      cov_matrix <- tryCatch(
        cov(feature_matrix, use = "pairwise.complete.obs"),
        error = function(e) diag(p)
      )

      # Add small ridge to avoid singularity
      cov_matrix <- cov_matrix + 0.01 * diag(p)

      # Calculate Mahalanobis distances
      distances <- mahalanobis(feature_matrix, center, cov_matrix)

      # Normalize to [0, 1] using chi-squared CDF
      pchisq(distances, df = p)
    },

    #' Compute SHAP-like feature contributions via permutation
    #' For each feature, replace with median and measure score change
    .compute_contributions = function(feature_matrix, feature_names, iso_model, baseline_scores) {
      n_obs <- nrow(feature_matrix)
      n_features <- length(feature_names)

      # Compute medians for perturbation baseline
      medians <- apply(feature_matrix, 2, median, na.rm = TRUE)

      # For each observation, compute contribution of each feature
      contributions <- lapply(seq_len(n_obs), function(i) {
        private$.compute_single_contribution(
          i, feature_matrix, feature_names, iso_model, baseline_scores[i], medians
        )
      })

      contributions
    },

    #' Compute contributions for single observation
    .compute_single_contribution = function(idx, feature_matrix, feature_names,
                                             iso_model, baseline_score, medians) {
      n_features <- length(feature_names)
      contributions <- numeric(n_features)
      names(contributions) <- feature_names

      for (j in seq_len(n_features)) {
        # Create perturbed observation (replace feature j with median)
        perturbed <- feature_matrix[idx, , drop = FALSE]
        perturbed[1, j] <- medians[j]

        # Get perturbed score
        perturbed_score <- predict(iso_model, perturbed)

        # Contribution = baseline - perturbed
        # Positive: feature increases anomaly score
        # Negative: feature decreases anomaly score
        contributions[j] <- baseline_score - perturbed_score
      }

      contributions
    },

    #' Generate SHAP-based explanations for flagged participants
    .generate_shap_explanations = function(features_df, feature_names, is_anomaly, contributions) {
      explanations <- character(nrow(features_df))

      for (i in which(is_anomaly)) {
        if (!is.null(contributions) && !is.null(contributions[[i]])) {
          explanations[i] <- private$.format_shap_explanation(
            features_df$id[i], contributions[[i]], feature_names, features_df, i
          )
        } else {
          # Fallback to z-score if no contributions
          explanations[i] <- private$.format_zscore_explanation(
            features_df$id[i], feature_names, features_df, i
          )
        }
      }

      explanations
    },

    #' Format SHAP-based explanation string
    .format_shap_explanation = function(id, contributions, feature_names, features_df, idx) {
      sorted_idx <- order(abs(contributions), decreasing = TRUE)
      top_idx <- head(sorted_idx, 3)

      reasons <- sapply(top_idx, function(j) {
        contrib <- contributions[j]
        feature <- feature_names[j]
        direction <- if (contrib > 0) "increases" else "decreases"
        magnitude <- abs(contrib)
        sprintf("%s %s score by %.3f", feature, direction, magnitude)
      })

      paste(id, "flagged:", paste(reasons, collapse = "; "))
    },

    #' Format z-score fallback explanation
    .format_zscore_explanation = function(id, feature_names, features_df, idx) {
      feature_means <- sapply(feature_names, function(f) mean(features_df[[f]], na.rm = TRUE))
      feature_sds <- sapply(feature_names, function(f) {
        s <- sd(features_df[[f]], na.rm = TRUE)
        if (is.na(s) || s < 1e-10) 1 else s
      })

      z_scores <- sapply(seq_along(feature_names), function(j) {
        val <- as.numeric(features_df[[feature_names[j]]][idx])
        if (is.na(val)) return(0)
        (val - feature_means[j]) / feature_sds[j]
      })
      names(z_scores) <- feature_names

      sorted_idx <- order(abs(z_scores), decreasing = TRUE)
      top_idx <- head(sorted_idx, 3)

      reasons <- sapply(top_idx, function(j) {
        z <- z_scores[j]
        direction <- if (z > 0) "HIGH" else "LOW"
        sprintf("%s %s (z=%.1f)", feature_names[j], direction, z)
      })

      paste(id, "flagged:", paste(reasons, collapse = "; "))
    }
  )
)
