# R/calculators/anomaly_detector.R
# Service: Multivariate Anomaly Detection using Isolation Forest
#
# =============================================================================
# PRINCIPLES APPLIED (from CLAUDE.md)
# =============================================================================
#
# NAMING PRINCIPLES:
# - Understandability: Problem domain terms (anomaly, contamination, threshold)
# - Consistency: All detection methods start with "detect_", plots with "plot_"
# - Distinguishability: Specific names (raw_data vs cv_residual vs random_effects)
# - Conciseness: Clear but complete (contamination_rate not contamination_rate_proportion)
#
# SOLID PRINCIPLES:
# - SRP: Single responsibility - only detects anomalies
# - OCP: Open for extension via new detection stages
# - DIP: Depends on abstractions (data frames), not isotree internals
#
# CUPID PRINCIPLES:
# - Composable: Each detection method is independent
# - Unix: Each method does one thing well
# - Predictable: Same inputs + same seed -> same outputs
# - Idiomatic: Follows R6/R conventions
# - Domain-based: Names reflect statistical concepts
#
# SCIENTIFIC VALIDITY:
# - Isolation Forest is appropriate for multivariate anomaly detection
# - Contamination rate controls false positive rate
# - Three-stage pipeline catches different anomaly types
# =============================================================================

box::use(
  R6[R6Class],
  isotree[isolation.forest, predict.isolation_forest],
  ggplot2[
    ggplot, aes, geom_histogram, geom_vline, geom_point, geom_segment,
    geom_hline, labs, theme_minimal, theme, element_text, element_blank,
    scale_fill_manual, scale_color_manual, after_stat, coord_flip
  ],
  dplyr[mutate, arrange, group_by, summarize, ungroup, n, .data],
  stats[predict, setNames, quantile, sd]
)

#' Anomaly Result Value Object
#'
#' Stores results from anomaly detection on observations.
#'
#' @export
AnomalyResult <- R6Class(
  classname = "AnomalyResult",
  cloneable = FALSE,

  public = list(
    #' @field scores Numeric vector [0,1] per observation (higher = more anomalous)
    scores = NULL,
    #' @field is_anomaly Logical vector flagging anomalies
    is_anomaly = NULL,
    #' @field threshold Cutoff score used for flagging
    threshold = NULL,
    #' @field n_anomalies Count of flagged observations
    n_anomalies = NULL,
    #' @field anomaly_indices Which rows are anomalies
    anomaly_indices = NULL,
    #' @field contamination_rate Expected contamination rate used
    contamination_rate = NULL,

    #' @description Create AnomalyResult
    #' @param scores Numeric vector of anomaly scores
    #' @param threshold Cutoff for flagging
    #' @param contamination_rate Expected contamination rate
    initialize = function(scores, threshold, contamination_rate) {
      self$scores <- scores
      self$threshold <- threshold
      self$contamination_rate <- contamination_rate
      self$is_anomaly <- scores >= threshold
      self$n_anomalies <- sum(self$is_anomaly)
      self$anomaly_indices <- which(self$is_anomaly)
    },

    #' @description Summary of anomaly detection
    #' @return Named list with key metrics
    summarize = function() {
      list(
        n_total = length(self$scores),
        n_anomalies = self$n_anomalies,
        anomaly_rate = self$n_anomalies / length(self$scores),
        threshold = self$threshold,
        score_mean = mean(self$scores),
        score_sd = sd(self$scores),
        score_max = max(self$scores)
      )
    }
  )
)

#' Participant Anomaly Result Value Object
#'
#' Stores per-participant anomaly detection results.
#'
#' @export
ParticipantAnomalyResult <- R6Class(
  classname = "ParticipantAnomalyResult",
  cloneable = FALSE,

  public = list(
    #' @field participant_id Character vector of participant IDs
    participant_id = NULL,
    #' @field mean_score Mean anomaly score per participant
    mean_score = NULL,
    #' @field max_score Max anomaly score per participant
    max_score = NULL,
    #' @field n_anomalies Number of anomalous observations per participant
    n_anomalies = NULL,
    #' @field n_observations Total observations per participant
    n_observations = NULL,
    #' @field is_anomaly Logical per participant (based on aggregated scores)
    is_anomaly = NULL,
    #' @field interpretation "normal", "borderline", "anomalous"
    interpretation = NULL,
    #' @field threshold Cutoff used for participant-level flagging
    threshold = NULL,

    #' @description Create ParticipantAnomalyResult
    #' @param participant_summary Data frame with per-participant metrics
    #' @param threshold Cutoff for participant-level flagging
    initialize = function(participant_summary, threshold) {
      self$participant_id <- participant_summary$participant_id
      self$mean_score <- participant_summary$mean_score
      self$max_score <- participant_summary$max_score
      self$n_anomalies <- participant_summary$n_anomalies
      self$n_observations <- participant_summary$n_observations
      self$threshold <- threshold

      # Flag participants based on mean score
      self$is_anomaly <- self$mean_score >= threshold

      # Interpretation based on score bands
      self$interpretation <- ifelse(
        self$mean_score >= threshold, "anomalous",
        ifelse(self$mean_score >= threshold * 0.7, "borderline", "normal")
      )
    },

    #' @description Get anomalous participants
    #' @return Character vector of anomalous participant IDs
    get_anomalous_ids = function() {
      self$participant_id[self$is_anomaly]
    },

    #' @description Summary of participant-level detection
    #' @return Named list with key metrics
    summarize = function() {
      list(
        n_participants = length(self$participant_id),
        n_anomalous = sum(self$is_anomaly),
        n_borderline = sum(self$interpretation == "borderline"),
        n_normal = sum(self$interpretation == "normal"),
        anomalous_ids = self$get_anomalous_ids()
      )
    }
  )
)

#' Anomaly Detector
#'
#' R6 class for multivariate anomaly detection using Isolation Forest.
#' Implements three-stage pipeline: pre-modeling, CV residuals, random effects.
#'
#' @section Scientific Background:
#' Isolation Forest detects anomalies by isolating observations.
#' Anomalies are easier to isolate (shorter path length in tree).
#' Score normalized to [0,1] where higher = more anomalous.
#'
#' @export
AnomalyDetector <- R6Class(

  classname = "AnomalyDetector",
  cloneable = FALSE,

  public = list(
    #' @description Create a new AnomalyDetector instance
    #' @param random_state Seed for reproducibility (default: 42)
    #' @param n_trees Number of trees in forest (default: 100)
    initialize = function(random_state = 42, n_trees = 100) {
      private$.random_state <- random_state
      private$.n_trees <- n_trees
    },

    #' @description Stage 1: Pre-modeling raw data screening
    #'
    #' Detects multivariate anomalies in raw data before modeling.
    #' Useful for identifying measurement errors, data entry issues,

    #' or unusual combinations of features.
    #'
    #' @param data Data frame with numeric features
    #' @param features Character vector of column names to use
    #' @param contamination Expected proportion of anomalies (default: 0.05)
    #' @return AnomalyResult object
    detect_raw_data_anomalies = function(data, features, contamination = 0.05) {
      stopifnot(
        "data must be a data frame" = is.data.frame(data),
        "features must be character vector" = is.character(features),
        "all features must exist in data" = all(features %in% names(data)),
        "contamination must be between 0 and 0.5" = contamination > 0 && contamination <= 0.5
      )

      # Extract feature matrix
      feature_matrix <- as.matrix(data[, features, drop = FALSE])

      # Check for missing values
      if (any(is.na(feature_matrix))) {
        warning("Missing values detected. Rows with NA will have NA scores.")
      }

      # Fit Isolation Forest and get scores
      scores <- private$.fit_and_score(feature_matrix)

      # Calculate threshold based on contamination rate
      threshold <- private$.calculate_threshold(scores, contamination)

      AnomalyResult$new(scores, threshold, contamination)
    },

    #' @description Stage 2: CV residual anomaly detection
    #'
    #' Identifies participants whose predictions are consistently poor
    #' based on cross-validation residuals.
    #'
    #' @param residuals Numeric vector of CV residuals
    #' @param participant_ids Character/factor vector of participant IDs
    #' @param contamination Expected proportion of anomalous participants (default: 0.1)
    #' @return ParticipantAnomalyResult object
    detect_cv_residual_anomalies = function(residuals, participant_ids, contamination = 0.1) {
      stopifnot(
        "residuals must be numeric" = is.numeric(residuals),
        "participant_ids must match residuals length" = length(participant_ids) == length(residuals),
        "contamination must be between 0 and 0.5" = contamination > 0 && contamination <= 0.5
      )

      # Aggregate residuals by participant
      resid_df <- data.frame(
        participant_id = as.character(participant_ids),
        residual = residuals,
        abs_residual = abs(residuals)
      )

      participant_stats <- resid_df |>
        group_by(.data$participant_id) |>
        summarize(
          mean_abs_residual = mean(.data$abs_residual),
          max_abs_residual = max(.data$abs_residual),
          sd_residual = sd(.data$residual),
          n_observations = n(),
          .groups = "drop"
        )

      # Use residual statistics as features for Isolation Forest
      feature_matrix <- as.matrix(participant_stats[, c("mean_abs_residual", "max_abs_residual", "sd_residual")])

      # Fit and score
      scores <- private$.fit_and_score(feature_matrix)

      # Calculate threshold
      threshold <- private$.calculate_threshold(scores, contamination)

      # Build participant summary for result object
      participant_summary <- data.frame(
        participant_id = participant_stats$participant_id,
        mean_score = scores,
        max_score = scores,  # Using same score for simplicity
        n_anomalies = ifelse(scores >= threshold, 1, 0),
        n_observations = participant_stats$n_observations
      )

      ParticipantAnomalyResult$new(participant_summary, threshold)
    },

    #' @description Stage 3: Random effects anomaly detection
    #'
    #' Identifies participants with unusual random effects (intercepts/slopes).
    #' Useful for detecting athletes with atypical velocity-RIR relationships.
    #'
    #' @param model An lme4 lmerMod object
    #' @param contamination Expected proportion of anomalous participants (default: 0.1)
    #' @return ParticipantAnomalyResult object
    detect_random_effects_anomalies = function(model, contamination = 0.1) {
      stopifnot(
        "model must be an lmerMod object" = inherits(model, "lmerMod"),
        "contamination must be between 0 and 0.5" = contamination > 0 && contamination <= 0.5
      )

      # Extract random effects
      re <- lme4::ranef(model)[[1]]
      re_df <- as.data.frame(re)
      re_df$participant_id <- rownames(re)

      # Use random effects as features
      feature_cols <- setdiff(names(re_df), "participant_id")
      feature_matrix <- as.matrix(re_df[, feature_cols, drop = FALSE])

      # Fit and score
      scores <- private$.fit_and_score(feature_matrix)

      # Calculate threshold
      threshold <- private$.calculate_threshold(scores, contamination)

      # Build participant summary
      participant_summary <- data.frame(
        participant_id = re_df$participant_id,
        mean_score = scores,
        max_score = scores,
        n_anomalies = ifelse(scores >= threshold, 1, 0),
        n_observations = NA_integer_  # Not applicable for random effects
      )

      ParticipantAnomalyResult$new(participant_summary, threshold)
    },

    #' @description Compare Isolation Forest results with influence diagnostics
    #'
    #' Shows agreement/disagreement between IF anomaly detection
    #' and traditional influence measures (Cook's D).
    #'
    #' @param anomaly_result AnomalyResult or ParticipantAnomalyResult
    #' @param influence_flags Logical vector from Cook's D (TRUE = influential)
    #' @return List with agreement metrics
    compare_with_influence_diagnostics = function(anomaly_result, influence_flags) {
      stopifnot(
        "influence_flags must be logical" = is.logical(influence_flags)
      )

      if (inherits(anomaly_result, "AnomalyResult")) {
        if_flags <- anomaly_result$is_anomaly
      } else if (inherits(anomaly_result, "ParticipantAnomalyResult")) {
        if_flags <- anomaly_result$is_anomaly
      } else {
        stop("anomaly_result must be AnomalyResult or ParticipantAnomalyResult")
      }

      stopifnot(
        "lengths must match" = length(if_flags) == length(influence_flags)
      )

      # Calculate agreement metrics
      both_flagged <- sum(if_flags & influence_flags)
      only_if <- sum(if_flags & !influence_flags)
      only_cooks <- sum(!if_flags & influence_flags)
      neither <- sum(!if_flags & !influence_flags)

      total <- length(if_flags)
      agreement_rate <- (both_flagged + neither) / total

      list(
        both_flagged = both_flagged,
        only_isolation_forest = only_if,
        only_cooks_distance = only_cooks,
        neither_flagged = neither,
        agreement_rate = agreement_rate,
        jaccard_index = if (both_flagged + only_if + only_cooks > 0) {
          both_flagged / (both_flagged + only_if + only_cooks)
        } else {
          1.0  # Perfect agreement if neither method flags anything
        },
        interpretation = private$.interpret_comparison(both_flagged, only_if, only_cooks)
      )
    },

    #' @description Plot anomaly score distribution
    #'
    #' @param anomaly_result AnomalyResult object
    #' @param title Plot title
    #' @return ggplot object
    plot_anomaly_scores = function(anomaly_result, title = "Anomaly Score Distribution") {
      stopifnot(
        "anomaly_result must be AnomalyResult" = inherits(anomaly_result, "AnomalyResult")
      )

      df <- data.frame(
        score = anomaly_result$scores,
        is_anomaly = factor(
          ifelse(anomaly_result$is_anomaly, "Anomaly", "Normal"),
          levels = c("Normal", "Anomaly")
        )
      )

      ggplot(df, aes(x = .data$score, fill = .data$is_anomaly)) +
        geom_histogram(bins = 30, alpha = 0.7, position = "identity") +
        geom_vline(
          xintercept = anomaly_result$threshold,
          linetype = "dashed",
          color = "red",
          linewidth = 1
        ) +
        scale_fill_manual(values = c("Normal" = "#2E86AB", "Anomaly" = "#E63946")) +
        labs(
          title = title,
          subtitle = sprintf(
            "Threshold = %.3f | %d anomalies flagged (%.1f%%)",
            anomaly_result$threshold,
            anomaly_result$n_anomalies,
            100 * anomaly_result$n_anomalies / length(anomaly_result$scores)
          ),
          x = "Anomaly Score (higher = more anomalous)",
          y = "Count",
          fill = "Classification"
        ) +
        theme_minimal(base_size = 12) +
        theme(
          plot.title = element_text(face = "bold"),
          legend.position = "bottom"
        )
    },

    #' @description Plot participant-level anomaly results
    #'
    #' @param participant_result ParticipantAnomalyResult object
    #' @param title Plot title
    #' @return ggplot object
    plot_participant_anomalies = function(participant_result, title = "Participant Anomaly Scores") {
      stopifnot(
        "participant_result must be ParticipantAnomalyResult" =
          inherits(participant_result, "ParticipantAnomalyResult")
      )

      df <- data.frame(
        participant_id = participant_result$participant_id,
        score = participant_result$mean_score,
        interpretation = factor(
          participant_result$interpretation,
          levels = c("normal", "borderline", "anomalous")
        )
      )

      # Order by score
      df <- df[order(df$score), ]
      df$order <- seq_len(nrow(df))

      ggplot(df, aes(x = .data$score, y = .data$order, color = .data$interpretation)) +
        geom_vline(
          xintercept = participant_result$threshold,
          linetype = "dashed",
          color = "gray50"
        ) +
        geom_point(size = 3) +
        scale_color_manual(values = c(
          "normal" = "#2E86AB",
          "borderline" = "#F77F00",
          "anomalous" = "#E63946"
        )) +
        labs(
          title = title,
          subtitle = sprintf(
            "Threshold = %.3f | %d anomalous, %d borderline",
            participant_result$threshold,
            sum(participant_result$is_anomaly),
            sum(participant_result$interpretation == "borderline")
          ),
          x = "Anomaly Score",
          y = "Participant (ordered by score)",
          color = "Classification"
        ) +
        theme_minimal(base_size = 12) +
        theme(
          plot.title = element_text(face = "bold"),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          legend.position = "bottom"
        )
    }
  ),

  private = list(
    .random_state = NULL,
    .n_trees = NULL,

    #' @description Fit Isolation Forest and return anomaly scores
    #' @param feature_matrix Numeric matrix
    #' @return Numeric vector of scores [0, 1]
    .fit_and_score = function(feature_matrix) {
      # Handle single-column case
      if (ncol(feature_matrix) == 1) {
        feature_matrix <- cbind(feature_matrix, feature_matrix)
      }

      # Fit Isolation Forest
      model <- isolation.forest(
        feature_matrix,
        ntrees = private$.n_trees,
        seed = private$.random_state,
        nthreads = 1  # Reproducibility
      )

      # Get anomaly scores (depth-based, normalized)
      scores <- predict(model, feature_matrix)

      # Ensure scores are in [0, 1]
      scores <- pmin(pmax(scores, 0), 1)

      scores
    },

    #' @description Calculate threshold based on contamination rate
    #' @param scores Numeric vector of anomaly scores
    #' @param contamination Expected contamination rate
    #' @return Numeric threshold value
    .calculate_threshold = function(scores, contamination) {
      # Threshold is the (1 - contamination) quantile
      # Higher scores = more anomalous, so we want top `contamination` fraction
      quantile(scores, probs = 1 - contamination, na.rm = TRUE)
    },

    #' @description Interpret comparison between methods
    .interpret_comparison = function(both, only_if, only_cooks) {
      total_flagged <- both + only_if + only_cooks

      if (total_flagged == 0) {
        return("No anomalies detected by either method")
      }

      if (both > 0 && only_if == 0 && only_cooks == 0) {
        return("Perfect agreement: both methods flag the same observations")
      }

      if (only_if > only_cooks) {
        return("Isolation Forest flags more observations (may detect multivariate patterns)")
      }

      if (only_cooks > only_if) {
        return("Cook's distance flags more observations (may be more sensitive to leverage)")
      }

      "Methods show partial agreement"
    }
  )
)

#' Convenience function for raw data anomaly detection
#'
#' @param data Data frame with numeric features
#' @param features Character vector of column names to use
#' @param contamination Expected proportion of anomalies (default: 0.05)
#' @param random_state Seed for reproducibility
#' @return AnomalyResult object
#' @export
detect_raw_anomalies <- function(data, features, contamination = 0.05, random_state = 42) {
  detector <- AnomalyDetector$new(random_state = random_state)
  detector$detect_raw_data_anomalies(data, features, contamination)
}

#' Convenience function for random effects anomaly detection
#'
#' @param model An lme4 lmerMod object
#' @param contamination Expected proportion of anomalous participants
#' @param random_state Seed for reproducibility
#' @return ParticipantAnomalyResult object
#' @export
detect_re_anomalies <- function(model, contamination = 0.1, random_state = 42) {
  detector <- AnomalyDetector$new(random_state = random_state)
  detector$detect_random_effects_anomalies(model, contamination)
}
