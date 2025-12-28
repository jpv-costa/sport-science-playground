# R/calculators/threshold_selector.R
# Service: Configurable Threshold Selection for Anomaly Detection
#
# =============================================================================
# PRINCIPLES APPLIED (from CLAUDE.md)
# =============================================================================
#
# NAMING PRINCIPLES:
# - Understandability: Strategy names describe their statistical approach
# - Consistency: All strategies follow same interface (scores -> threshold)
# - Distinguishability: Clear differences between strategy names
# - Conciseness: Short but meaningful strategy identifiers
#
# SOLID PRINCIPLES:
# - SRP: Single responsibility - only selects thresholds
# - OCP: Open for extension with new strategies via switch
# - DIP: Depends on abstractions (numeric vectors), not implementations
#
# CUPID PRINCIPLES:
# - Composable: Can be used with any anomaly score vector
# - Unix: Each method does one thing well
# - Predictable: Same inputs -> same outputs (deterministic)
# - Idiomatic: Follows R conventions, proper roxygen docs
# - Domain-based: Names reflect statistical concepts
#
# =============================================================================

box::use(
  R6[R6Class],
  stats[quantile, median, mad, kmeans, pchisq],
  ggplot2[
    ggplot, aes, geom_point, geom_hline, geom_text, labs,
    theme_minimal, theme, element_text, scale_color_manual, .data
  ]
)

#' Threshold Selection Strategy Constants
#'
#' Available strategies for selecting anomaly score thresholds.
#'
#' @export
ThresholdStrategy <- list(

FIXED_PERCENTILE = "fixed_percentile",
  ELBOW_DETECTION = "elbow_detection",
  GAP_STATISTICS = "gap_statistics",
  KNEE_POINT = "knee_point",
  MAD_BASED = "mad_based"
)

#' Threshold Selection Result Value Object
#'
#' Stores threshold selection results with diagnostics.
#'
#' @export
ThresholdResult <- R6Class(
  classname = "ThresholdResult",
  cloneable = FALSE,

  public = list(
    #' @field threshold Selected threshold value
    threshold = NULL,
    #' @field strategy Strategy used for selection
    strategy = NULL,
    #' @field diagnostics List of diagnostic information
    diagnostics = NULL,

    #' @description Create ThresholdResult
    #' @param threshold Numeric threshold value
    #' @param strategy Strategy name used
    #' @param diagnostics List of diagnostic metrics
    initialize = function(threshold, strategy, diagnostics) {
      self$threshold <- threshold
      self$strategy <- strategy
      self$diagnostics <- diagnostics
    },

    #' @description Summary of threshold selection
    #' @return Named list with key metrics
    summarize = function() {
      list(
        threshold = self$threshold,
        strategy = self$strategy,
        n_flagged = self$diagnostics$n_flagged,
        pct_flagged = self$diagnostics$pct_flagged
      )
    }
  )
)

#' Threshold Selector
#'
#' Selects optimal anomaly score threshold using configurable strategies.
#' Supports both fixed percentile and data-driven approaches.
#'
#' @section Strategies:
#' - FIXED_PERCENTILE: Traditional quantile-based (default: 95th percentile)
#' - ELBOW_DETECTION: Maximum perpendicular distance to diagonal
#' - GAP_STATISTICS: Largest gap in sorted scores
#' - KNEE_POINT: Kneedle algorithm for automatic knee detection
#' - MAD_BASED: Median + k * MAD (robust to outliers)
#'
#' @export
ThresholdSelector <- R6Class(

  classname = "ThresholdSelector",
  cloneable = FALSE,

  public = list(

    #' @description Select threshold using specified strategy
    #' @param scores Numeric vector of anomaly scores
    #' @param strategy One of ThresholdStrategy values (default: FIXED_PERCENTILE)
    #' @param params Named list of strategy-specific parameters:
    #'   - For FIXED_PERCENTILE: percentile (default: 0.95)
    #'   - For MAD_BASED: multiplier (default: 3)
    #' @return ThresholdResult object
    select_threshold = function(scores,
                                 strategy = ThresholdStrategy$FIXED_PERCENTILE,
                                 params = list()) {
      stopifnot(
        "scores must be numeric" = is.numeric(scores),
        "scores must have length > 1" = length(scores) > 1,
        "strategy must be valid" = strategy %in% unlist(ThresholdStrategy)
      )

      # Remove NAs for threshold calculation
      valid_scores <- scores[!is.na(scores)]

      if (length(valid_scores) < 2) {
        stop("Need at least 2 non-NA scores for threshold selection")
      }

      # Dispatch to appropriate strategy
      threshold <- switch(
        strategy,
        "fixed_percentile" = private$.percentile_threshold(valid_scores, params),
        "elbow_detection" = private$.elbow_threshold(valid_scores, params),
        "gap_statistics" = private$.gap_threshold(valid_scores, params),
        "knee_point" = private$.knee_point_threshold(valid_scores, params),
        "mad_based" = private$.mad_threshold(valid_scores, params),
        stop("Unknown threshold strategy: ", strategy)
      )

      # Build diagnostics
      diagnostics <- private$.build_diagnostics(scores, threshold, strategy, params)

      ThresholdResult$new(threshold, strategy, diagnostics)
    },

    #' @description Compare multiple threshold strategies
    #' @param scores Numeric vector of anomaly scores
    #' @return Data frame comparing all strategies
    compare_strategies = function(scores) {
      stopifnot(
        "scores must be numeric" = is.numeric(scores),
        "scores must have length > 1" = length(scores) > 1
      )

      strategies <- c(
        ThresholdStrategy$FIXED_PERCENTILE,
        ThresholdStrategy$ELBOW_DETECTION,
        ThresholdStrategy$GAP_STATISTICS,
        ThresholdStrategy$KNEE_POINT,
        ThresholdStrategy$MAD_BASED
      )

      results <- lapply(strategies, function(s) {
        result <- tryCatch(
          self$select_threshold(scores, s),
          error = function(e) NULL
        )
        if (!is.null(result)) {
          data.frame(
            strategy = s,
            threshold = result$threshold,
            n_flagged = result$diagnostics$n_flagged,
            pct_flagged = result$diagnostics$pct_flagged,
            stringsAsFactors = FALSE
          )
        }
      })

      do.call(rbind, results[!sapply(results, is.null)])
    },

    #' @description Plot score distribution with multiple thresholds
    #' @param scores Numeric vector of anomaly scores
    #' @param title Plot title
    #' @return ggplot object showing all strategy thresholds
    plot_threshold_comparison = function(scores, title = "Threshold Strategy Comparison") {
      comparison <- self$compare_strategies(scores)

      # Create data frame for plotting
      df <- data.frame(
        score = scores,
        index = seq_along(scores)
      )
      df <- df[order(df$score), ]
      df$rank <- seq_len(nrow(df))

      p <- ggplot(df, aes(x = .data$rank, y = .data$score)) +
        geom_point(alpha = 0.5, size = 2, color = "#2E86AB")

      # Add threshold lines for each strategy
      colors <- c(
        "fixed_percentile" = "#E63946",
        "elbow_detection" = "#F77F00",
        "gap_statistics" = "#2A9D8F",
        "knee_point" = "#9B59B6",
        "mad_based" = "#1ABC9C"
      )

      for (i in seq_len(nrow(comparison))) {
        strategy <- comparison$strategy[i]
        threshold <- comparison$threshold[i]
        p <- p + geom_hline(
          yintercept = threshold,
          linetype = "dashed",
          color = colors[strategy],
          linewidth = 0.8
        )
      }

      # Add labels
      label_df <- data.frame(
        strategy = comparison$strategy,
        threshold = comparison$threshold,
        label = sprintf("%s: %.3f (%d)", comparison$strategy,
                        comparison$threshold, comparison$n_flagged),
        x = max(df$rank) * 0.02
      )

      p <- p +
        geom_text(
          data = label_df,
          aes(x = .data$x, y = .data$threshold, label = .data$label),
          hjust = 0, vjust = -0.5, size = 3
        ) +
        labs(
          title = title,
          subtitle = "Dashed lines show thresholds from different strategies",
          x = "Observation (ranked by score)",
          y = "Anomaly Score"
        ) +
        theme_minimal(base_size = 12) +
        theme(plot.title = element_text(face = "bold"))

      p
    }
  ),

  private = list(

    #' Fixed percentile threshold
    .percentile_threshold = function(scores, params) {
      percentile <- params$percentile %||% 0.95
      stopifnot(
        "percentile must be between 0 and 1" = percentile > 0 && percentile < 1
      )
      as.numeric(quantile(scores, probs = percentile, na.rm = TRUE))
    },

    #' Elbow detection threshold
    #' Finds point of maximum perpendicular distance from diagonal
    .elbow_threshold = function(scores, params) {
      # Sort scores descending (highest first)
      sorted <- sort(scores, decreasing = TRUE)
      n <- length(sorted)

      if (n < 3) {
        return(median(sorted))
      }

      # Create normalized coordinates
      x <- seq_len(n)
      y <- sorted

      # Line from (1, y[1]) to (n, y[n])
      line_start <- c(1, y[1])
      line_end <- c(n, y[n])

      # Calculate perpendicular distance from each point to the line
      distances <- sapply(seq_len(n), function(i) {
        private$.point_to_line_distance(c(i, y[i]), line_start, line_end)
      })

      # Elbow is point of maximum distance
      elbow_idx <- which.max(distances)
      sorted[elbow_idx]
    },

    #' Gap statistics threshold
    #' Finds largest gap in sorted scores
    .gap_threshold = function(scores, params) {
      sorted <- sort(scores)
      n <- length(sorted)

      if (n < 3) {
        return(median(sorted))
      }

      # Calculate gaps between consecutive scores
      gaps <- diff(sorted)

      # Find index of largest gap
      max_gap_idx <- which.max(gaps)

      # Threshold is midpoint of largest gap
      (sorted[max_gap_idx] + sorted[max_gap_idx + 1]) / 2
    },

    #' Knee point threshold (Kneedle algorithm variant)
    .knee_point_threshold = function(scores, params) {
      # Sort scores descending
      sorted <- sort(scores, decreasing = TRUE)
      n <- length(sorted)

      if (n < 3) {
        return(median(sorted))
      }

      # Normalize to [0, 1]
      x_norm <- (seq_len(n) - 1) / (n - 1)

      score_range <- max(sorted) - min(sorted)
      if (score_range < 1e-10) {
        return(median(sorted))
      }
      y_norm <- (sorted - min(sorted)) / score_range

      # Compute difference curve (distance from diagonal)
      diff_curve <- y_norm - x_norm

      # Knee is maximum of difference curve
      knee_idx <- which.max(diff_curve)
      sorted[knee_idx]
    },

    #' MAD-based threshold
    #' Uses Median Absolute Deviation for robust threshold
    .mad_threshold = function(scores, params) {
      multiplier <- params$multiplier %||% 3
      stopifnot(
        "multiplier must be positive" = multiplier > 0
      )

      med <- median(scores, na.rm = TRUE)
      mad_val <- mad(scores, na.rm = TRUE)

      # If MAD is 0 (all same values), use IQR-based estimate
      if (mad_val < 1e-10) {
        iqr_val <- diff(quantile(scores, c(0.25, 0.75), na.rm = TRUE))
        mad_val <- iqr_val / 1.349  # Convert IQR to MAD equivalent
      }

      med + multiplier * mad_val
    },

    #' Calculate perpendicular distance from point to line
    .point_to_line_distance = function(point, line_start, line_end) {
      dx <- line_end[1] - line_start[1]
      dy <- line_end[2] - line_start[2]

      if (abs(dx) < 1e-10 && abs(dy) < 1e-10) {
        # Line has no length, return Euclidean distance to start
        return(sqrt(sum((point - line_start)^2)))
      }

      numerator <- abs(dy * point[1] - dx * point[2] +
                       line_end[1] * line_start[2] -
                       line_end[2] * line_start[1])
      denominator <- sqrt(dx^2 + dy^2)

      numerator / denominator
    },

    #' Build diagnostics for threshold result
    .build_diagnostics = function(scores, threshold, strategy, params) {
      n_total <- length(scores)
      n_flagged <- sum(scores >= threshold, na.rm = TRUE)

      list(
        n_total = n_total,
        n_flagged = n_flagged,
        pct_flagged = 100 * n_flagged / n_total,
        score_range = range(scores, na.rm = TRUE),
        score_median = median(scores, na.rm = TRUE),
        score_mad = mad(scores, na.rm = TRUE),
        params_used = params
      )
    }
  )
)

#' Null-coalescing operator
#' @keywords internal
`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}
