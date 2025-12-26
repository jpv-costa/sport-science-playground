# R/calculators/correlation_analyzer.R
# Service: Correlation Analysis
#
# SOLID Principles Applied:
# - SRP: Single responsibility - calculate correlations (actor: statistician)
# - OCP: Open for extension with different correlation methods

box::use(
  R6[R6Class]
)

#' Correlation Result
#'
#' Value object containing a correlation result.
#'
#' @export
CorrelationResult <- R6Class(

  classname = "CorrelationResult",
  cloneable = FALSE,

  public = list(

    #' @description Create correlation result
    #' @param r Correlation coefficient
    #' @param n Sample size
    #' @param group_id Optional group identifier
    initialize = function(r, n, group_id = NULL) {
      private$.r <- r
      private$.n <- n
      private$.group_id <- group_id
      private$.r_squared <- r^2
    },

    #' @description Convert to list
    to_list = function() {
      list(
        r = private$.r,
        r_squared = private$.r_squared,
        n = private$.n,
        group_id = private$.group_id
      )
    }
  ),

  active = list(
    r = function() private$.r,
    r_squared = function() private$.r_squared,
    n = function() private$.n,
    group_id = function() private$.group_id
  ),

  private = list(
    .r = NULL,
    .r_squared = NULL,
    .n = NULL,
    .group_id = NULL
  )
)


#' Correlation Summary
#'
#' Value object containing summary statistics for multiple correlations.
#'
#' @export
CorrelationSummary <- R6Class(

  classname = "CorrelationSummary",
  cloneable = FALSE,

  public = list(

    #' @description Create correlation summary
    #' @param correlations Vector of correlation coefficients
    #' @param group_ids Optional vector of group identifiers
    initialize = function(correlations, group_ids = NULL) {
      private$.correlations <- correlations
      private$.group_ids <- group_ids
      private$.n_groups <- length(correlations)
      private$.mean_r <- mean(correlations, na.rm = TRUE)
      private$.sd_r <- stats::sd(correlations, na.rm = TRUE)
      private$.min_r <- min(correlations, na.rm = TRUE)
      private$.max_r <- max(correlations, na.rm = TRUE)
      private$.mean_r_squared <- mean(correlations^2, na.rm = TRUE)
    },

    #' @description Convert to list
    to_list = function() {
      list(
        n_groups = private$.n_groups,
        mean_r = private$.mean_r,
        sd_r = private$.sd_r,
        min_r = private$.min_r,
        max_r = private$.max_r,
        mean_r_squared = private$.mean_r_squared
      )
    }
  ),

  active = list(
    n_groups = function() private$.n_groups,
    mean_r = function() private$.mean_r,
    sd_r = function() private$.sd_r,
    min_r = function() private$.min_r,
    max_r = function() private$.max_r,
    mean_r_squared = function() private$.mean_r_squared,
    correlations = function() private$.correlations
  ),

  private = list(
    .correlations = NULL,
    .group_ids = NULL,
    .n_groups = NULL,
    .mean_r = NULL,
    .sd_r = NULL,
    .min_r = NULL,
    .max_r = NULL,
    .mean_r_squared = NULL
  )
)


#' Correlation Analyzer
#'
#' Analyzes correlations between velocity and RIR.
#' Supports overall, per-participant, and per-exercise analyses.
#'
#' @export
CorrelationAnalyzer <- R6Class(

  classname = "CorrelationAnalyzer",

  public = list(

    #' @description Create analyzer
    #' @param method Correlation method ("pearson", "spearman", "kendall")
    initialize = function(method = "pearson") {
      private$.method <- method
    },

    #' @description Calculate overall correlation
    #' @param x First variable
    #' @param y Second variable
    #' @return CorrelationResult object
    calculate_overall = function(x, y) {
      r <- stats::cor(x, y, use = "complete.obs", method = private$.method)
      n <- sum(stats::complete.cases(x, y))

      CorrelationResult$new(r = r, n = n)
    },

    #' @description Calculate correlations by group
    #' @param data Data frame with x, y, and group columns
    #' @param x_col Name of x column
    #' @param y_col Name of y column
    #' @param group_col Name of grouping column
    #' @return CorrelationSummary object
    calculate_by_group = function(data, x_col, y_col, group_col) {
      groups <- unique(data[[group_col]])
      correlations <- numeric(length(groups))

      for (i in seq_along(groups)) {
        subset_data <- data[data[[group_col]] == groups[i], ]
        correlations[i] <- stats::cor(
          subset_data[[x_col]],
          subset_data[[y_col]],
          use = "complete.obs",
          method = private$.method
        )
      }

      CorrelationSummary$new(correlations = correlations, group_ids = groups)
    },

    #' @description Calculate correlations for multiple subgroups
    #' @param data Data frame
    #' @param x_col Name of x column
    #' @param y_col Name of y column
    #' @param subgroup_col Name of subgroup column (e.g., "exercise")
    #' @return Named list of CorrelationResult objects
    calculate_by_subgroup = function(data, x_col, y_col, subgroup_col) {
      subgroups <- unique(data[[subgroup_col]])
      results <- list()

      for (subgroup in subgroups) {
        subset_data <- data[data[[subgroup_col]] == subgroup, ]
        r <- stats::cor(
          subset_data[[x_col]],
          subset_data[[y_col]],
          use = "complete.obs",
          method = private$.method
        )
        n <- sum(stats::complete.cases(subset_data[[x_col]], subset_data[[y_col]]))
        results[[subgroup]] <- CorrelationResult$new(r = r, n = n, group_id = subgroup)
      }

      results
    },

    #' @description Full analysis with all metrics
    #' @param data Data frame
    #' @param velocity_col Velocity column name
    #' @param rir_col RIR column name
    #' @param participant_col Participant ID column name
    #' @param exercise_col Exercise column name
    #' @return Named list with all analysis results
    analyze_velocity_rir = function(data,
                                    velocity_col = "mean_velocity",
                                    rir_col = "perceived_rir",
                                    participant_col = "id",
                                    exercise_col = "exercise") {

      overall <- self$calculate_overall(data[[velocity_col]], data[[rir_col]])

      by_participant <- self$calculate_by_group(
        data, velocity_col, rir_col, participant_col
      )

      by_exercise <- self$calculate_by_subgroup(
        data, velocity_col, rir_col, exercise_col
      )

      list(
        overall = overall,
        by_participant = by_participant,
        by_exercise = by_exercise
      )
    }
  ),

  private = list(
    .method = NULL
  )
)
