# R/domain/effect_size.R
# Domain Entity: EffectSize
#
# SOLID Principles Applied:
# - SRP: Represents a calculated effect size result (actor: meta-analyst)
# - OCP: Different effect size types use same interface

box::use(
  R6[R6Class],
  stats[qnorm]
)

#' Effect Size Entity
#'
#' Immutable value object representing a calculated effect size.
#' Contains the effect estimate, variance, and metadata about calculation.
#'
#' @description
#' Domain object encapsulating effect size calculation results.
#' Follows value object pattern - created once, never modified.
#'
#' @export
EffectSize <- R6Class(

  classname = "EffectSize",
  cloneable = FALSE,  # Immutable value object

  public = list(

    # =========================================================================
    # Constructor
    # =========================================================================

    #' @description Create a new EffectSize value object
    #' @param effect_estimate Standardized effect size (e.g., Hedge's g)
    #' @param sampling_variance Variance of the effect size
    #' @param sample_size Sample size used in calculation
    #' @param effect_type Type of effect size ("smcr", "romc")
    #' @param pre_post_correlation Correlation used in calculation
    #' @param bias_correction Whether small sample correction applied
    initialize = function(effect_estimate,
                          sampling_variance,
                          sample_size,
                          effect_type,
                          pre_post_correlation,
                          bias_correction = TRUE) {

      private$.effect_estimate <- as.numeric(effect_estimate)
      private$.sampling_variance <- as.numeric(sampling_variance)
      private$.sample_size <- as.integer(sample_size)
      private$.effect_type <- effect_type
      private$.pre_post_correlation <- as.numeric(pre_post_correlation)
      private$.bias_correction <- bias_correction
      private$.standard_error <- sqrt(sampling_variance)
    },

    # =========================================================================
    # Public Methods - Queries
    # =========================================================================

    #' @description Calculate confidence interval
    #' @param confidence_level Confidence level (default: 0.95)
    #' @return Named list with lower and upper bounds
    calculate_confidence_interval = function(confidence_level = 0.95) {
      z_critical <- qnorm(1 - (1 - confidence_level) / 2)
      margin <- z_critical * private$.standard_error

      list(
        lower = private$.effect_estimate - margin,
        upper = private$.effect_estimate + margin,
        confidence_level = confidence_level
      )
    },

    #' @description Calculate weight for meta-analysis (inverse variance)
    #' @return Numeric weight value
    calculate_weight = function() {
      1 / private$.sampling_variance
    },

    #' @description Check if effect size is statistically significant

    #' @param alpha Significance level (default: 0.05)
    #' @return Logical
    is_statistically_significant = function(alpha = 0.05) {
      ci <- self$calculate_confidence_interval(1 - alpha)
      # Significant if CI does not include zero
      !(ci$lower <= 0 && ci$upper >= 0)
    },

    #' @description Convert to list for serialization/data frame row
    #' @return Named list representation
    to_list = function() {
      ci <- self$calculate_confidence_interval()

      list(
        yi = private$.effect_estimate,       # metafor convention
        vi = private$.sampling_variance,     # metafor convention
        sei = private$.standard_error,
        ni = private$.sample_size,
        effect_type = private$.effect_type,
        correlation = private$.pre_post_correlation,
        bias_corrected = private$.bias_correction,
        ci_lower = ci$lower,
        ci_upper = ci$upper,
        weight = self$calculate_weight()
      )
    }
  ),

  # ===========================================================================
  # Active Bindings - Read-only access
  # ===========================================================================

  active = list(

    effect_estimate = function() private$.effect_estimate,

    sampling_variance = function() private$.sampling_variance,

    standard_error = function() private$.standard_error,

    sample_size = function() private$.sample_size,

    effect_type = function() private$.effect_type,

    pre_post_correlation = function() private$.pre_post_correlation,

    is_bias_corrected = function() private$.bias_correction,

    # metafor-compatible accessors
    yi = function() private$.effect_estimate,
    vi = function() private$.sampling_variance
  ),

  # ===========================================================================
  # Private
  # ===========================================================================

  private = list(
    .effect_estimate = NULL,
    .sampling_variance = NULL,
    .standard_error = NULL,
    .sample_size = NULL,
    .effect_type = NULL,
    .pre_post_correlation = NULL,
    .bias_correction = NULL
  )
)
