# R/calculators/effect_size_calculator.R
# Service: Effect Size Calculation
#
# SOLID Principles Applied:
# - SRP: Single responsibility - calculate effect sizes (actor: statistician)
# - OCP: Open for extension via strategy pattern for different ES types
# - DIP: Depends on EffectSize abstraction, not concrete implementation
# - ISP: Minimal interface - one public method

box::use(
  R6[R6Class],
  ../domain/effect_size[EffectSize],
  ../domain/treatment_group[TreatmentGroup]
)

# lgamma is a base function, always available

#' SMCR Effect Size Calculator
#'
#' Calculates Standardized Mean Change using Raw score standardization.
#' Applies Hedge's g correction for small sample bias.
#'
#' @description
#' Service class for computing SMCR effect sizes from TreatmentGroup data.
#' Named after domain concept: "Standardized Mean Change, Raw score".
#'
#' @export
SMCRCalculator <- R6Class(

  classname = "SMCRCalculator",

  public = list(

    # =========================================================================
    # Constructor - Dependency Injection
    # =========================================================================

    #' @description Create calculator with configuration
    #' @param default_correlation Default pre-post correlation when not provided
    #' @param apply_bias_correction Apply Hedge's g small sample correction
    initialize = function(default_correlation = 0.5,
                          apply_bias_correction = TRUE) {

      private$.default_correlation <- default_correlation
      private$.apply_bias_correction <- apply_bias_correction
    },

    # =========================================================================
    # Public Methods - Commands
    # =========================================================================

    #' @description Calculate effect size from a TreatmentGroup
    #' @param treatment_group TreatmentGroup object with measurements
    #' @return EffectSize object
    calculate = function(treatment_group) {
      private$.validate_input(treatment_group)

      correlation <- private$.get_correlation(treatment_group)
      raw_effect <- private$.compute_raw_effect(treatment_group)
      corrected_effect <- private$.apply_hedges_correction(
        raw_effect,
        treatment_group$sample_size
      )
      variance <- private$.compute_variance(
        corrected_effect,
        treatment_group$sample_size,
        correlation
      )

      EffectSize$new(
        effect_estimate = corrected_effect,
        sampling_variance = variance,
        sample_size = treatment_group$sample_size,
        effect_type = "smcr",
        pre_post_correlation = correlation,
        bias_correction = private$.apply_bias_correction
      )
    },

    #' @description Calculate effect sizes for multiple groups
    #' @param treatment_groups List of TreatmentGroup objects
    #' @return List of EffectSize objects
    calculate_batch = function(treatment_groups) {
      lapply(treatment_groups, self$calculate)
    }
  ),

  # ===========================================================================
  # Private - Internal computation methods
  # ===========================================================================

  private = list(

    .default_correlation = NULL,
    .apply_bias_correction = NULL,

    .validate_input = function(treatment_group) {
      if (!inherits(treatment_group, "TreatmentGroup")) {
        stop("Input must be a TreatmentGroup object")
      }

      validation <- treatment_group$validate()
      if (!validation$is_valid) {
        stop(paste("Invalid treatment group:", paste(validation$errors, collapse = "; ")))
      }
    },

    .get_correlation = function(treatment_group) {
      if (treatment_group$has_correlation) {
        treatment_group$pre_post_correlation
      } else {
        private$.default_correlation
      }
    },

    .compute_raw_effect = function(treatment_group) {
      mean_change <- treatment_group$calculate_mean_change()
      mean_change / treatment_group$standard_deviation_pre
    },

    .apply_hedges_correction = function(raw_effect, sample_size) {
      if (!private$.apply_bias_correction) {
        return(raw_effect)
      }

      correction_factor <- private$.compute_bias_correction(sample_size)
      raw_effect * correction_factor
    },

    .compute_bias_correction = function(sample_size) {
      # Exact bias correction factor using gamma function
      # This matches metafor::escalc() formula for SMCR
      # cm = exp(lgamma(df/2) - log(sqrt(df/2)) - lgamma((df-1)/2))
      degrees_of_freedom <- sample_size - 1
      exp(lgamma(degrees_of_freedom / 2) -
          log(sqrt(degrees_of_freedom / 2)) -
          lgamma((degrees_of_freedom - 1) / 2))
    },

    .compute_variance = function(effect_size, sample_size, correlation) {
      # Variance formula matching metafor::escalc(measure = "SMCR")
      # vi = 2*(1-r)/n + yi^2/(2*n)
      # where yi is the bias-corrected effect size
      term_correlation <- 2 * (1 - correlation) / sample_size
      term_effect <- effect_size^2 / (2 * sample_size)
      term_correlation + term_effect
    }
  )
)

#' ROMC Effect Size Calculator
#'
#' Calculates Response Ratio using Mean Change (log scale).
#' Used for percentage-based outcomes.
#'
#' @export
ROMCCalculator <- R6Class(

  classname = "ROMCCalculator",

  public = list(

    #' @description Create calculator with configuration
    #' @param default_correlation Default pre-post correlation
    initialize = function(default_correlation = 0.5) {
      private$.default_correlation <- default_correlation
    },

    #' @description Calculate log response ratio from TreatmentGroup
    #' @param treatment_group TreatmentGroup object
    #' @return EffectSize object
    calculate = function(treatment_group) {
      private$.validate_input(treatment_group)

      correlation <- private$.get_correlation(treatment_group)
      log_ratio <- private$.compute_log_ratio(treatment_group)
      variance <- private$.compute_variance(treatment_group, correlation)

      EffectSize$new(
        effect_estimate = log_ratio,
        sampling_variance = variance,
        sample_size = treatment_group$sample_size,
        effect_type = "romc",
        pre_post_correlation = correlation,
        bias_correction = FALSE
      )
    },

    #' @description Calculate for multiple groups
    #' @param treatment_groups List of TreatmentGroup objects
    #' @return List of EffectSize objects
    calculate_batch = function(treatment_groups) {
      lapply(treatment_groups, self$calculate)
    }
  ),

  private = list(

    .default_correlation = NULL,

    .validate_input = function(treatment_group) {
      if (!inherits(treatment_group, "TreatmentGroup")) {
        stop("Input must be a TreatmentGroup object")
      }
      if (treatment_group$mean_pre <= 0) {
        stop("Pre-intervention mean must be positive for log ratio")
      }
    },

    .get_correlation = function(treatment_group) {
      if (treatment_group$has_correlation) {
        treatment_group$pre_post_correlation
      } else {
        private$.default_correlation
      }
    },

    .compute_log_ratio = function(treatment_group) {
      log(treatment_group$mean_post / treatment_group$mean_pre)
    },

    .compute_variance = function(treatment_group, correlation) {
      n <- treatment_group$sample_size
      cv_pre <- treatment_group$standard_deviation_pre / treatment_group$mean_pre

      # Approximate variance for log response ratio
      (cv_pre^2 / n) * 2 * (1 - correlation)
    }
  )
)
