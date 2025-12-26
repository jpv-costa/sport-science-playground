# R/domain/treatment_group.R
# Domain Entity: TreatmentGroup
#
# SOLID Principles Applied:
# - SRP: Represents a single treatment group within a study (actor: meta-analyst)
# - DIP: Depends on domain concepts, not external libraries

box::use(
  R6[R6Class]
)

#' Treatment Group Entity
#'
#' Represents a single experimental or control group within a study.
#' Contains all measurements needed for effect size calculation.
#'
#' @description
#' Domain object modeling a treatment group with pre/post measurements.
#' Mirrors the concept of an "arm" or "condition" in RCT studies.
#'
#' @export
TreatmentGroup <- R6Class(

  classname = "TreatmentGroup",

  public = list(

    # =========================================================================
    # Constructor
    # =========================================================================

    #' @description Create a new TreatmentGroup entity
    #' @param group_id Unique identifier for this group
    #' @param repetitions_in_reserve Target RIR for this group (0 = failure)
    #' @param mean_pre Pre-intervention mean
    #' @param mean_post Post-intervention mean
    #' @param standard_deviation_pre Pre-intervention standard deviation
    #' @param sample_size Number of participants
    #' @param pre_post_correlation Correlation between pre/post (optional)
    initialize = function(group_id,
                          repetitions_in_reserve,
                          mean_pre,
                          mean_post,
                          standard_deviation_pre,
                          sample_size,
                          pre_post_correlation = NULL) {

      private$.group_id <- group_id
      private$.repetitions_in_reserve <- as.numeric(repetitions_in_reserve)
      private$.mean_pre <- as.numeric(mean_pre)
      private$.mean_post <- as.numeric(mean_post)
      private$.standard_deviation_pre <- as.numeric(standard_deviation_pre)
      private$.sample_size <- as.integer(sample_size)
      private$.pre_post_correlation <- pre_post_correlation
    },

    # =========================================================================
    # Public Methods - Queries
    # =========================================================================

    #' @description Calculate raw mean change (post - pre)
    #' @return Numeric mean change value
    calculate_mean_change = function() {
      private$.mean_post - private$.mean_pre
    },

    #' @description Check if group data is valid for effect size calculation
    #' @return Named list with is_valid and errors
    validate = function() {
      errors <- character()
      errors <- private$.validate_rir(errors)
      errors <- private$.validate_standard_deviation(errors)
      errors <- private$.validate_sample_size(errors)

      list(
        is_valid = length(errors) == 0,
        errors = errors
      )
    },

    #' @description Convert to list for serialization
    #' @return Named list representation
    to_list = function() {
      list(
        group_id = private$.group_id,
        repetitions_in_reserve = private$.repetitions_in_reserve,
        mean_pre = private$.mean_pre,
        mean_post = private$.mean_post,
        standard_deviation_pre = private$.standard_deviation_pre,
        sample_size = private$.sample_size,
        pre_post_correlation = private$.pre_post_correlation,
        mean_change = self$calculate_mean_change()
      )
    }
  ),

  # ===========================================================================
  # Active Bindings
  # ===========================================================================

  active = list(

    group_id = function() private$.group_id,

    repetitions_in_reserve = function() private$.repetitions_in_reserve,

    mean_pre = function() private$.mean_pre,

    mean_post = function() private$.mean_post,

    standard_deviation_pre = function() private$.standard_deviation_pre,

    sample_size = function() private$.sample_size,

    pre_post_correlation = function() private$.pre_post_correlation,

    is_training_to_failure = function() private$.repetitions_in_reserve == 0,

    has_correlation = function() !is.null(private$.pre_post_correlation)
  ),

  # ===========================================================================
  # Private
  # ===========================================================================

  private = list(

    .group_id = NULL,
    .repetitions_in_reserve = NULL,
    .mean_pre = NULL,
    .mean_post = NULL,
    .standard_deviation_pre = NULL,
    .sample_size = NULL,
    .pre_post_correlation = NULL,

    .validate_rir = function(errors) {
      rir <- private$.repetitions_in_reserve
      if (is.null(rir) || length(rir) == 0 || is.na(rir) || rir < 0 || rir > 20) {
        errors <- c(errors, "Repetitions in reserve must be between 0 and 20")
      }
      errors
    },

    .validate_standard_deviation = function(errors) {
      sd_value <- private$.standard_deviation_pre
      if (is.null(sd_value) || length(sd_value) == 0 || is.na(sd_value) || sd_value <= 0) {
        errors <- c(errors, "Standard deviation must be positive")
      }
      errors
    },

    .validate_sample_size = function(errors) {
      n <- private$.sample_size
      if (is.null(n) || length(n) == 0 || is.na(n) || n < 2) {
        errors <- c(errors, "Sample size must be at least 2")
      }
      errors
    }
  )
)
