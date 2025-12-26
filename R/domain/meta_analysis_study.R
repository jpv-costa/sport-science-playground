# R/domain/meta_analysis_study.R
# Domain Entity: MetaAnalysisStudy
#
# SOLID Principles Applied:
# - SRP: Represents a single study/effect in meta-analysis (actor: meta-analyst)
# - OCP: Can be extended for different study types without modification

box::use(
  R6[R6Class]
)

#' Meta-Analysis Study Entity
#'
#' Represents a single effect size entry in a meta-analysis dataset.
#' Contains study-level metadata and calculated effect sizes.
#'
#' @description
#' Domain object encapsulating a meta-analysis data row with all required
#' fields for multi-level meta-regression.
#'
#' @export
MetaAnalysisStudy <- R6Class(

  classname = "MetaAnalysisStudy",

  public = list(

    # =========================================================================
    # Constructor
    # =========================================================================

    #' @description Create a new MetaAnalysisStudy entity
    #' @param study_id Unique study identifier
    #' @param group_id Group identifier within study
    #' @param observation_id Observation identifier within group
    #' @param author_year Study citation (e.g., "Smith 2020")
    #' @param outcome Outcome type ("Strength" or "Hypertrophy")
    #' @param effect_estimate Effect size (yi)
    #' @param sampling_variance Variance of effect size (vi)
    #' @param average_rir Average repetitions in reserve
    #' @param moderators Named list of moderator variables
    initialize = function(study_id,
                          group_id,
                          observation_id,
                          author_year,
                          outcome,
                          effect_estimate,
                          sampling_variance,
                          average_rir,
                          moderators = list()) {

      private$.study_id <- study_id
      private$.group_id <- group_id
      private$.observation_id <- observation_id
      private$.author_year <- author_year
      private$.outcome <- outcome
      private$.effect_estimate <- as.numeric(effect_estimate)
      private$.sampling_variance <- as.numeric(sampling_variance)
      private$.average_rir <- as.numeric(average_rir)
      private$.moderators <- moderators
    },

    # =========================================================================
    # Public Methods - Queries
    # =========================================================================

    #' @description Get a moderator value by name
    #' @param name Moderator variable name
    #' @return Moderator value or NULL if not found
    get_moderator = function(name) {
      private$.moderators[[name]]
    },

    #' @description Calculate inverse variance weight
    #' @return Numeric weight
    calculate_weight = function() {
      1 / private$.sampling_variance
    },

    #' @description Convert to list for data frame creation
    #' @return Named list representation
    to_list = function() {
      base_list <- list(
        study = private$.study_id,
        group = private$.group_id,
        obs = private$.observation_id,
        author_year = private$.author_year,
        outcome = private$.outcome,
        yi = private$.effect_estimate,
        vi = private$.sampling_variance,
        avg_rir = private$.average_rir,
        weight = self$calculate_weight()
      )
      c(base_list, private$.moderators)
    }
  ),

  # ===========================================================================
  # Active Bindings - Read-only access
  # ===========================================================================

  active = list(

    study_id = function() private$.study_id,
    group_id = function() private$.group_id,
    observation_id = function() private$.observation_id,
    author_year = function() private$.author_year,
    outcome = function() private$.outcome,
    effect_estimate = function() private$.effect_estimate,
    sampling_variance = function() private$.sampling_variance,
    average_rir = function() private$.average_rir,
    moderators = function() private$.moderators,

    # metafor-compatible accessors
    yi = function() private$.effect_estimate,
    vi = function() private$.sampling_variance,

    is_strength_outcome = function() private$.outcome == "Strength",
    is_hypertrophy_outcome = function() private$.outcome == "Hypertrophy"
  ),

  # ===========================================================================
  # Private
  # ===========================================================================

  private = list(
    .study_id = NULL,
    .group_id = NULL,
    .observation_id = NULL,
    .author_year = NULL,
    .outcome = NULL,
    .effect_estimate = NULL,
    .sampling_variance = NULL,
    .average_rir = NULL,
    .moderators = NULL
  )
)
