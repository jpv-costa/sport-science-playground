# R/domain/study.R
# Domain Entity: Study
#
# SOLID Principles Applied:
# - SRP: Represents a single study entity (actor: meta-analyst)
# - OCP: Extended via composition with TreatmentGroup objects
# - LSP: All study subclasses must fulfill Study contract
# - DIP: Depends on abstract validation, not concrete validators

box::use(
  R6[R6Class]
)

#' Study Entity
#'
#' Represents a single research study in the meta-analysis.
#' Domain object modeling a published study examining RIR and training outcomes.
#'
#' @description
#' This class encapsulates all study-level metadata and treatment groups.
#' It mirrors the business domain concept of a "study" in systematic reviews.
#'
#' @export
Study <- R6Class(

  classname = "Study",

  public = list(

    # =========================================================================
    # Constructor
    # =========================================================================

    #' @description Create a new Study entity
    #' @param study_id Unique identifier for the study
    #' @param first_author First author's surname
    #' @param publication_year Year of publication
    #' @param outcome_category "strength" or "hypertrophy"
    #' @param treatment_groups List of TreatmentGroup objects (optional)
    initialize = function(study_id,
                          first_author,
                          publication_year,
                          outcome_category,
                          treatment_groups = list()) {

      private$.study_id <- study_id
      private$.first_author <- first_author
      private$.publication_year <- as.integer(publication_year)
      private$.outcome_category <- outcome_category
      private$.treatment_groups <- treatment_groups
      private$.created_at <- Sys.time()
    },

    # =========================================================================
    # Public Methods - Commands (modify state)
    # =========================================================================

    #' @description Add a treatment group to this study
    #' @param treatment_group TreatmentGroup object to add
    #' @return self (enables method chaining)
    add_treatment_group = function(treatment_group) {
      private$.treatment_groups <- c(private$.treatment_groups, list(treatment_group))
      invisible(self)
    },

    # =========================================================================
    # Public Methods - Queries (return data, no side effects)
    # =========================================================================

    #' @description Check if study data is valid
    #' @return Named list with is_valid (logical) and errors (character vector)
    validate = function() {
      errors <- character()
      errors <- private$.validate_publication_year(errors)
      errors <- private$.validate_outcome_category(errors)
      errors <- private$.validate_treatment_groups(errors)

      list(
        is_valid = length(errors) == 0,
        errors = errors
      )
    },

    #' @description Convert study to a list for serialization
    #' @return Named list representation
    to_list = function() {
      list(
        study_id = private$.study_id,
        first_author = private$.first_author,
        publication_year = private$.publication_year,
        outcome_category = private$.outcome_category,
        treatment_group_count = self$treatment_group_count,
        created_at = private$.created_at
      )
    },

    #' @description Format study as citation string
    #' @return Character string like "Smith et al. (2023)"
    format_citation = function() {
      sprintf("%s et al. (%d)", private$.first_author, private$.publication_year)
    }
  ),

  # ===========================================================================
  # Active Bindings - Properties (read-only access to private fields)
  # ===========================================================================

  active = list(

    study_id = function() private$.study_id,

    first_author = function() private$.first_author,

    publication_year = function() private$.publication_year,

    outcome_category = function() private$.outcome_category,

    treatment_groups = function() private$.treatment_groups,

    treatment_group_count = function() length(private$.treatment_groups),

    is_strength_study = function() private$.outcome_category == "strength",

    is_hypertrophy_study = function() private$.outcome_category == "hypertrophy"
  ),

  # ===========================================================================
  # Private - Internal state and helper methods
  # ===========================================================================

  private = list(

    # State
    .study_id = NULL,
    .first_author = NULL,
    .publication_year = NULL,
    .outcome_category = NULL,
    .treatment_groups = NULL,
    .created_at = NULL,

    # Validation helpers - each returns accumulated errors
    .validate_publication_year = function(errors) {
      year <- private$.publication_year
      if (is.na(year) || year < 1980 || year > 2030) {
        errors <- c(errors, "Publication year must be between 1980 and 2030")
      }
      errors
    },

    .validate_outcome_category = function(errors) {
      valid_categories <- c("strength", "hypertrophy")
      if (!private$.outcome_category %in% valid_categories) {
        errors <- c(errors, sprintf(
          "Outcome category must be one of: %s",
          paste(valid_categories, collapse = ", ")
        ))
      }
      errors
    },

    .validate_treatment_groups = function(errors) {
      for (group in private$.treatment_groups) {
        if (!inherits(group, "TreatmentGroup")) {
          errors <- c(errors, "All treatment groups must be TreatmentGroup objects")
          break
        }
      }
      errors
    }
  )
)
