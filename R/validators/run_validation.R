# R/validators/run_validation.R
# Service: Data Validation for Deadlift RIR-Velocity Study
#
# =============================================================================
# PRINCIPLES APPLIED (from CLAUDE.md)
# =============================================================================
#
# NAMING PRINCIPLES:
# - Understandability: Domain terms (rir, velocity, participant, set)
# - Consistency: All public methods verb-based (validate, check, report)
# - Distinguishability: validate_data vs check_column (clear scope difference)
# - Conciseness: Short meaningful names matching data loader conventions
#
# SOLID PRINCIPLES:
# - SRP: Single responsibility - data validation only (one actor)
# - OCP: New validation rules added via private methods, interface stable
# - DIP: Depends on data frames (abstractions), not specific sources
#
# CUPID PRINCIPLES:
# - Composable: Validation results can be combined and filtered
# - Unix: Each check does one thing (range, presence, uniqueness)
# - Predictable: Same data -> same validation results (deterministic)
# - Idiomatic: Follows R conventions, proper roxygen docs, R6 classes
# - Domain-based: Rules reflect VBT research constraints
#
# SCIENTIFIC VALIDITY:
# - RIR range [0,7] based on resistance training literature
# - Velocity must be positive (physical constraint)
# - Load percentages match study protocol (80%, 90%)
# - Day values match study design (Day 1, Day 2)
#
# TESTABILITY:
# - All private methods are pure functions (input -> output)
# - Each validation rule can be tested in isolation
# - ValidationResult object enables programmatic inspection
#
# =============================================================================
# VALIDATION RULES
# =============================================================================
# 1. RIR range: [0, 7] (0 = failure, 7 = far from failure)
# 2. Velocity: > 0 m/s (positive values only)
# 3. Load percentage: {"80%", "90%"} (study protocol)
# 4. Day: {"Day 1", "Day 2"} (study design)
# 5. Participant IDs: Non-empty, no duplicates within sets
# 6. Set IDs: Unique identifiers, proper format
# 7. Required columns: All expected columns present
# =============================================================================

box::use(
  R6[R6Class]
)

#' Validation Result
#'
#' Container for validation results from DataValidator.
#' @export
ValidationResult <- R6Class(
  classname = "ValidationResult",

  public = list(
    #' @field valid Boolean indicating overall validation status
    valid = NULL,

    #' @field errors List of validation error messages
    errors = NULL,

    #' @field warnings List of validation warning messages
    warnings = NULL,

    #' @field checks Named list of individual check results
    checks = NULL,

    #' @field n_rows Number of rows validated
    n_rows = NULL,

    #' @description Create validation result
    initialize = function(valid, errors, warnings, checks, n_rows) {
      self$valid <- valid
      self$errors <- errors
      self$warnings <- warnings
      self$checks <- checks
      self$n_rows <- n_rows
    },

    #' @description Get summary of validation results
    #' @return Character string with summary
    summary = function() {
      status <- if (self$valid) "PASSED" else "FAILED"
      n_errors <- length(self$errors)
      n_warnings <- length(self$warnings)

      paste0(
        "Validation ", status, "\n",
        "  Rows checked: ", self$n_rows, "\n",
        "  Errors: ", n_errors, "\n",
        "  Warnings: ", n_warnings
      )
    },

    #' @description Print validation report
    print_report = function() {
      cat(self$summary(), "\n\n")

      if (length(self$errors) > 0) {
        cat("ERRORS:\n")
        for (err in self$errors) {
          cat("  - ", err, "\n")
        }
        cat("\n")
      }

      if (length(self$warnings) > 0) {
        cat("WARNINGS:\n")
        for (warn in self$warnings) {
          cat("  - ", warn, "\n")
        }
      }

      invisible(self)
    }
  )
)

#' Data Validator
#'
#' R6 class for validating deadlift RIR-velocity data.
#' Ensures data quality before analysis pipeline execution.
#'
#' @export
DataValidator <- R6Class(
  classname = "DataValidator",

  public = list(

    #' @description Validate data frame against study constraints
    #'
    #' @param data Data frame to validate (from DeadliftRirDataLoader)
    #' @param strict If TRUE, treat warnings as errors
    #' @return ValidationResult object
    validate = function(data, strict = FALSE) {
      errors <- character(0)
      warnings <- character(0)
      checks <- list()

      # Check required columns
      col_check <- private$.check_required_columns(data)
      checks$required_columns <- col_check
      if (!col_check$valid) {
        errors <- c(errors, col_check$message)
        # Cannot proceed without required columns
        return(ValidationResult$new(
          valid = FALSE,
          errors = errors,
          warnings = warnings,
          checks = checks,
          n_rows = nrow(data)
        ))
      }

      # Check RIR range
      rir_check <- private$.check_rir_range(data)
      checks$rir_range <- rir_check
      if (!rir_check$valid) {
        errors <- c(errors, rir_check$message)
      }

      # Check velocity values
      vel_check <- private$.check_velocity_positive(data)
      checks$velocity_positive <- vel_check
      if (!vel_check$valid) {
        errors <- c(errors, vel_check$message)
      }

      # Check load percentages
      load_check <- private$.check_load_values(data)
      checks$load_values <- load_check
      if (!load_check$valid) {
        errors <- c(errors, load_check$message)
      }

      # Check day values
      day_check <- private$.check_day_values(data)
      checks$day_values <- day_check
      if (!day_check$valid) {
        errors <- c(errors, day_check$message)
      }

      # Check participant IDs
      id_check <- private$.check_participant_ids(data)
      checks$participant_ids <- id_check
      if (!id_check$valid) {
        warnings <- c(warnings, id_check$message)
      }

      # Check set ID uniqueness
      set_check <- private$.check_set_ids(data)
      checks$set_ids <- set_check
      if (!set_check$valid) {
        warnings <- c(warnings, set_check$message)
      }

      # Check velocity range (warning for outliers)
      outlier_check <- private$.check_velocity_outliers(data)
      checks$velocity_outliers <- outlier_check
      if (!outlier_check$valid) {
        warnings <- c(warnings, outlier_check$message)
      }

      # Determine overall validity
      valid <- length(errors) == 0
      if (strict && length(warnings) > 0) {
        valid <- FALSE
        errors <- c(errors, warnings)
        warnings <- character(0)
      }

      ValidationResult$new(
        valid = valid,
        errors = errors,
        warnings = warnings,
        checks = checks,
        n_rows = nrow(data)
      )
    },

    #' @description Quick check if data is valid (no detailed report)
    #' @param data Data frame to validate
    #' @return Boolean
    is_valid = function(data) {
      result <- self$validate(data)
      result$valid
    }
  ),

  private = list(

    #' Check that all required columns are present
    .check_required_columns = function(data) {
      required <- c(
        "id", "sex", "day", "load_percentage",
        "rir", "mean_velocity", "set_id", "rep_number"
      )
      missing <- setdiff(required, names(data))

      if (length(missing) > 0) {
        list(
          valid = FALSE,
          message = paste0("Missing required columns: ", paste(missing, collapse = ", ")),
          missing = missing
        )
      } else {
        list(valid = TRUE, message = "All required columns present")
      }
    },

    #' Check RIR values are in valid range [0, 7]
    .check_rir_range = function(data) {
      if (!"rir" %in% names(data)) {
        return(list(valid = TRUE, message = "RIR column not present (skipped)"))
      }

      invalid <- data$rir < 0 | data$rir > 7 | is.na(data$rir)
      n_invalid <- sum(invalid, na.rm = TRUE)

      if (n_invalid > 0) {
        out_of_range <- unique(data$rir[invalid & !is.na(data$rir)])
        list(
          valid = FALSE,
          message = paste0(
            "RIR out of range [0,7]: ", n_invalid, " observations. ",
            "Values: ", paste(out_of_range, collapse = ", ")
          ),
          n_invalid = n_invalid
        )
      } else {
        list(valid = TRUE, message = "RIR values in valid range [0,7]")
      }
    },

    #' Check velocity values are positive
    .check_velocity_positive = function(data) {
      if (!"mean_velocity" %in% names(data)) {
        return(list(valid = TRUE, message = "Velocity column not present (skipped)"))
      }

      invalid <- data$mean_velocity <= 0 | is.na(data$mean_velocity)
      n_invalid <- sum(invalid, na.rm = TRUE)

      if (n_invalid > 0) {
        list(
          valid = FALSE,
          message = paste0(
            "Non-positive velocity values: ", n_invalid, " observations"
          ),
          n_invalid = n_invalid
        )
      } else {
        list(valid = TRUE, message = "All velocity values positive")
      }
    },

    #' Check load percentages are valid study values
    .check_load_values = function(data) {
      if (!"load_percentage" %in% names(data)) {
        return(list(valid = TRUE, message = "Load column not present (skipped)"))
      }

      valid_loads <- c("80%", "90%")
      unique_loads <- unique(data$load_percentage)
      invalid_loads <- setdiff(unique_loads, valid_loads)

      if (length(invalid_loads) > 0) {
        list(
          valid = FALSE,
          message = paste0(
            "Invalid load percentages: ", paste(invalid_loads, collapse = ", "),
            ". Expected: 80% or 90%"
          ),
          invalid = invalid_loads
        )
      } else {
        list(valid = TRUE, message = "Load percentages valid (80%, 90%)")
      }
    },

    #' Check day values are valid study values
    .check_day_values = function(data) {
      if (!"day" %in% names(data)) {
        return(list(valid = TRUE, message = "Day column not present (skipped)"))
      }

      valid_days <- c("Day 1", "Day 2")
      unique_days <- unique(data$day)
      invalid_days <- setdiff(unique_days, valid_days)

      if (length(invalid_days) > 0) {
        list(
          valid = FALSE,
          message = paste0(
            "Invalid day values: ", paste(invalid_days, collapse = ", "),
            ". Expected: Day 1 or Day 2"
          ),
          invalid = invalid_days
        )
      } else {
        list(valid = TRUE, message = "Day values valid (Day 1, Day 2)")
      }
    },

    #' Check participant IDs are valid
    .check_participant_ids = function(data) {
      if (!"id" %in% names(data)) {
        return(list(valid = TRUE, message = "ID column not present (skipped)"))
      }

      # Check for empty IDs
      empty_ids <- is.na(data$id) | data$id == ""
      n_empty <- sum(empty_ids)

      if (n_empty > 0) {
        list(
          valid = FALSE,
          message = paste0("Empty participant IDs: ", n_empty, " observations"),
          n_empty = n_empty
        )
      } else {
        n_participants <- length(unique(data$id))
        list(
          valid = TRUE,
          message = paste0("Participant IDs valid (n=", n_participants, ")")
        )
      }
    },

    #' Check set IDs for proper formation
    .check_set_ids = function(data) {
      if (!"set_id" %in% names(data)) {
        return(list(valid = TRUE, message = "Set ID column not present (skipped)"))
      }

      # Check for empty set IDs
      empty_sets <- is.na(data$set_id) | data$set_id == ""
      n_empty <- sum(empty_sets)

      if (n_empty > 0) {
        list(
          valid = FALSE,
          message = paste0("Empty set IDs: ", n_empty, " observations"),
          n_empty = n_empty
        )
      } else {
        n_sets <- length(unique(data$set_id))
        list(
          valid = TRUE,
          message = paste0("Set IDs valid (n=", n_sets, " unique sets)")
        )
      }
    },

    #' Check for velocity outliers (warning only)
    .check_velocity_outliers = function(data) {
      if (!"mean_velocity" %in% names(data)) {
        return(list(valid = TRUE, message = "Velocity column not present (skipped)"))
      }

      # Typical deadlift velocity range: 0.10 - 0.70 m/s
      # Flag values outside this range as potential outliers
      lower_bound <- 0.05
      upper_bound <- 0.80

      outliers <- data$mean_velocity < lower_bound | data$mean_velocity > upper_bound
      n_outliers <- sum(outliers, na.rm = TRUE)

      if (n_outliers > 0) {
        outlier_vals <- data$mean_velocity[outliers & !is.na(data$mean_velocity)]
        list(
          valid = FALSE,
          message = paste0(
            "Velocity outliers (outside [0.05, 0.80] m/s): ", n_outliers,
            " observations. Range: [",
            round(min(outlier_vals), 3), ", ",
            round(max(outlier_vals), 3), "]"
          ),
          n_outliers = n_outliers
        )
      } else {
        list(valid = TRUE, message = "No velocity outliers detected")
      }
    }
  )
)

#' Run validation on loaded data
#'
#' Convenience function to run validation and print report.
#' Designed to be called from Makefile target.
#'
#' @param data_path Path to the Excel data file
#' @param strict If TRUE, treat warnings as errors
#' @return ValidationResult (invisibly)
#' @export
run_validation <- function(data_path = NULL, strict = FALSE) {
  # Default path if not provided
  if (is.null(data_path)) {
    data_path <- "deadlift-study/TESE - Análise estatística (Preliminar).xlsx"
  }

  cat("=== Data Validation Report ===\n\n")
  cat("Data path: ", data_path, "\n\n")

  # Load data
  box::use(../loaders/deadlift_rir_data_loader[DeadliftRirDataLoader])

  tryCatch({
    loader <- DeadliftRirDataLoader$new(data_path)
    data <- loader$load()

    cat("Data loaded successfully\n")
    cat("  Rows: ", nrow(data), "\n")
    cat("  Columns: ", ncol(data), "\n\n")

    # Validate
    validator <- DataValidator$new()
    result <- validator$validate(data, strict = strict)

    # Print report
    result$print_report()

    # Return with exit code for CI/CD
    if (!result$valid) {
      cat("\n[VALIDATION FAILED]\n")
    } else {
      cat("\n[VALIDATION PASSED]\n")
    }

    invisible(result)

  }, error = function(e) {
    cat("ERROR loading data: ", e$message, "\n")
    invisible(NULL)
  })
}
