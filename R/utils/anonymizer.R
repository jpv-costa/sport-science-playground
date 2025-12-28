# R/utils/anonymizer.R
# Service: Participant Pseudoanonymization
#
# =============================================================================
# PRINCIPLES APPLIED (from CLAUDE.md)
# =============================================================================
#
# NAMING PRINCIPLES:
# - Understandability: Domain terms (anonymize, pseudonym, participant)
# - Consistency: All public methods verb-based (anonymize_X, get_mapping)
# - Distinguishability: anonymize_id vs anonymize_set_id (clear purpose)
# - Conciseness: Short meaningful names
#
# SOLID PRINCIPLES:
# - SRP: Single responsibility - name mapping (one actor: data privacy)
# - OCP: Open for extension (new anonymization patterns), closed for modification
# - DIP: Depends on abstractions (character vectors), not concrete implementations
#
# CUPID PRINCIPLES:
# - Composable: Can chain with data loader (load |> anonymize)
# - Unix: Each method does one thing (map ID, map set_id, map data)
# - Predictable: Same input -> same output (deterministic mapping)
# - Idiomatic: Follows R conventions, proper roxygen docs
# - Domain-based: Names reflect research data privacy concepts
#
# =============================================================================

box::use(
  R6[R6Class]
)

#' Participant Anonymizer
#'
#' Maps real participant names to pseudonyms (P01-P19) for data privacy.
#' Maintains consistent mapping throughout the analysis pipeline.
#'
#' @details
#' The mapping is fixed and deterministic to ensure reproducibility:
#' - Original Excel sheet names map to P01-P19
#' - Set IDs containing names are also anonymized
#' - Mapping is alphabetically ordered by original name
#'
#' @export
Anonymizer <- R6Class(

  classname = "Anonymizer",

  public = list(

    #' @description Create anonymizer with fixed mapping
    #' @return Anonymizer instance with initialized mapping
    initialize = function() {
      private$.init_mapping()
    },

    #' @description Anonymize a single participant ID
    #' @param original_id Original participant name
    #' @return Pseudonymized ID (e.g., "P01")
    #' @examples
    #' anon <- Anonymizer$new()
    #' anon$anonymize_id("Sebastião")  # Returns "P17"
    anonymize_id = function(original_id) {
      if (length(original_id) == 0) {
        return(character(0))
      }

      if (length(original_id) > 1) {
        result <- vapply(original_id, self$anonymize_id, character(1))
        return(unname(result))
      }

      pseudonym <- private$.mapping[[original_id]]

      if (is.null(pseudonym)) {
        stop(sprintf("Unknown participant ID: '%s'. Valid IDs: %s",
                     original_id, paste(names(private$.mapping), collapse = ", ")))
      }

      pseudonym
    },

    #' @description Anonymize a set_id (format: "Name_Day1_90pct_S1")
    #' @param set_id Original set ID containing participant name
    #' @return Anonymized set ID (e.g., "P17_Day1_90pct_S1")
    anonymize_set_id = function(set_id) {
      if (length(set_id) == 0) {
        return(character(0))
      }

      if (length(set_id) > 1) {
        result <- vapply(set_id, self$anonymize_set_id, character(1))
        return(unname(result))
      }

      # Extract participant name (first part before underscore)
      parts <- strsplit(set_id, "_", fixed = TRUE)[[1]]

      if (length(parts) < 2) {
        stop(sprintf("Invalid set_id format: '%s'. Expected format: 'Name_Day_Load_Set'",
                     set_id))
      }

      original_name <- parts[1]
      pseudonym <- self$anonymize_id(original_name)

      # Reconstruct set_id with pseudonym
      parts[1] <- pseudonym
      paste(parts, collapse = "_")
    },

    #' @description Anonymize entire data frame
    #' @param data Data frame containing participant data
    #' @param id_col Column name for participant ID (default: "id")
    #' @param set_id_col Column name for set ID (default: "set_id")
    #' @return Data frame with anonymized IDs
    anonymize_data = function(data, id_col = "id", set_id_col = "set_id") {
      if (!is.data.frame(data)) {
        stop("Input must be a data frame")
      }

      result <- data

      # Anonymize participant ID column
      if (id_col %in% names(result)) {
        result[[id_col]] <- self$anonymize_id(result[[id_col]])
      }

      # Anonymize set_id column if present
      if (set_id_col %in% names(result)) {
        result[[set_id_col]] <- self$anonymize_set_id(result[[set_id_col]])
      }

      result
    },

    #' @description Get full mapping table
    #' @return Data frame with original and pseudonym columns
    get_mapping = function() {
      data.frame(
        original = names(private$.mapping),
        pseudonym = unlist(private$.mapping),
        row.names = NULL,
        stringsAsFactors = FALSE
      )
    },

    #' @description Reverse lookup: get original name from pseudonym
    #' @param pseudonym Pseudonymized ID (e.g., "P01")
    #' @return Original participant name
    deanonymize_id = function(pseudonym) {
      if (length(pseudonym) == 0) {
        return(character(0))
      }

      if (length(pseudonym) > 1) {
        result <- vapply(pseudonym, self$deanonymize_id, character(1))
        return(unname(result))
      }

      reverse_mapping <- names(private$.mapping)
      names(reverse_mapping) <- unlist(private$.mapping)

      # Check if pseudonym exists in mapping
      if (!(pseudonym %in% names(reverse_mapping))) {
        stop(sprintf("Unknown pseudonym: '%s'. Valid pseudonyms: P01-P19", pseudonym))
      }

      original <- reverse_mapping[[pseudonym]]
      original
    },

    #' @description Check if a name is already anonymized
    #' @param id ID to check
    #' @return TRUE if already a pseudonym (P01-P19 format)
    is_anonymized = function(id) {
      grepl("^P[0-9]{2}$", id)
    }
  ),

  private = list(

    .mapping = NULL,

    #' Initialize the fixed mapping (alphabetically ordered originals)
    .init_mapping = function() {
      # Original participant names from Excel sheets (sorted alphabetically)
      originals <- c(
        "Alexia",
        "Chaves",
        "Filipa",
        "Leonardo",
        "Manuel M.",
        "Mário",
        "MartaM.",
        "MartaR.",
        "Miguel B.",
        "Miguel M.",
        "Miguel R.",
        "Mikhael",
        "Pedro B.",
        "Ricardo G.",
        "Ricardo L.",
        "Samuel",
        "Sebastião",
        "Simão",
        "Toni"
      )

      # Create pseudonyms P01-P19
      pseudonyms <- sprintf("P%02d", seq_along(originals))

      # Build mapping as named list
      private$.mapping <- as.list(pseudonyms)
      names(private$.mapping) <- originals
    }
  )
)
