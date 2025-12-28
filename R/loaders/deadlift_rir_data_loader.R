# R/loaders/deadlift_rir_data_loader.R
# Service: Deadlift RIR-Velocity Data Loading and Preprocessing
#
# =============================================================================
# PRINCIPLES APPLIED (from CLAUDE.md)
# =============================================================================
#
# NAMING PRINCIPLES:
# - Understandability: Problem domain terms (rir, velocity, set_id, day)
# - Consistency: All public methods verb-based (load, filter_by_X, summarize)
# - Distinguishability: filter_by_day vs filter_by_load (clear differences)
# - Conciseness: Short meaningful names (rir for Repetitions in Reserve)
#
# SOLID PRINCIPLES:
# - SRP: Single responsibility - Excel parsing and preprocessing (one actor)
# - OCP: Open for extension with new filter methods, closed for modification
# - DIP: Depends on abstractions (data frames), not concrete implementations
#
# CUPID PRINCIPLES:
# - Composable: Filter methods can be chained (filter_by_day |> filter_by_load)
# - Unix: Each method does one thing (load, filter, summarize)
# - Predictable: Same data file -> same output (deterministic)
# - Idiomatic: Follows R conventions, proper roxygen docs
# - Domain-based: Names reflect VBT research concepts (RIR, velocity, sets)
#
# SCIENTIFIC VALIDITY:
# - RIR=0 represents muscular failure (validated in biomechanics literature)
# - Set ID formation ensures traceability from raw Excel to processed data
# - Rep numbering maintains temporal sequence for within-set analysis
#
# TESTABILITY:
# - All private methods are pure functions (input -> output)
# - No side effects in data transformation
# - Each method can be tested in isolation
#
# =============================================================================
# OUTPUT COLUMNS
# =============================================================================
# - id: Participant identifier (string)
# - sex: "male" or "female"
# - day: "Day 1" or "Day 2"
# - load_percentage: "80%" or "90%"
# - weight_kg: Absolute load in kg
# - rir: Repetitions in Reserve (0 at failure, higher = further from failure)
# - mean_velocity: Mean concentric velocity (m/s)
# - set_id: Unique set identifier (e.g., "Alexia_Day1_90pct_S1")
# - rep_number: Position in set (1, 2, 3, ... n)
# - reps_to_failure: Total reps completed in that set
# =============================================================================

box::use(
  R6[R6Class],
  ../utils/anonymizer[Anonymizer]
)

#' Deadlift RIR-Velocity Data Loader
#'
#' Loads and preprocesses thesis data on deadlift RIR-velocity relationships.
#' Handles multi-sheet Excel format with one sheet per participant.
#'
#' @export
DeadliftRirDataLoader <- R6Class(

  classname = "DeadliftRirDataLoader",

  public = list(

    #' @description Create loader with data path
    #' @param data_path Path to the Excel data file with participant sheets
    initialize = function(data_path) {
      private$.data_path <- data_path
    },

    #' @description Load and preprocess all participant data
    #' @param anonymize Logical. If TRUE, replace participant names with
    #'   pseudonyms (P01-P19) for data privacy. Default: FALSE.
    #' @return Data frame with cleaned and combined data in long format
    load = function(anonymize = FALSE) {
      if (!requireNamespace("readxl", quietly = TRUE)) {
        stop("Package 'readxl' is required to load Excel data")
      }

      sheets <- readxl::excel_sheets(private$.data_path)
      all_data <- lapply(sheets, private$.parse_participant_sheet)
      combined <- do.call(rbind, all_data[!sapply(all_data, is.null)])

      # Clean up row names
      rownames(combined) <- NULL

      # Add reps_to_failure (total reps per set)
      combined <- private$.add_reps_to_failure(combined)

      # Anonymize participant IDs if requested
      if (anonymize) {
        anon <- Anonymizer$new()
        combined <- anon$anonymize_data(combined, id_col = "id", set_id_col = "set_id")
      }

      combined
    },

    #' @description Filter data by day
    #' @param data Loaded data frame
    #' @param day Day to filter ("Day 1" or "Day 2")
    #' @return Filtered data frame
    filter_by_day = function(data, day) {
      data[data$day == day, ]
    },

    #' @description Filter data by load percentage
    #' @param data Loaded data frame
    #' @param load Load percentage ("80%" or "90%")
    #' @return Filtered data frame
    filter_by_load = function(data, load) {
      data[data$load_percentage == load, ]
    },

    #' @description Get summary statistics
    #' @param data Loaded data frame
    #' @return Named list with summary stats
    summarize = function(data) {
      list(
        n_observations = nrow(data),
        n_participants = length(unique(data$id)),
        n_sets = length(unique(data$set_id)),
        n_male = sum(data$sex == "male") / nrow(data) * length(unique(data$id)),
        n_female = sum(data$sex == "female") / nrow(data) * length(unique(data$id)),
        load_types = unique(data$load_percentage),
        days = unique(data$day),
        velocity_range = range(data$mean_velocity, na.rm = TRUE),
        rir_range = range(data$rir, na.rm = TRUE),
        weight_range = range(data$weight_kg, na.rm = TRUE),
        reps_per_set_range = range(data$reps_to_failure, na.rm = TRUE),
        mean_reps_per_set = mean(
          tapply(data$reps_to_failure, data$set_id, function(x) x[1]),
          na.rm = TRUE
        )
      )
    }
  ),

  private = list(

    .data_path = NULL,

    #' Parse a single participant sheet into long format
    #' Orchestrates sheet parsing using smaller focused methods
    .parse_participant_sheet = function(sheet_name) {
      raw_data <- private$.read_sheet_data(sheet_name)
      conditions <- private$.parse_headers(as.character(raw_data[1, ]))

      if (length(conditions) == 0) return(NULL)

      participant <- private$.extract_participant_info(sheet_name)
      data_rows <- raw_data[3:nrow(raw_data), ]

      observations <- private$.parse_all_conditions(
        conditions, data_rows, participant
      )

      if (length(observations) == 0) return(NULL)
      do.call(rbind, observations)
    },

    #' Read raw Excel sheet data
    .read_sheet_data = function(sheet_name) {
      as.data.frame(
        readxl::read_excel(
          private$.data_path,
          sheet = sheet_name,
          col_names = FALSE
        )
      )
    },

    #' Extract participant ID and sex from sheet name
    .extract_participant_info = function(sheet_name) {
      list(
        sex = ifelse(grepl("^Mulheres_", sheet_name), "female", "male"),
        id = gsub("^Mulheres_", "", sheet_name)
      )
    },

    #' Parse all conditions for a participant
    .parse_all_conditions = function(conditions, data_rows, participant) {
      observations <- list()

      for (condition in conditions) {
        condition_obs <- private$.parse_single_condition(
          condition, data_rows, participant
        )
        if (!is.null(condition_obs)) {
          observations <- c(observations, list(condition_obs))
        }
      }

      observations
    },

    #' Parse a single condition (set) into observations
    .parse_single_condition = function(condition, data_rows, participant) {
      values <- private$.extract_rir_velocity_values(condition, data_rows)

      if (values$n_valid == 0) return(NULL)

      set_id <- private$.build_set_id(participant$id, condition)

      data.frame(
        id = participant$id,
        sex = participant$sex,
        day = condition$day,
        load_percentage = condition$load_pct,
        weight_kg = condition$weight_kg,
        rir = values$rir,
        mean_velocity = values$velocity,
        set_id = set_id,
        rep_number = seq_len(values$n_valid),
        stringsAsFactors = FALSE
      )
    },

    #' Extract RIR and velocity values from data columns
    .extract_rir_velocity_values = function(condition, data_rows) {
      rir_col <- condition$col_start + 1
      mcv_col <- condition$col_start + 2

      rir_values <- as.numeric(data_rows[, rir_col])
      mcv_values <- as.numeric(data_rows[, mcv_col])

      valid_rows <- !is.na(rir_values) & !is.na(mcv_values)

      list(
        rir = rir_values[valid_rows],
        velocity = mcv_values[valid_rows],
        n_valid = sum(valid_rows)
      )
    },

    #' Build unique set identifier
    .build_set_id = function(participant_id, condition) {
      day_code <- gsub("Day ", "Day", condition$day)
      load_code <- gsub("%", "pct", condition$load_pct)
      series_code <- paste0("S", condition$series_num)
      paste(participant_id, day_code, load_code, series_code, sep = "_")
    },

    #' Parse header row to extract condition information
    .parse_headers = function(headers) {
      conditions <- list()
      col_idx <- 1

      while (col_idx <= length(headers)) {
        header <- headers[col_idx]

        if (!is.na(header) && grepl("Série", header)) {
          # Parse header: "Série X - YY% (ZZZ kg) - Dia N"
          parsed <- private$.extract_condition_from_header(header)

          if (!is.null(parsed)) {
            parsed$col_start <- col_idx
            conditions <- c(conditions, list(parsed))
          }
        }

        col_idx <- col_idx + 1
      }

      conditions
    },

    #' Extract series number, load, weight, and day from header string
    .extract_condition_from_header = function(header) {
      # Pattern: "Série X - YY% (ZZZ kg) - Dia N"
      series_match <- regmatches(header, regexpr("Série\\s*\\d+", header))
      load_match <- regmatches(header, regexpr("\\d+%", header))
      weight_match <- regmatches(header, regexpr("\\d+\\.?\\d*\\s*kg", header))
      day_match <- regmatches(header, regexpr("Dia\\s*\\d+", header))

      if (length(load_match) == 0 || length(day_match) == 0) {
        return(NULL)
      }

      # Extract series number (default to 1 if not found)
      series_num <- 1
      if (length(series_match) > 0) {
        series_num <- as.integer(gsub("Série\\s*", "", series_match))
      }

      # Extract weight value
      weight_kg <- NA_real_
      if (length(weight_match) > 0) {
        weight_kg <- as.numeric(gsub("\\s*kg", "", weight_match))
      }

      # Format day
      day_num <- gsub("Dia\\s*", "", day_match)
      day <- paste("Day", day_num)

      list(
        series_num = series_num,
        load_pct = load_match,
        weight_kg = weight_kg,
        day = day
      )
    },

    #' Add reps_to_failure column (total reps per set)
    .add_reps_to_failure = function(data) {
      # Calculate max rep_number per set_id
      reps_per_set <- tapply(data$rep_number, data$set_id, max)

      # Map back to each observation
      data$reps_to_failure <- reps_per_set[data$set_id]

      data
    }
  )
)
