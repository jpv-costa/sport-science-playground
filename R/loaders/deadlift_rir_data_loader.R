# R/loaders/deadlift_rir_data_loader.R
# Loader: Deadlift RIR-Velocity Data (Thesis Research)
#
# SOLID Principles Applied:
# - SRP: Single responsibility - load and preprocess deadlift RIR-velocity data
#
# Output columns:
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

box::use(
  R6[R6Class]
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
    #' @return Data frame with cleaned and combined data in long format
    load = function() {
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
    .parse_participant_sheet = function(sheet_name) {
      raw_data <- as.data.frame(
        readxl::read_excel(
          private$.data_path,
          sheet = sheet_name,
          col_names = FALSE
        )
      )

      # Extract condition info from row 1 headers
      headers <- as.character(raw_data[1, ])
      conditions <- private$.parse_headers(headers)

      if (length(conditions) == 0) {
        return(NULL)
      }

      # Data starts at row 3 (row 1 = headers, row 2 = column labels)
      data_rows <- raw_data[3:nrow(raw_data), ]

      # Determine participant sex from sheet name
      sex <- ifelse(grepl("^Mulheres_", sheet_name), "female", "male")

      # Clean participant name
      id <- gsub("^Mulheres_", "", sheet_name)

      # Parse each condition (set of 3 columns)
      all_observations <- list()

      for (i in seq_along(conditions)) {
        condition <- conditions[[i]]
        col_start <- condition$col_start

        # Extract RIR and MCV columns
        rir_col <- col_start + 1
        mcv_col <- col_start + 2

        rir_values <- as.numeric(data_rows[, rir_col])
        mcv_values <- as.numeric(data_rows[, mcv_col])

        # Filter to rows with actual data
        valid_rows <- !is.na(rir_values) & !is.na(mcv_values)
        n_valid <- sum(valid_rows)

        if (n_valid == 0) {
          next
        }

        # Create unique set_id: "id_DayX_YYpct_SZ"
        day_code <- gsub("Day ", "Day", condition$day)
        load_code <- gsub("%", "pct", condition$load_pct)
        series_code <- paste0("S", condition$series_num)
        set_id <- paste(id, day_code, load_code, series_code, sep = "_")

        # Rep numbers are 1, 2, 3, ... within the set
        # The first row in the data corresponds to the first rep
        rep_numbers <- seq_len(n_valid)

        condition_data <- data.frame(
          id = id,
          sex = sex,
          day = condition$day,
          load_percentage = condition$load_pct,
          weight_kg = condition$weight_kg,
          rir = rir_values[valid_rows],
          mean_velocity = mcv_values[valid_rows],
          set_id = set_id,
          rep_number = rep_numbers,
          stringsAsFactors = FALSE
        )

        all_observations <- c(all_observations, list(condition_data))
      }

      if (length(all_observations) == 0) {
        return(NULL)
      }

      do.call(rbind, all_observations)
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
