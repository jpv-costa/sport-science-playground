# R/loaders/velocity_data_loader.R
# Loader: PeerJ Velocity-RIR Data
#
# SOLID Principles Applied:
# - SRP: Single responsibility - load and preprocess velocity data (actor: data engineer)

box::use(
  R6[R6Class]
)

#' Velocity Data Loader
#'
#' Loads and preprocesses data from Paulsen et al. (2025) PeerJ study
#' on velocity-RIR relationships.
#'
#' @export
VelocityDataLoader <- R6Class(

  classname = "VelocityDataLoader",

  public = list(

    #' @description Create loader with data path
    #' @param data_path Path to the PeerJ Excel data file
    initialize = function(data_path) {
      private$.data_path <- data_path
    },

    #' @description Load and preprocess velocity data
    #' @return Data frame with cleaned column names
    load = function() {
      if (!requireNamespace("readxl", quietly = TRUE)) {
        stop("Package 'readxl' is required to load Excel data")
      }

      raw_data <- as.data.frame(readxl::read_xlsx(private$.data_path, sheet = 1))
      cleaned <- private$.clean_column_names(raw_data)
      cleaned
    },

    #' @description Get summary statistics
    #' @param data Loaded data frame
    #' @return Named list with summary stats
    summarize = function(data) {
      list(
        n_observations = nrow(data),
        n_participants = length(unique(data$id)),
        exercises = unique(data$exercise),
        velocity_range = range(data$mean_velocity, na.rm = TRUE),
        rir_range = range(data$perceived_rir, na.rm = TRUE)
      )
    }
  ),

  private = list(

    .data_path = NULL,

    .clean_column_names = function(data) {
      # Standardize column names
      names(data) <- tolower(names(data))
      names(data) <- gsub(" ", "_", names(data))
      names(data) <- gsub("\\.", "_", names(data))

      # Remove unit annotations from column names
      names(data) <- gsub("\\[.*\\]", "", names(data))
      names(data) <- gsub("\\$.*\\$", "", names(data))
      names(data) <- trimws(names(data))

      # Rename key columns for clarity
      if ("perceived_reps_in_reserve_" %in% names(data)) {
        names(data)[names(data) == "perceived_reps_in_reserve_"] <- "perceived_rir"
      }
      if ("mean_velocity_" %in% names(data)) {
        names(data)[names(data) == "mean_velocity_"] <- "mean_velocity"
      }
      if ("velocity_loss_" %in% names(data)) {
        names(data)[names(data) == "velocity_loss_"] <- "velocity_loss"
      }

      data
    }
  )
)
