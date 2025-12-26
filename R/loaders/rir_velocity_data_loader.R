# R/loaders/rir_velocity_data_loader.R
# Loader: RIR-Velocity Relationship Data (Jukic et al. 2024)
#
# SOLID Principles Applied:
# - SRP: Single responsibility - load and preprocess RIR-velocity data

box::use(
  R6[R6Class]
)

#' RIR-Velocity Data Loader
#'
#' Loads and preprocesses data from Jukic et al. (2024) study on
#' modeling the RIR-velocity relationship.
#'
#' @export
RirVelocityDataLoader <- R6Class(

  classname = "RirVelocityDataLoader",

  public = list(

    #' @description Create loader with data path
    #' @param data_path Path to the Excel data file
    initialize = function(data_path) {
      private$.data_path <- data_path
    },

    #' @description Load and preprocess RIR-velocity data
    #' @return Data frame with cleaned and derived variables
    load = function() {
      if (!requireNamespace("readxl", quietly = TRUE)) {
        stop("Package 'readxl' is required to load Excel data")
      }

      raw_data <- as.data.frame(readxl::read_xlsx(private$.data_path, sheet = 1))
      cleaned <- private$.preprocess(raw_data)
      cleaned
    },

    #' @description Filter data by day
    #' @param data Loaded data frame
    #' @param day Day to filter ("Day 1" or "Day 2")
    #' @return Filtered data frame
    filter_by_day = function(data, day) {
      data[data$day == day, ]
    },

    #' @description Filter data by load (Set.Type)
    #' @param data Loaded data frame
    #' @param load Load percentage ("RTF70", "RTF80", "RTF90")
    #' @return Filtered data frame
    filter_by_load = function(data, load) {
      data[data$set_type == load, ]
    },

    #' @description Get summary statistics
    #' @param data Loaded data frame
    #' @return Named list with summary stats
    summarize = function(data) {
      list(
        n_observations = nrow(data),
        n_participants = length(unique(data$id)),
        n_male = sum(data$sex == "male") / nrow(data) * length(unique(data$id)),
        n_female = sum(data$sex == "female") / nrow(data) * length(unique(data$id)),
        load_types = unique(data$set_type),
        days = unique(data$day),
        velocity_range = range(data$mean_velocity, na.rm = TRUE),
        rir_range = range(data$rir, na.rm = TRUE),
        mean_relative_strength = mean(data$relative_strength, na.rm = TRUE)
      )
    }
  ),

  private = list(

    .data_path = NULL,

    .preprocess = function(data) {
      # Standardize column names
      names(data) <- tolower(names(data))
      names(data) <- gsub(" ", "_", names(data))
      names(data) <- gsub("\\.", "_", names(data))

      # Rename key columns for clarity
      names(data)[names(data) == "right_gym_mv"] <- "mean_velocity"
      names(data)[names(data) == "right_gym_pv"] <- "peak_velocity"
      names(data)[names(data) == "1rm/bm"] <- "relative_strength"
      names(data)[names(data) == "set_type"] <- "set_type"

      # Derive velocity loss and max velocity per set
      data <- private$.calculate_velocity_metrics(data)

      # Convert categorical variables
      data$id <- as.factor(data$id)
      data$sex <- as.factor(data$sex)
      data$set_type <- as.factor(data$set_type)
      data$day <- as.factor(data$day)

      data
    },

    .calculate_velocity_metrics = function(data) {
      # Group by ID, Day, Set.Type to calculate velocity metrics
      groups <- split(data, list(data$id, data$day, data$set_type))

      processed <- lapply(groups, function(group) {
        if (nrow(group) == 0) return(NULL)

        group$max_velocity <- max(group$mean_velocity, na.rm = TRUE)
        group$velocity_last <- group$mean_velocity[nrow(group)]
        group$total_velocity_loss <- (group$velocity_last - group$max_velocity) /
          group$max_velocity * -100
        group$velocity_percentage <- (group$mean_velocity - group$max_velocity) /
          group$max_velocity * -100
        group$max_reps <- max(group$rep, na.rm = TRUE)

        group
      })

      do.call(rbind, processed[!sapply(processed, is.null)])
    }
  )
)
