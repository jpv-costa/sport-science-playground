# R/calculators/rir_velocity_modeler.R
# Service: RIR-Velocity Relationship Modeling
#
# SOLID Principles Applied:
# - SRP: Single responsibility - fit and evaluate RIR-velocity models
# - OCP: Open for extension with different model types

box::use(
  R6[R6Class]
)

#' RIR-Velocity Model Result
#'
#' Value object containing model fit results.
#'
#' @export
RirVelocityModelResult <- R6Class(

  classname = "RirVelocityModelResult",
  cloneable = FALSE,

  public = list(

    #' @description Create model result
    #' @param model_type Type of model ("linear" or "polynomial")
    #' @param r_squared R-squared value
    #' @param rse Residual standard error
    #' @param coefficients Model coefficients
    #' @param n Sample size
    #' @param group_id Optional group identifier
    initialize = function(model_type, r_squared, rse, coefficients, n, group_id = NULL) {
      private$.model_type <- model_type
      private$.r_squared <- r_squared
      private$.rse <- rse
      private$.coefficients <- coefficients
      private$.n <- n
      private$.group_id <- group_id
    },

    #' @description Convert to list
    to_list = function() {
      list(
        model_type = private$.model_type,
        r_squared = private$.r_squared,
        rse = private$.rse,
        coefficients = private$.coefficients,
        n = private$.n,
        group_id = private$.group_id
      )
    }
  ),

  active = list(
    model_type = function() private$.model_type,
    r_squared = function() private$.r_squared,
    rse = function() private$.rse,
    coefficients = function() private$.coefficients,
    n = function() private$.n,
    group_id = function() private$.group_id
  ),

  private = list(
    .model_type = NULL,
    .r_squared = NULL,
    .rse = NULL,
    .coefficients = NULL,
    .n = NULL,
    .group_id = NULL
  )
)


#' RIR-Velocity Model Summary
#'
#' Value object containing summary statistics across multiple models.
#'
#' @export
RirVelocityModelSummary <- R6Class(

  classname = "RirVelocityModelSummary",
  cloneable = FALSE,

  public = list(

    #' @description Create model summary
    #' @param results List of RirVelocityModelResult objects
    #' @param model_type Model type being summarized
    initialize = function(results, model_type) {
      private$.results <- results
      private$.model_type <- model_type

      r_squared_values <- sapply(results, function(r) r$r_squared)
      rse_values <- sapply(results, function(r) r$rse)

      # Remove NA values
      r_squared_values <- r_squared_values[!is.na(r_squared_values)]
      rse_values <- rse_values[!is.na(rse_values)]

      private$.n_models <- length(results)
      private$.median_r_squared <- stats::median(r_squared_values, na.rm = TRUE)
      private$.mean_r_squared <- mean(r_squared_values, na.rm = TRUE)
      private$.min_r_squared <- min(r_squared_values, na.rm = TRUE)
      private$.max_r_squared <- max(r_squared_values, na.rm = TRUE)
      private$.median_rse <- stats::median(rse_values, na.rm = TRUE)
      private$.mean_rse <- mean(rse_values, na.rm = TRUE)
      private$.min_rse <- min(rse_values, na.rm = TRUE)
      private$.max_rse <- max(rse_values, na.rm = TRUE)
    },

    #' @description Convert to list
    to_list = function() {
      list(
        model_type = private$.model_type,
        n_models = private$.n_models,
        median_r_squared = private$.median_r_squared,
        mean_r_squared = private$.mean_r_squared,
        min_r_squared = private$.min_r_squared,
        max_r_squared = private$.max_r_squared,
        median_rse = private$.median_rse,
        mean_rse = private$.mean_rse,
        min_rse = private$.min_rse,
        max_rse = private$.max_rse
      )
    }
  ),

  active = list(
    model_type = function() private$.model_type,
    n_models = function() private$.n_models,
    median_r_squared = function() private$.median_r_squared,
    mean_r_squared = function() private$.mean_r_squared,
    min_r_squared = function() private$.min_r_squared,
    max_r_squared = function() private$.max_r_squared,
    median_rse = function() private$.median_rse,
    mean_rse = function() private$.mean_rse,
    min_rse = function() private$.min_rse,
    max_rse = function() private$.max_rse,
    results = function() private$.results
  ),

  private = list(
    .results = NULL,
    .model_type = NULL,
    .n_models = NULL,
    .median_r_squared = NULL,
    .mean_r_squared = NULL,
    .min_r_squared = NULL,
    .max_r_squared = NULL,
    .median_rse = NULL,
    .mean_rse = NULL,
    .min_rse = NULL,
    .max_rse = NULL
  )
)


#' RIR-Velocity Modeler
#'
#' Fits and evaluates RIR-velocity relationship models.
#' Supports linear and polynomial (quadratic) regression.
#'
#' @export
RirVelocityModeler <- R6Class(

  classname = "RirVelocityModeler",

  public = list(

    #' @description Create modeler
    initialize = function() {
      # No configuration needed
    },

    #' @description Fit linear model
    #' @param data Data frame with rir and mean_velocity columns
    #' @param group_id Optional group identifier
    #' @return RirVelocityModelResult object
    fit_linear = function(data, group_id = NULL) {
      if (nrow(data) < 3) {
        return(RirVelocityModelResult$new(
          model_type = "linear",
          r_squared = NA_real_,
          rse = NA_real_,
          coefficients = NA,
          n = nrow(data),
          group_id = group_id
        ))
      }

      model <- stats::lm(rir ~ mean_velocity, data = data)
      summary_model <- summary(model)

      RirVelocityModelResult$new(
        model_type = "linear",
        r_squared = summary_model$r.squared,
        rse = summary_model$sigma,
        coefficients = stats::coef(model),
        n = nrow(data),
        group_id = group_id
      )
    },

    #' @description Fit polynomial (quadratic) model
    #' @param data Data frame with rir and mean_velocity columns
    #' @param group_id Optional group identifier
    #' @return RirVelocityModelResult object
    fit_polynomial = function(data, group_id = NULL) {
      if (nrow(data) < 4) {
        return(RirVelocityModelResult$new(
          model_type = "polynomial",
          r_squared = NA_real_,
          rse = NA_real_,
          coefficients = NA,
          n = nrow(data),
          group_id = group_id
        ))
      }

      model <- stats::lm(rir ~ stats::poly(mean_velocity, 2, raw = TRUE), data = data)
      summary_model <- summary(model)

      RirVelocityModelResult$new(
        model_type = "polynomial",
        r_squared = summary_model$r.squared,
        rse = summary_model$sigma,
        coefficients = stats::coef(model),
        n = nrow(data),
        group_id = group_id
      )
    },

    #' @description Fit both linear and polynomial models
    #' @param data Data frame with rir and mean_velocity columns
    #' @param group_id Optional group identifier
    #' @return Named list with linear and polynomial results
    fit_both = function(data, group_id = NULL) {
      list(
        linear = self$fit_linear(data, group_id),
        polynomial = self$fit_polynomial(data, group_id)
      )
    },

    #' @description Fit general (pooled) models by load
    #' @param data Full data frame
    #' @param load_col Name of load/set_type column
    #' @return Named list of results by load
    fit_general_by_load = function(data, load_col = "set_type") {
      loads <- unique(data[[load_col]])
      results <- list()

      for (load in loads) {
        subset_data <- data[data[[load_col]] == load, ]
        results[[as.character(load)]] <- self$fit_both(subset_data, group_id = load)
      }

      results
    },

    #' @description Fit individual models by participant and load
    #' @param data Full data frame
    #' @param id_col Name of participant ID column
    #' @param load_col Name of load/set_type column
    #' @return Nested list of results
    fit_individual = function(data, id_col = "id", load_col = "set_type") {
      ids <- unique(data[[id_col]])
      loads <- unique(data[[load_col]])
      results <- list()

      for (id in ids) {
        results[[as.character(id)]] <- list()
        for (load in loads) {
          subset_data <- data[data[[id_col]] == id & data[[load_col]] == load, ]
          if (nrow(subset_data) >= 3) {
            results[[as.character(id)]][[as.character(load)]] <- self$fit_both(
              subset_data,
              group_id = paste(id, load, sep = "_")
            )
          }
        }
      }

      results
    },

    #' @description Summarize individual model results
    #' @param individual_results Results from fit_individual
    #' @param model_type "linear" or "polynomial"
    #' @return RirVelocityModelSummary object
    summarize_individual = function(individual_results, model_type = "polynomial") {
      all_results <- list()

      for (id in names(individual_results)) {
        for (load in names(individual_results[[id]])) {
          if (!is.null(individual_results[[id]][[load]][[model_type]])) {
            all_results <- c(all_results, list(individual_results[[id]][[load]][[model_type]]))
          }
        }
      }

      RirVelocityModelSummary$new(all_results, model_type)
    },

    #' @description Calculate prediction accuracy
    #' @param train_data Training data (e.g., Day 1)
    #' @param test_data Test data (e.g., Day 2)
    #' @param model_type "linear" or "polynomial"
    #' @return Data frame with actual RIR, predicted RIR, and error
    calculate_prediction_accuracy = function(train_data, test_data, model_type = "polynomial") {
      if (model_type == "polynomial") {
        model <- stats::lm(rir ~ stats::poly(mean_velocity, 2, raw = TRUE), data = train_data)
      } else {
        model <- stats::lm(rir ~ mean_velocity, data = train_data)
      }

      predictions <- stats::predict(model, newdata = test_data)

      data.frame(
        actual_rir = test_data$rir,
        predicted_rir = predictions,
        error = predictions - test_data$rir,
        absolute_error = abs(predictions - test_data$rir)
      )
    }
  )
)
