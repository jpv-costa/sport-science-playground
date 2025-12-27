# R/calculators/conformal_predictor.R
# Service: Conformal Prediction for Distribution-Free Intervals
#
# =============================================================================
# PRINCIPLES APPLIED (from CLAUDE.md)
# =============================================================================
#
# NAMING PRINCIPLES:
# - Understandability: Domain terms (conformal, calibration, coverage)
# - Consistency: All public methods verb-based (fit, predict, calibrate)
# - Distinguishability: predict_interval vs predict_point
# - Conciseness: cp for Conformal Prediction, ci for Confidence Interval
#
# SOLID PRINCIPLES:
# - SRP: Single responsibility - generate distribution-free prediction intervals
# - OCP: Open for extension with different conformity score functions
# - DIP: Depends on model abstractions (predict function), not implementations
#
# CUPID PRINCIPLES:
# - Composable: fit -> calibrate -> predict pipeline
# - Unix: Each method does one thing (fit, calibrate, predict)
# - Predictable: Given calibration data -> deterministic intervals
# - Idiomatic: Follows R conventions for prediction methods
# - Domain-based: Names reflect statistical concepts (conformal, nonconformity)
#
# SCIENTIFIC VALIDITY:
# - Split conformal prediction: Vovk et al. (2005)
# - Guaranteed coverage: P(Y ∈ C(X)) ≥ 1-α for any distribution
# - Nonconformity score: |y - ŷ| (absolute residual)
# - Calibration: Use held-out calibration set to determine quantile
#
# TESTABILITY:
# - Coverage guarantees can be verified empirically
# - Calibration is deterministic given residuals
# =============================================================================

box::use(
  R6[R6Class]
)

#' Conformal Prediction Result
#'
#' Value object containing prediction with conformal interval.
#'
#' @export
ConformalPredictionResult <- R6Class(

  classname = "ConformalPredictionResult",
  cloneable = FALSE,

  public = list(

    #' @description Create conformal prediction result
    #' @param predictions Data frame with point predictions and intervals
    #' @param coverage_target Target coverage (1 - alpha)
    #' @param q_hat Calibrated quantile of nonconformity scores
    #' @param n_calibration Number of calibration observations
    initialize = function(predictions, coverage_target, q_hat, n_calibration) {
      private$.predictions <- predictions
      private$.coverage_target <- coverage_target
      private$.q_hat <- q_hat
      private$.n_calibration <- n_calibration
    },

    #' @description Convert to list
    to_list = function() {
      list(
        predictions = private$.predictions,
        coverage_target = private$.coverage_target,
        q_hat = private$.q_hat,
        n_calibration = private$.n_calibration
      )
    }
  ),

  active = list(
    predictions = function() private$.predictions,
    coverage_target = function() private$.coverage_target,
    q_hat = function() private$.q_hat,
    n_calibration = function() private$.n_calibration
  ),

  private = list(
    .predictions = NULL,
    .coverage_target = NULL,
    .q_hat = NULL,
    .n_calibration = NULL
  )
)


#' Coverage Comparison Result
#'
#' Value object containing comparison between parametric and conformal intervals.
#'
#' @export
CoverageComparisonResult <- R6Class(

  classname = "CoverageComparisonResult",
  cloneable = FALSE,

  public = list(

    #' @description Create coverage comparison result
    #' @param parametric_coverage Empirical coverage of parametric intervals
    #' @param conformal_coverage Empirical coverage of conformal intervals
    #' @param parametric_width Mean width of parametric intervals
    #' @param conformal_width Mean width of conformal intervals
    #' @param target_coverage Target coverage level
    initialize = function(parametric_coverage,
                          conformal_coverage,
                          parametric_width,
                          conformal_width,
                          target_coverage) {
      private$.parametric_coverage <- parametric_coverage
      private$.conformal_coverage <- conformal_coverage
      private$.parametric_width <- parametric_width
      private$.conformal_width <- conformal_width
      private$.target_coverage <- target_coverage
    },

    #' @description Convert to list
    to_list = function() {
      list(
        parametric_coverage = private$.parametric_coverage,
        conformal_coverage = private$.conformal_coverage,
        parametric_width = private$.parametric_width,
        conformal_width = private$.conformal_width,
        target_coverage = private$.target_coverage,
        parametric_coverage_error = abs(private$.parametric_coverage - private$.target_coverage),
        conformal_coverage_error = abs(private$.conformal_coverage - private$.target_coverage)
      )
    }
  ),

  active = list(
    parametric_coverage = function() private$.parametric_coverage,
    conformal_coverage = function() private$.conformal_coverage,
    parametric_width = function() private$.parametric_width,
    conformal_width = function() private$.conformal_width,
    target_coverage = function() private$.target_coverage
  ),

  private = list(
    .parametric_coverage = NULL,
    .conformal_coverage = NULL,
    .parametric_width = NULL,
    .conformal_width = NULL,
    .target_coverage = NULL
  )
)


#' Conformal Predictor
#'
#' Implements split conformal prediction for distribution-free prediction intervals.
#' Uses absolute residual as the nonconformity score.
#'
#' @references
#' Lei, J., G'Sell, M., Rinaldo, A., Tibshirani, R. J., & Wasserman, L. (2018).
#' Distribution-Free Predictive Inference for Regression.
#' Journal of the American Statistical Association, 113(523), 1094-1111.
#'
#' @export
ConformalPredictor <- R6Class(

  classname = "ConformalPredictor",

  public = list(

    #' @description Create conformal predictor
    #' @param alpha Miscoverage rate (default 0.05 for 95% intervals)
    initialize = function(alpha = 0.05) {
      private$.alpha <- alpha
      private$.fitted <- FALSE
    },

    #' @description Fit conformal predictor using calibration data
    #' @param model Fitted model object (must have predict method)
    #' @param calibration_data Data frame for calibration
    #' @param response_var Name of the response variable
    fit = function(model, calibration_data, response_var = "mean_velocity") {
      private$.model <- model
      private$.response_var <- response_var

      # Get predictions on calibration set
      calibration_predictions <- private$.get_predictions(calibration_data)

      # Calculate nonconformity scores (absolute residuals)
      actual <- calibration_data[[response_var]]
      private$.nonconformity_scores <- abs(actual - calibration_predictions)

      # Calculate calibrated quantile
      # Use (1 - alpha)(1 + 1/n) quantile for finite sample guarantee
      n <- length(private$.nonconformity_scores)
      quantile_level <- (1 - private$.alpha) * (1 + 1 / n)
      quantile_level <- min(quantile_level, 1)  # Cap at 1

      private$.q_hat <- stats::quantile(private$.nonconformity_scores, probs = quantile_level)
      private$.n_calibration <- n
      private$.fitted <- TRUE

      invisible(self)
    },

    #' @description Generate prediction intervals for new data
    #' @param newdata Data frame for prediction
    #' @return ConformalPredictionResult object
    predict_interval = function(newdata) {
      if (!private$.fitted) {
        stop("ConformalPredictor must be fitted before prediction")
      }

      # Get point predictions
      predictions <- private$.get_predictions(newdata)

      # Create conformal intervals
      lower <- predictions - private$.q_hat
      upper <- predictions + private$.q_hat

      result_df <- data.frame(
        prediction = predictions,
        lower = lower,
        upper = upper,
        interval_width = upper - lower
      )

      ConformalPredictionResult$new(
        predictions = result_df,
        coverage_target = 1 - private$.alpha,
        q_hat = private$.q_hat,
        n_calibration = private$.n_calibration
      )
    },

    #' @description Calculate empirical coverage on test data
    #' @param test_data Data frame with actual values
    #' @param conformal_result ConformalPredictionResult from predict_interval
    #' @return Numeric coverage rate
    calculate_coverage = function(test_data, conformal_result) {
      actual <- test_data[[private$.response_var]]
      predictions <- conformal_result$predictions

      covered <- (actual >= predictions$lower) & (actual <= predictions$upper)
      mean(covered, na.rm = TRUE)
    },

    #' @description Compare conformal with parametric prediction intervals
    #' @param model_result LmmModelResult object for parametric intervals
    #' @param test_data Data frame for testing
    #' @param parametric_level Confidence level for parametric intervals (default 0.95)
    #' @return CoverageComparisonResult object
    compare_with_parametric = function(model_result, test_data, parametric_level = 0.95) {
      if (!private$.fitted) {
        stop("ConformalPredictor must be fitted before comparison")
      }

      # Get conformal intervals
      conformal_result <- self$predict_interval(test_data)

      # Get parametric intervals
      parametric_result <- private$.get_parametric_intervals(
        model_result$model,
        test_data,
        parametric_level
      )

      # Calculate coverages
      actual <- test_data[[private$.response_var]]

      conformal_coverage <- self$calculate_coverage(test_data, conformal_result)

      parametric_covered <- (actual >= parametric_result$lower) &
        (actual <= parametric_result$upper)
      parametric_coverage <- mean(parametric_covered, na.rm = TRUE)

      # Calculate widths
      conformal_width <- mean(conformal_result$predictions$interval_width, na.rm = TRUE)
      parametric_width <- mean(parametric_result$upper - parametric_result$lower, na.rm = TRUE)

      CoverageComparisonResult$new(
        parametric_coverage = parametric_coverage,
        conformal_coverage = conformal_coverage,
        parametric_width = parametric_width,
        conformal_width = conformal_width,
        target_coverage = 1 - private$.alpha
      )
    },

    #' @description Get distribution of interval widths
    #' @param conformal_result ConformalPredictionResult object
    #' @return List with summary statistics of interval widths
    get_interval_width_distribution = function(conformal_result) {
      widths <- conformal_result$predictions$interval_width

      list(
        mean = mean(widths, na.rm = TRUE),
        median = stats::median(widths, na.rm = TRUE),
        sd = stats::sd(widths, na.rm = TRUE),
        min = min(widths, na.rm = TRUE),
        max = max(widths, na.rm = TRUE),
        q25 = stats::quantile(widths, 0.25, na.rm = TRUE),
        q75 = stats::quantile(widths, 0.75, na.rm = TRUE)
      )
    },

    #' @description Get the calibrated quantile
    #' @return Numeric q_hat value
    get_q_hat = function() {
      if (!private$.fitted) {
        stop("ConformalPredictor must be fitted first")
      }
      private$.q_hat
    },

    #' @description Get nonconformity scores from calibration
    #' @return Vector of nonconformity scores
    get_nonconformity_scores = function() {
      if (!private$.fitted) {
        stop("ConformalPredictor must be fitted first")
      }
      private$.nonconformity_scores
    }
  ),

  private = list(

    .alpha = 0.05,
    .model = NULL,
    .response_var = NULL,
    .q_hat = NULL,
    .n_calibration = NULL,
    .nonconformity_scores = NULL,
    .fitted = FALSE,

    .get_predictions = function(data) {
      # Handle different model types
      model <- private$.model

      if (inherits(model, "lmerMod")) {
        # LMM model - predict with random effects
        stats::predict(model, newdata = data, allow.new.levels = TRUE)
      } else if (inherits(model, "LmmModelResult")) {
        # Our wrapper class
        stats::predict(model$model, newdata = data, allow.new.levels = TRUE)
      } else if (inherits(model, "lm")) {
        # Simple linear model
        stats::predict(model, newdata = data)
      } else {
        # Try generic predict
        stats::predict(model, newdata = data)
      }
    },

    .get_parametric_intervals = function(model, newdata, level) {
      # Get parametric prediction intervals from LMM
      if (inherits(model, "lmerMod")) {
        # For LMM, approximate using fixed effects SE + residual SE
        pred <- stats::predict(model, newdata = newdata, re.form = NA)
        residual_se <- stats::sigma(model)

        # Use normal approximation
        z <- stats::qnorm((1 + level) / 2)

        data.frame(
          prediction = pred,
          lower = pred - z * residual_se,
          upper = pred + z * residual_se
        )
      } else if (inherits(model, "lm")) {
        # For lm, use built-in prediction intervals
        pred <- stats::predict(model, newdata = newdata, interval = "prediction", level = level)
        data.frame(
          prediction = pred[, "fit"],
          lower = pred[, "lwr"],
          upper = pred[, "upr"]
        )
      } else {
        stop("Unsupported model type for parametric intervals")
      }
    }
  )
)
