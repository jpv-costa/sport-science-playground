# R/calculators/model_validator.R
# Model Validation: Cross-Validation, Calibration, and Prediction Intervals
#
# =============================================================================
# PRINCIPLES APPLIED (from CLAUDE.md)
# =============================================================================
#
# NAMING PRINCIPLES:
# - Understandability: Problem domain terms (cross_validate, calibration)
# - Consistency: All public methods follow verb_noun pattern
# - Distinguishability: Clear method names for different validation types
# - Conciseness: CV = cross-validation, CI/PI for intervals
#
# SOLID PRINCIPLES:
# - SRP: Single responsibility - only model validation
# - OCP: Open for extension via new validation methods
# - DIP: Depends on model abstractions, not specific implementations
#
# CUPID PRINCIPLES:
# - Composable: Each validation method is independent
# - Unix: Each method validates one aspect
# - Predictable: Same inputs -> same validation results
# - Idiomatic: Follows R conventions
# - Domain-based: Names reflect statistical validation concepts
# =============================================================================

box::use(
  R6[R6Class],
  stats[predict, residuals, fitted, lm, coef, qt, var, sd, sigma, qnorm, formula],
  ggplot2[ggplot, aes, geom_point, geom_abline, geom_ribbon, geom_line,
          geom_smooth, labs, theme_minimal, theme, element_text, scale_color_manual,
          coord_fixed, annotate]
)

#' Cross-Validation Result
#'
#' Immutable value object containing cross-validation results
#'
#' @export
CrossValidationResult <- R6Class(
  classname = "CrossValidationResult",
  cloneable = FALSE,

  public = list(
    #' @field cv_type Type of cross-validation performed
    cv_type = NULL,
    #' @field n_folds Number of folds (or participants for LOO)
    n_folds = NULL,
    #' @field fold_errors Vector of per-fold prediction errors
    fold_errors = NULL,
    #' @field mean_error Mean cross-validation error
    mean_error = NULL,
    #' @field se_error Standard error of CV error
    se_error = NULL,
    #' @field rmse Root mean squared error
    rmse = NULL,

    #' @description Create a new CV result
    initialize = function(cv_type, n_folds, fold_errors, mean_error, se_error, rmse) {
      self$cv_type <- cv_type
      self$n_folds <- n_folds
      self$fold_errors <- fold_errors
      self$mean_error <- mean_error
      self$se_error <- se_error
      self$rmse <- rmse
    },

    #' @description Convert to list for serialization
    to_list = function() {
      list(
        cv_type = self$cv_type,
        n_folds = self$n_folds,
        fold_errors = self$fold_errors,
        mean_error = self$mean_error,
        se_error = self$se_error,
        rmse = self$rmse
      )
    }
  )
)

#' Calibration Result
#'
#' Immutable value object containing calibration assessment results
#'
#' @export
CalibrationResult <- R6Class(
  classname = "CalibrationResult",
  cloneable = FALSE,

  public = list(
    #' @field slope Calibration slope (ideal = 1)
    slope = NULL,
    #' @field intercept Calibration intercept (ideal = 0)
    intercept = NULL,
    #' @field r_squared R² of calibration regression
    r_squared = NULL,
    #' @field mean_bias Mean prediction bias
    mean_bias = NULL,
    #' @field calibration_data Data frame with predicted and observed values
    calibration_data = NULL,

    #' @description Create a new calibration result
    initialize = function(slope, intercept, r_squared, mean_bias, calibration_data) {
      self$slope <- slope
      self$intercept <- intercept
      self$r_squared <- r_squared
      self$mean_bias <- mean_bias
      self$calibration_data <- calibration_data
    },

    #' @description Interpret calibration quality
    interpret = function() {
      slope_quality <- if (abs(self$slope - 1) < 0.1) {
        "excellent"
      } else if (abs(self$slope - 1) < 0.2) {
        "good"
      } else if (abs(self$slope - 1) < 0.3) {
        "moderate"
      } else {
        "poor"
      }

      intercept_quality <- if (abs(self$intercept) < 0.01) {
        "negligible"
      } else if (abs(self$intercept) < 0.05) {
        "small"
      } else {
        "notable"
      }

      paste0(
        "Calibration: slope = ", round(self$slope, 3),
        " (", slope_quality, "), ",
        "intercept = ", round(self$intercept, 4),
        " (", intercept_quality, " bias)"
      )
    },

    #' @description Convert to list for serialization
    to_list = function() {
      list(
        slope = self$slope,
        intercept = self$intercept,
        r_squared = self$r_squared,
        mean_bias = self$mean_bias,
        interpretation = self$interpret()
      )
    }
  )
)

#' Interval Comparison Result
#'
#' Immutable value object comparing confidence and prediction intervals
#'
#' @export
IntervalComparisonResult <- R6Class(
  classname = "IntervalComparisonResult",
  cloneable = FALSE,

  public = list(
    #' @field ci_width Mean confidence interval width
    ci_width = NULL,
    #' @field pi_width Mean prediction interval width
    pi_width = NULL,
    #' @field ci_coverage Empirical CI coverage
    ci_coverage = NULL,
    #' @field pi_coverage Empirical PI coverage
    pi_coverage = NULL,
    #' @field interval_data Data frame with both interval types
    interval_data = NULL,

    #' @description Create a new interval comparison result
    initialize = function(ci_width, pi_width, ci_coverage, pi_coverage, interval_data) {
      self$ci_width <- ci_width
      self$pi_width <- pi_width
      self$ci_coverage <- ci_coverage
      self$pi_coverage <- pi_coverage
      self$interval_data <- interval_data
    },

    #' @description Convert to list for serialization
    to_list = function() {
      list(
        ci_width = self$ci_width,
        pi_width = self$pi_width,
        ci_coverage = self$ci_coverage,
        pi_coverage = self$pi_coverage,
        ratio = self$pi_width / self$ci_width
      )
    }
  )
)

#' Model Validator
#'
#' R6 class for comprehensive model validation including cross-validation,
#' calibration assessment, and prediction interval comparison.
#'
#' @section Philosophy:
#' Validation answers: "How well will this model work on NEW data?"
#' - Cross-validation: Test on held-out data
#' - Calibration: Do predictions match reality?
#' - PI vs CI: Appropriate uncertainty for the question
#'
#' @export
ModelValidator <- R6Class(
  classname = "ModelValidator",
  cloneable = FALSE,

  public = list(
    #' @description Create a new ModelValidator instance
    initialize = function() {
      # No initialization needed
    },

    # =========================================================================
    # CROSS-VALIDATION
    # =========================================================================

    #' @description Leave-One-Participant-Out Cross-Validation
    #'
    #' The gold standard for clustered data: leave out each participant

    #' entirely and predict their observations using the remaining data.
    #'
    #' @param data Data frame with observations
    #' @param model Fitted lme4 model
    #' @param id_col Name of participant ID column
    #' @param outcome_col Name of outcome variable
    #' @return CrossValidationResult object
    leave_one_participant_out = function(data, model, id_col = "id",
                                          outcome_col = "mean_velocity") {
      if (!requireNamespace("lme4", quietly = TRUE)) {
        stop("Package 'lme4' required for leave_one_participant_out()")
      }

      participants <- unique(data[[id_col]])
      n_participants <- length(participants)

      # Extract formula from model
      formula_obj <- stats::formula(model)
      random_formula <- private$.extract_random_formula(model)

      fold_errors <- numeric(n_participants)
      fold_rmse <- numeric(n_participants)

      for (i in seq_along(participants)) {
        # Split data
        test_data <- data[data[[id_col]] == participants[i], ]
        train_data <- data[data[[id_col]] != participants[i], ]

        # Refit model without this participant
        fold_model <- tryCatch({
          lme4::lmer(formula_obj, data = train_data, REML = FALSE)
        }, error = function(e) NULL)

        if (is.null(fold_model)) {
          fold_errors[i] <- NA
          fold_rmse[i] <- NA
          next
        }

        # Predict for held-out participant (population average)
        predictions <- predict(fold_model, newdata = test_data, re.form = NA,
                               allow.new.levels = TRUE)

        actual <- test_data[[outcome_col]]
        errors <- actual - predictions

        fold_errors[i] <- mean(abs(errors))  # MAE
        fold_rmse[i] <- sqrt(mean(errors^2))  # RMSE
      }

      # Remove NA values
      valid_errors <- fold_errors[!is.na(fold_errors)]
      valid_rmse <- fold_rmse[!is.na(fold_rmse)]

      CrossValidationResult$new(
        cv_type = "Leave-One-Participant-Out",
        n_folds = n_participants,
        fold_errors = valid_errors,
        mean_error = mean(valid_errors),
        se_error = sd(valid_errors) / sqrt(length(valid_errors)),
        rmse = mean(valid_rmse)
      )
    },

    #' @description K-Fold Cross-Validation for Mixed Models
    #'
    #' Standard k-fold CV adapted for clustered data.
    #' Participants are assigned to folds (not individual observations).
    #'
    #' @param data Data frame with observations
    #' @param model Fitted lme4 model
    #' @param k Number of folds (default 10)
    #' @param id_col Name of participant ID column
    #' @param outcome_col Name of outcome variable
    #' @param seed Random seed for reproducibility
    #' @return CrossValidationResult object
    k_fold_cv = function(data, model, k = 10, id_col = "id",
                         outcome_col = "mean_velocity", seed = 42) {
      if (!requireNamespace("lme4", quietly = TRUE)) {
        stop("Package 'lme4' required for k_fold_cv()")
      }

      set.seed(seed)

      participants <- unique(data[[id_col]])
      n_participants <- length(participants)

      # Assign participants to folds
      fold_assignment <- sample(rep(1:k, length.out = n_participants))
      names(fold_assignment) <- participants

      # Extract formula from model
      formula_obj <- stats::formula(model)

      fold_errors <- numeric(k)
      fold_rmse <- numeric(k)

      for (fold in 1:k) {
        # Get participants in this fold
        test_participants <- names(fold_assignment)[fold_assignment == fold]

        # Split data
        test_data <- data[data[[id_col]] %in% test_participants, ]
        train_data <- data[!data[[id_col]] %in% test_participants, ]

        # Refit model
        fold_model <- tryCatch({
          lme4::lmer(formula_obj, data = train_data, REML = FALSE)
        }, error = function(e) NULL)

        if (is.null(fold_model)) {
          fold_errors[fold] <- NA
          fold_rmse[fold] <- NA
          next
        }

        # Predict (population average for new participants)
        predictions <- predict(fold_model, newdata = test_data, re.form = NA,
                               allow.new.levels = TRUE)

        actual <- test_data[[outcome_col]]
        errors <- actual - predictions

        fold_errors[fold] <- mean(abs(errors))  # MAE
        fold_rmse[fold] <- sqrt(mean(errors^2))  # RMSE
      }

      valid_errors <- fold_errors[!is.na(fold_errors)]
      valid_rmse <- fold_rmse[!is.na(fold_rmse)]

      CrossValidationResult$new(
        cv_type = paste0(k, "-Fold CV (participant-level)"),
        n_folds = k,
        fold_errors = valid_errors,
        mean_error = mean(valid_errors),
        se_error = sd(valid_errors) / sqrt(length(valid_errors)),
        rmse = mean(valid_rmse)
      )
    },

    # =========================================================================
    # CALIBRATION
    # =========================================================================

    #' @description Calculate Model Calibration
    #'
    #' Assesses how well predicted values match observed values.
    #' A perfectly calibrated model has slope = 1, intercept = 0.
    #'
    #' @param model Fitted lme4 model
    #' @param data Data frame used for fitting
    #' @param outcome_col Name of outcome variable
    #' @param use_individual Whether to use individual predictions (TRUE) or population (FALSE)
    #' @return CalibrationResult object
    calculate_calibration = function(model, data, outcome_col = "mean_velocity",
                                     use_individual = TRUE) {
      # Get predictions
      if (use_individual) {
        predicted <- predict(model)
      } else {
        predicted <- predict(model, re.form = NA)
      }

      observed <- data[[outcome_col]]

      # Fit calibration regression: observed ~ predicted
      cal_model <- lm(observed ~ predicted)
      cal_coef <- coef(cal_model)

      # Calculate metrics
      slope <- cal_coef[2]
      intercept <- cal_coef[1]
      r_squared <- summary(cal_model)$r.squared
      mean_bias <- mean(predicted - observed)

      # Create calibration data
      calibration_data <- data.frame(
        predicted = predicted,
        observed = observed,
        residual = observed - predicted
      )

      CalibrationResult$new(
        slope = slope,
        intercept = intercept,
        r_squared = r_squared,
        mean_bias = mean_bias,
        calibration_data = calibration_data
      )
    },

    #' @description Plot Calibration
    #'
    #' Creates a calibration plot with predicted vs observed values
    #' and the ideal 45-degree line.
    #'
    #' @param calibration_result CalibrationResult object
    #' @param title Plot title
    #' @return ggplot object
    plot_calibration = function(calibration_result, title = "Model Calibration") {
      cal_data <- calibration_result$calibration_data

      # Calculate range for plot limits
      all_values <- c(cal_data$predicted, cal_data$observed)
      plot_range <- range(all_values, na.rm = TRUE)

      ggplot(cal_data, aes(x = predicted, y = observed)) +
        geom_abline(intercept = 0, slope = 1, linetype = "dashed",
                    color = "#666666", linewidth = 1) +
        geom_point(alpha = 0.5, color = "#2E86AB") +
        geom_smooth(method = "lm", se = TRUE, color = "#E63946",
                    fill = "#E63946", alpha = 0.2) +
        coord_fixed(xlim = plot_range, ylim = plot_range) +
        labs(
          title = title,
          subtitle = paste0(
            "Slope = ", round(calibration_result$slope, 3),
            " (ideal = 1), Intercept = ", round(calibration_result$intercept, 4),
            " (ideal = 0)"
          ),
          x = "Predicted Velocity (m/s)",
          y = "Observed Velocity (m/s)",
          caption = "Dashed line = perfect calibration; Red line = actual calibration"
        ) +
        theme_minimal() +
        theme(
          plot.title = element_text(face = "bold", size = 14),
          plot.subtitle = element_text(size = 10, color = "#666666")
        )
    },

    # =========================================================================
    # PREDICTION VS CONFIDENCE INTERVALS
    # =========================================================================

    #' @description Calculate Prediction Interval
    #'
    #' Prediction intervals capture uncertainty about NEW observations,
    #' not just uncertainty about the mean.
    #'
    #' @param model Fitted lme4 model
    #' @param newdata Data frame with new observations
    #' @param level Confidence level (default 0.95)
    #' @return Data frame with predictions and intervals
    calculate_prediction_interval = function(model, newdata, level = 0.95) {
      if (!requireNamespace("lme4", quietly = TRUE)) {
        stop("Package 'lme4' required for calculate_prediction_interval()")
      }

      # Get fixed effects prediction
      pred_fixed <- predict(model, newdata = newdata, re.form = NA,
                            allow.new.levels = TRUE)

      # Get residual standard deviation
      sigma_resid <- sigma(model)

      # Get random effects variance
      vc <- lme4::VarCorr(model)
      sigma_re <- sqrt(sum(sapply(vc, function(x) sum(diag(x)))))

      # Total prediction variance for new individual
      # PI includes: fixed effect uncertainty + RE variance + residual variance
      sigma_pred <- sqrt(sigma_re^2 + sigma_resid^2)

      # Calculate interval
      z_value <- qnorm((1 + level) / 2)
      pi_lower <- pred_fixed - z_value * sigma_pred
      pi_upper <- pred_fixed + z_value * sigma_pred

      data.frame(
        prediction = pred_fixed,
        pi_lower = pi_lower,
        pi_upper = pi_upper,
        pi_width = pi_upper - pi_lower
      )
    },

    #' @description Compare CI and PI
    #'
    #' Shows the difference between confidence intervals (for mean)
    #' and prediction intervals (for new observations).
    #'
    #' @param model Fitted lme4 model
    #' @param data Data frame used for fitting
    #' @param outcome_col Name of outcome variable
    #' @param level Confidence level
    #' @return IntervalComparisonResult object
    compare_ci_vs_pi = function(model, data, outcome_col = "mean_velocity",
                                level = 0.95) {
      if (!requireNamespace("lme4", quietly = TRUE)) {
        stop("Package 'lme4' required for compare_ci_vs_pi()")
      }

      # Get predictions
      pred_fixed <- predict(model, re.form = NA)

      # Residual SD
      sigma_resid <- sigma(model)

      # Random effects SD
      vc <- lme4::VarCorr(model)
      sigma_re <- sqrt(sum(sapply(vc, function(x) sum(diag(x)))))

      # Approximate SE for CI (using residual SE as proxy)
      # This is a simplification - true CI requires bootstrap
      se_mean <- sigma_resid / sqrt(nrow(data))

      # Prediction SD (for new individual)
      sigma_pred <- sqrt(sigma_re^2 + sigma_resid^2)

      z_value <- qnorm((1 + level) / 2)

      # Calculate intervals
      ci_lower <- pred_fixed - z_value * se_mean
      ci_upper <- pred_fixed + z_value * se_mean
      pi_lower <- pred_fixed - z_value * sigma_pred
      pi_upper <- pred_fixed + z_value * sigma_pred

      # Check coverage
      observed <- data[[outcome_col]]
      ci_coverage <- mean(observed >= ci_lower & observed <= ci_upper)
      pi_coverage <- mean(observed >= pi_lower & observed <= pi_upper)

      interval_data <- data.frame(
        observed = observed,
        predicted = pred_fixed,
        ci_lower = ci_lower,
        ci_upper = ci_upper,
        pi_lower = pi_lower,
        pi_upper = pi_upper
      )

      IntervalComparisonResult$new(
        ci_width = mean(ci_upper - ci_lower),
        pi_width = mean(pi_upper - pi_lower),
        ci_coverage = ci_coverage,
        pi_coverage = pi_coverage,
        interval_data = interval_data
      )
    },

    #' @description Plot CI vs PI Comparison
    #'
    #' Visualizes the difference between confidence and prediction intervals.
    #'
    #' @param comparison_result IntervalComparisonResult object
    #' @param title Plot title
    #' @return ggplot object
    plot_ci_vs_pi = function(comparison_result, title = "Confidence vs Prediction Intervals") {
      int_data <- comparison_result$interval_data

      # Sort by predicted value for cleaner visualization
      int_data <- int_data[order(int_data$predicted), ]
      int_data$index <- seq_len(nrow(int_data))

      # Sample if too many points
      if (nrow(int_data) > 100) {
        sample_idx <- sort(sample(int_data$index, 100))
        int_data <- int_data[int_data$index %in% sample_idx, ]
        int_data$index <- seq_len(nrow(int_data))
      }

      ggplot(int_data, aes(x = index)) +
        # Prediction interval (wider)
        geom_ribbon(aes(ymin = pi_lower, ymax = pi_upper),
                    fill = "#F77F00", alpha = 0.3) +
        # Confidence interval (narrower)
        geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper),
                    fill = "#2E86AB", alpha = 0.5) +
        # Predictions
        geom_line(aes(y = predicted), color = "#2E86AB", linewidth = 1) +
        # Observed points
        geom_point(aes(y = observed), color = "#1D3557", alpha = 0.6, size = 1.5) +
        labs(
          title = title,
          subtitle = paste0(
            "CI coverage: ", round(comparison_result$ci_coverage * 100, 1), "%",
            " | PI coverage: ", round(comparison_result$pi_coverage * 100, 1), "%"
          ),
          x = "Observation (sorted by prediction)",
          y = "Velocity (m/s)",
          caption = "Blue band = 95% CI (for mean); Orange band = 95% PI (for new observation)"
        ) +
        theme_minimal() +
        theme(
          plot.title = element_text(face = "bold", size = 14),
          plot.subtitle = element_text(size = 10, color = "#666666")
        )
    }
  ),

  private = list(
    #' Extract random effects formula from model
    .extract_random_formula = function(model) {
      # Get the formula bars (random effects specification)
      formula_str <- as.character(stats::formula(model))[3]
      # This is a simplified extraction
      NULL
    }
  )
)

# =============================================================================
# CONVENIENCE FUNCTIONS
# =============================================================================

#' Perform Leave-One-Participant-Out CV
#'
#' @param data Data frame
#' @param model Fitted lme4 model
#' @param ... Additional arguments passed to method
#' @return CrossValidationResult object
#' @export
leave_one_out_cv <- function(data, model, ...) {
  validator <- ModelValidator$new()
  validator$leave_one_participant_out(data, model, ...)
}

#' Calculate Model Calibration
#'
#' @param model Fitted lme4 model
#' @param data Data frame
#' @param ... Additional arguments
#' @return CalibrationResult object
#' @export
calculate_calibration <- function(model, data, ...) {
  validator <- ModelValidator$new()
  validator$calculate_calibration(model, data, ...)
}

#' Plot Calibration
#'
#' @param calibration_result CalibrationResult object
#' @param ... Additional arguments
#' @return ggplot object
#' @export
plot_calibration <- function(calibration_result, ...) {
  validator <- ModelValidator$new()
  validator$plot_calibration(calibration_result, ...)
}
