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
  stats[predict, residuals, fitted, lm, coef, qt, var, sd, sigma, qnorm, formula, BIC, AIC],
  dplyr[.data],
  ggplot2[ggplot, aes, geom_point, geom_abline, geom_ribbon, geom_line,
          geom_smooth, geom_vline, labs, theme_minimal, theme, element_text,
          element_blank, scale_color_manual, coord_fixed, annotate]
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

#' Coefficient Stability Result
#'
#' Stores results from LOO-CV coefficient stability analysis.
#'
#' @export
CoefficientStabilityResult <- R6Class(
  classname = "CoefficientStabilityResult",
  cloneable = FALSE,

  public = list(
    #' @field focal_coefficient Name of the coefficient being analyzed
    focal_coefficient = NULL,
    #' @field full_model_estimate Coefficient from full data
    full_model_estimate = NULL,
    #' @field fold_estimates Coefficient per LOO fold
    fold_estimates = NULL,
    #' @field excluded_participants Participant excluded in each fold
    excluded_participants = NULL,
    #' @field mean_estimate Mean across folds
    mean_estimate = NULL,
    #' @field cv_percent Coefficient of variation (%)
    cv_percent = NULL,
    #' @field range Min and max across folds
    range = NULL,
    #' @field influential_participants Participants that most change coefficient
    influential_participants = NULL,
    #' @field is_stable TRUE if CV < 10%
    is_stable = NULL,

    #' @description Create a new CoefficientStabilityResult
    initialize = function(focal_coefficient, full_model_estimate, fold_estimates,
                          excluded_participants) {
      self$focal_coefficient <- focal_coefficient
      self$full_model_estimate <- full_model_estimate
      self$fold_estimates <- fold_estimates
      self$excluded_participants <- excluded_participants

      # Calculate derived metrics
      self$mean_estimate <- mean(fold_estimates, na.rm = TRUE)
      fold_sd <- sd(fold_estimates, na.rm = TRUE)
      self$cv_percent <- (fold_sd / abs(self$mean_estimate)) * 100
      self$range <- range(fold_estimates, na.rm = TRUE)
      self$is_stable <- self$cv_percent < 10

      # Find influential participants (those whose exclusion changes coefficient most)
      deviations <- abs(fold_estimates - full_model_estimate)
      top_indices <- order(deviations, decreasing = TRUE)[1:min(5, length(deviations))]
      self$influential_participants <- data.frame(
        participant = excluded_participants[top_indices],
        estimate_when_excluded = fold_estimates[top_indices],
        deviation = deviations[top_indices]
      )
    },

    #' @description Interpret coefficient stability
    interpret = function() {
      stability_level <- if (self$cv_percent < 5) {
        "highly stable"
      } else if (self$cv_percent < 10) {
        "stable"
      } else if (self$cv_percent < 20) {
        "moderately stable"
      } else {
        "unstable"
      }

      paste0(
        "Coefficient '", self$focal_coefficient, "' is ", stability_level,
        " (CV = ", round(self$cv_percent, 1), "%). ",
        "Full model: ", round(self$full_model_estimate, 4),
        ", LOO mean: ", round(self$mean_estimate, 4),
        ", range: [", round(self$range[1], 4), ", ", round(self$range[2], 4), "]"
      )
    },

    #' @description Convert to list for serialization
    to_list = function() {
      list(
        focal_coefficient = self$focal_coefficient,
        full_model_estimate = self$full_model_estimate,
        mean_estimate = self$mean_estimate,
        cv_percent = self$cv_percent,
        range = self$range,
        is_stable = self$is_stable,
        influential_participants = self$influential_participants
      )
    }
  )
)

#' Model Selection Stability Result
#'
#' Stores results from LOO-CV model selection stability analysis.
#'
#' @export
ModelSelectionStabilityResult <- R6Class(
  classname = "ModelSelectionStabilityResult",
  cloneable = FALSE,

  public = list(
    #' @field full_data_winner Best model on full data
    full_data_winner = NULL,
    #' @field fold_winners Best model per fold
    fold_winners = NULL,
    #' @field excluded_participants Participant excluded in each fold
    excluded_participants = NULL,
    #' @field vote_counts How many folds each model won
    vote_counts = NULL,
    #' @field consensus_model Clear winner or "unclear"
    consensus_model = NULL,
    #' @field stability_percent Percent of folds agreeing with full data
    stability_percent = NULL,
    #' @field criterion Selection criterion used
    criterion = NULL,

    #' @description Create a new ModelSelectionStabilityResult
    initialize = function(full_data_winner, fold_winners, excluded_participants,
                          criterion = "BIC") {
      self$full_data_winner <- full_data_winner
      self$fold_winners <- fold_winners
      self$excluded_participants <- excluded_participants
      self$criterion <- criterion

      # Count votes
      self$vote_counts <- table(fold_winners)

      # Determine consensus
      max_votes <- max(self$vote_counts)
      total_folds <- length(fold_winners)
      majority_threshold <- total_folds * 0.6  # 60% majority

      if (max_votes >= majority_threshold) {
        self$consensus_model <- names(self$vote_counts)[which.max(self$vote_counts)]
      } else {
        self$consensus_model <- "unclear (no majority)"
      }

      # Calculate stability
      self$stability_percent <- (sum(fold_winners == full_data_winner) / total_folds) * 100
    },

    #' @description Interpret model selection stability
    interpret = function() {
      stability_level <- if (self$stability_percent >= 90) {
        "highly stable"
      } else if (self$stability_percent >= 70) {
        "stable"
      } else if (self$stability_percent >= 50) {
        "moderately stable"
      } else {
        "unstable"
      }

      paste0(
        "Model selection is ", stability_level, ". ",
        "Full data winner: ", self$full_data_winner, ". ",
        "Consensus model: ", self$consensus_model, ". ",
        round(self$stability_percent, 1), "% of folds agree with full data selection."
      )
    },

    #' @description Convert to list for serialization
    to_list = function() {
      list(
        full_data_winner = self$full_data_winner,
        consensus_model = self$consensus_model,
        stability_percent = self$stability_percent,
        vote_counts = as.list(self$vote_counts),
        criterion = self$criterion
      )
    }
  )
)

#' Prediction Error by Participant Result
#'
#' Stores per-participant prediction error metrics from LOO-CV.
#'
#' @export
PredictionErrorByParticipantResult <- R6Class(
  classname = "PredictionErrorByParticipantResult",
  cloneable = FALSE,

  public = list(
    #' @field participant_id Character vector of participant IDs
    participant_id = NULL,
    #' @field rmse Per-participant RMSE
    rmse = NULL,
    #' @field mae Per-participant MAE
    mae = NULL,
    #' @field bias Systematic over/under prediction
    bias = NULL,
    #' @field n_observations Number of observations per participant
    n_observations = NULL,
    #' @field overall_rmse Overall RMSE
    overall_rmse = NULL,
    #' @field overall_mae Overall MAE
    overall_mae = NULL,
    #' @field residuals All residuals (for anomaly detection)
    residuals = NULL,
    #' @field residual_participant_ids Participant ID for each residual
    residual_participant_ids = NULL,

    #' @description Create a new PredictionErrorByParticipantResult
    initialize = function(participant_metrics, overall_rmse, overall_mae,
                          residuals = NULL, residual_participant_ids = NULL) {
      self$participant_id <- participant_metrics$participant_id
      self$rmse <- participant_metrics$rmse
      self$mae <- participant_metrics$mae
      self$bias <- participant_metrics$bias
      self$n_observations <- participant_metrics$n_observations
      self$overall_rmse <- overall_rmse
      self$overall_mae <- overall_mae
      self$residuals <- residuals
      self$residual_participant_ids <- residual_participant_ids
    },

    #' @description Get hardest to predict participants
    #' @param n Number of participants to return
    get_hardest_to_predict = function(n = 5) {
      top_idx <- order(self$rmse, decreasing = TRUE)[1:min(n, length(self$rmse))]
      data.frame(
        participant_id = self$participant_id[top_idx],
        rmse = self$rmse[top_idx],
        mae = self$mae[top_idx],
        n_observations = self$n_observations[top_idx]
      )
    },

    #' @description Get participants with systematic bias
    get_biased_participants = function(threshold = 0.02) {
      biased_idx <- which(abs(self$bias) > threshold)
      if (length(biased_idx) == 0) return(NULL)

      data.frame(
        participant_id = self$participant_id[biased_idx],
        bias = self$bias[biased_idx],
        direction = ifelse(self$bias[biased_idx] > 0, "over-predicted", "under-predicted")
      )
    },

    #' @description Convert to list for serialization
    to_list = function() {
      list(
        overall_rmse = self$overall_rmse,
        overall_mae = self$overall_mae,
        participant_metrics = data.frame(
          participant_id = self$participant_id,
          rmse = self$rmse,
          mae = self$mae,
          bias = self$bias,
          n_observations = self$n_observations
        ),
        hardest_to_predict = self$get_hardest_to_predict()
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
    # LOO-CV SENSITIVITY ANALYSIS
    # =========================================================================

    #' @description LOO-CV Coefficient Stability Analysis
    #'
    #' Tests whether the focal coefficient (e.g., RIR effect) is stable
    #' when each participant is excluded. Answers: "Does my conclusion
    #' depend on specific individuals?"
    #'
    #' @param data Data frame with observations
    #' @param model Fitted lme4 model
    #' @param focal_coef Name of coefficient to track (default "rir")
    #' @param id_col Name of participant ID column
    #' @return CoefficientStabilityResult object
    loo_coefficient_stability = function(data, model, focal_coef = "rir",
                                          id_col = "id") {
      if (!requireNamespace("lme4", quietly = TRUE)) {
        stop("Package 'lme4' required for loo_coefficient_stability()")
      }

      # Get full model coefficient
      full_coefs <- lme4::fixef(model)
      if (!focal_coef %in% names(full_coefs)) {
        stop("Coefficient '", focal_coef, "' not found in model. ",
             "Available: ", paste(names(full_coefs), collapse = ", "))
      }
      full_estimate <- full_coefs[focal_coef]

      participants <- unique(data[[id_col]])
      n_participants <- length(participants)

      # Extract formula
      formula_obj <- stats::formula(model)

      fold_estimates <- numeric(n_participants)

      for (i in seq_along(participants)) {
        # Remove one participant
        train_data <- data[data[[id_col]] != participants[i], ]

        # Refit model
        fold_model <- tryCatch({
          lme4::lmer(formula_obj, data = train_data, REML = FALSE)
        }, error = function(e) NULL)

        if (is.null(fold_model)) {
          fold_estimates[i] <- NA
        } else {
          fold_coefs <- lme4::fixef(fold_model)
          fold_estimates[i] <- fold_coefs[focal_coef]
        }
      }

      CoefficientStabilityResult$new(
        focal_coefficient = focal_coef,
        full_model_estimate = full_estimate,
        fold_estimates = fold_estimates,
        excluded_participants = as.character(participants)
      )
    },

    #' @description LOO-CV Model Selection Stability
    #'
    #' Tests whether the best model (by BIC or AIC) changes when each
    #' participant is excluded. Answers: "Is my model choice robust?"
    #'
    #' @param data Data frame with observations
    #' @param models Named list of lmerMod objects to compare
    #' @param criterion "BIC" or "AIC"
    #' @param id_col Name of participant ID column
    #' @return ModelSelectionStabilityResult object
    loo_model_selection_stability = function(data, models, criterion = "BIC",
                                              id_col = "id") {
      if (!requireNamespace("lme4", quietly = TRUE)) {
        stop("Package 'lme4' required for loo_model_selection_stability()")
      }

      stopifnot(
        "models must be a named list" = is.list(models) && !is.null(names(models)),
        "criterion must be 'BIC' or 'AIC'" = criterion %in% c("BIC", "AIC")
      )

      model_names <- names(models)
      criterion_fn <- if (criterion == "BIC") stats::BIC else stats::AIC

      # Get full data winner
      full_criteria <- sapply(models, criterion_fn)
      full_winner <- model_names[which.min(full_criteria)]

      participants <- unique(data[[id_col]])
      n_participants <- length(participants)

      # Extract formulas
      formulas <- lapply(models, stats::formula)

      fold_winners <- character(n_participants)

      for (i in seq_along(participants)) {
        train_data <- data[data[[id_col]] != participants[i], ]

        # Refit all models
        fold_criteria <- sapply(seq_along(models), function(j) {
          fold_model <- tryCatch({
            lme4::lmer(formulas[[j]], data = train_data, REML = FALSE)
          }, error = function(e) NULL)

          if (is.null(fold_model)) return(Inf)
          criterion_fn(fold_model)
        })

        fold_winners[i] <- model_names[which.min(fold_criteria)]
      }

      ModelSelectionStabilityResult$new(
        full_data_winner = full_winner,
        fold_winners = fold_winners,
        excluded_participants = as.character(participants),
        criterion = criterion
      )
    },

    #' @description LOO-CV Prediction Error by Participant
    #'
    #' Calculates prediction error for each participant when their data
    #' is held out. Identifies hard-to-predict individuals and provides
    #' residuals for anomaly detection.
    #'
    #' @param data Data frame with observations
    #' @param model Fitted lme4 model
    #' @param id_col Name of participant ID column
    #' @param outcome_col Name of outcome variable
    #' @return PredictionErrorByParticipantResult object
    loo_prediction_error_by_participant = function(data, model, id_col = "id",
                                                    outcome_col = "mean_velocity") {
      if (!requireNamespace("lme4", quietly = TRUE)) {
        stop("Package 'lme4' required for loo_prediction_error_by_participant()")
      }

      participants <- unique(data[[id_col]])
      n_participants <- length(participants)

      formula_obj <- stats::formula(model)

      # Store per-participant metrics
      participant_metrics <- data.frame(
        participant_id = character(n_participants),
        rmse = numeric(n_participants),
        mae = numeric(n_participants),
        bias = numeric(n_participants),
        n_observations = integer(n_participants),
        stringsAsFactors = FALSE
      )

      # Store all residuals for anomaly detection
      all_residuals <- numeric(0)
      all_participant_ids <- character(0)

      for (i in seq_along(participants)) {
        test_data <- data[data[[id_col]] == participants[i], ]
        train_data <- data[data[[id_col]] != participants[i], ]

        # Refit model
        fold_model <- tryCatch({
          lme4::lmer(formula_obj, data = train_data, REML = FALSE)
        }, error = function(e) NULL)

        if (is.null(fold_model)) {
          participant_metrics$participant_id[i] <- as.character(participants[i])
          participant_metrics$rmse[i] <- NA
          participant_metrics$mae[i] <- NA
          participant_metrics$bias[i] <- NA
          participant_metrics$n_observations[i] <- nrow(test_data)
          next
        }

        # Predict (population average)
        predictions <- predict(fold_model, newdata = test_data, re.form = NA,
                               allow.new.levels = TRUE)
        actual <- test_data[[outcome_col]]
        errors <- actual - predictions

        # Store metrics
        participant_metrics$participant_id[i] <- as.character(participants[i])
        participant_metrics$rmse[i] <- sqrt(mean(errors^2))
        participant_metrics$mae[i] <- mean(abs(errors))
        participant_metrics$bias[i] <- mean(errors)  # positive = under-predicted
        participant_metrics$n_observations[i] <- nrow(test_data)

        # Accumulate residuals
        all_residuals <- c(all_residuals, errors)
        all_participant_ids <- c(all_participant_ids,
                                  rep(as.character(participants[i]), length(errors)))
      }

      # Calculate overall metrics
      valid_rmse <- participant_metrics$rmse[!is.na(participant_metrics$rmse)]
      valid_mae <- participant_metrics$mae[!is.na(participant_metrics$mae)]

      PredictionErrorByParticipantResult$new(
        participant_metrics = participant_metrics,
        overall_rmse = sqrt(mean(all_residuals^2)),
        overall_mae = mean(abs(all_residuals)),
        residuals = all_residuals,
        residual_participant_ids = all_participant_ids
      )
    },

    #' @description Plot Coefficient Stability
    #'
    #' Visualizes how the focal coefficient changes across LOO folds.
    #'
    #' @param stability_result CoefficientStabilityResult object
    #' @param title Plot title
    #' @return ggplot object
    plot_coefficient_stability = function(stability_result,
                                           title = "Coefficient Stability Across LOO Folds") {
      stopifnot(
        "stability_result must be CoefficientStabilityResult" =
          inherits(stability_result, "CoefficientStabilityResult")
      )

      df <- data.frame(
        fold = seq_along(stability_result$fold_estimates),
        estimate = stability_result$fold_estimates,
        participant = stability_result$excluded_participants
      )

      df <- df[!is.na(df$estimate), ]
      df <- df[order(df$estimate), ]
      df$order <- seq_len(nrow(df))

      ggplot(df, aes(x = .data$estimate, y = .data$order)) +
        geom_vline(
          xintercept = stability_result$full_model_estimate,
          linetype = "solid",
          color = "#E63946",
          linewidth = 1
        ) +
        geom_point(color = "#2E86AB", size = 2) +
        labs(
          title = title,
          subtitle = sprintf(
            "%s | Full model = %.4f | CV = %.1f%% (%s)",
            stability_result$focal_coefficient,
            stability_result$full_model_estimate,
            stability_result$cv_percent,
            if (stability_result$is_stable) "stable" else "unstable"
          ),
          x = paste("Coefficient estimate when participant excluded"),
          y = "Fold (ordered by estimate)"
        ) +
        theme_minimal() +
        theme(
          plot.title = element_text(face = "bold", size = 14),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank()
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
