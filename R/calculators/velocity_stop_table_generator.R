# R/calculators/velocity_stop_table_generator.R
# Service: Velocity Stop Table Generation for VBT
#
# =============================================================================
# PRINCIPLES APPLIED (from CLAUDE.md)
# =============================================================================
#
# NAMING PRINCIPLES:
# - Understandability: Domain terms (velocity_threshold, rir_zone, stop_velocity)
# - Consistency: All public methods verb-based (generate, format, export)
# - Distinguishability: generate_table vs generate_zones
# - Conciseness: vbt for Velocity-Based Training
#
# SOLID PRINCIPLES:
# - SRP: Single responsibility - generate velocity thresholds for training
# - OCP: Open for extension with different table formats and export types
# - DIP: Depends on LmmAnalyzer abstraction, not specific model implementations
#
# CUPID PRINCIPLES:
# - Composable: Uses LmmAnalyzer output to generate practitioner tables
# - Unix: Each method does one thing (generate, format, export)
# - Predictable: Same model -> same velocity thresholds
# - Idiomatic: Follows R conventions for data frame output
# - Domain-based: Names reflect VBT practitioner concepts (zones, thresholds)
#
# SCIENTIFIC VALIDITY:
# - Velocity thresholds derived from population LMM fixed effects
# - Prediction intervals from model standard error
# - RIR zones based on literature (0-1: failure, 2-3: effective, 4+: submaximal)
#
# PRACTICAL APPLICATION:
# - Generates autoregulation tables for strength coaches
# - Provides velocity targets for each RIR level
# - Includes prediction intervals for practical guidance
# =============================================================================

box::use(
  R6[R6Class],
  ./lmm_analyzer[LmmAnalyzer, LmmModelResult]
)

#' Velocity Stop Table Result
#'
#' Value object containing velocity thresholds for each RIR level.
#'
#' @export
VelocityStopTable <- R6Class(

  classname = "VelocityStopTable",
  cloneable = FALSE,

  public = list(

    #' @description Create velocity stop table
    #' @param table_type "global" or "individual"
    #' @param id Participant ID (NULL for global)
    #' @param table Data frame with RIR, velocity, and prediction error
    #' @param model_used The LmmModelResult used to generate the table
    #' @param load_included Whether load percentage was included as predictor
    initialize = function(table_type, id, table, model_used, load_included) {
      private$.table_type <- table_type
      private$.id <- id
      private$.table <- table
      private$.model_used <- model_used
      private$.load_included <- load_included
    },

    #' @description Convert to list
    to_list = function() {
      list(
        table_type = private$.table_type,
        id = private$.id,
        table = private$.table,
        load_included = private$.load_included
      )
    },

    #' @description Format table for display
    format_for_display = function() {
      tbl <- private$.table
      # Round velocities to 3 decimals
      numeric_cols <- sapply(tbl, is.numeric)
      tbl[numeric_cols] <- lapply(tbl[numeric_cols], round, 3)
      tbl
    }
  ),

  active = list(
    table_type = function() private$.table_type,
    id = function() private$.id,
    table = function() private$.table,
    model_used = function() private$.model_used,
    load_included = function() private$.load_included
  ),

  private = list(
    .table_type = NULL,
    .id = NULL,
    .table = NULL,
    .model_used = NULL,
    .load_included = NULL
  )
)


#' Load Importance Test Result
#'
#' Value object containing results of testing whether load percentage matters.
#'
#' @export
LoadImportanceResult <- R6Class(

  classname = "LoadImportanceResult",
  cloneable = FALSE,

  public = list(

    #' @description Create load importance result
    #' @param load_significant Whether load significantly improves model
    #' @param lrt_p_value P-value from likelihood ratio test
    #' @param delta_aic AIC difference (positive = load model better)
    #' @param delta_bic BIC difference
    #' @param bayes_factor BF approximation
    #' @param bf_interpretation Interpretation of Bayes factor
    #' @param recommendation "global" or "load_specific"
    initialize = function(load_significant,
                          lrt_p_value,
                          delta_aic,
                          delta_bic,
                          bayes_factor,
                          bf_interpretation,
                          recommendation) {
      private$.load_significant <- load_significant
      private$.lrt_p_value <- lrt_p_value
      private$.delta_aic <- delta_aic
      private$.delta_bic <- delta_bic
      private$.bayes_factor <- bayes_factor
      private$.bf_interpretation <- bf_interpretation
      private$.recommendation <- recommendation
    },

    #' @description Convert to list
    to_list = function() {
      list(
        load_significant = private$.load_significant,
        lrt_p_value = private$.lrt_p_value,
        delta_aic = private$.delta_aic,
        delta_bic = private$.delta_bic,
        bayes_factor = private$.bayes_factor,
        bf_interpretation = private$.bf_interpretation,
        recommendation = private$.recommendation
      )
    }
  ),

  active = list(
    load_significant = function() private$.load_significant,
    lrt_p_value = function() private$.lrt_p_value,
    delta_aic = function() private$.delta_aic,
    delta_bic = function() private$.delta_bic,
    bayes_factor = function() private$.bayes_factor,
    bf_interpretation = function() private$.bf_interpretation,
    recommendation = function() private$.recommendation
  ),

  private = list(
    .load_significant = NULL,
    .lrt_p_value = NULL,
    .delta_aic = NULL,
    .delta_bic = NULL,
    .bayes_factor = NULL,
    .bf_interpretation = NULL,
    .recommendation = NULL
  )
)


#' Velocity Stop Table Generator
#'
#' Generates velocity stop tables for velocity-based training.
#' Tests whether load percentage matters, then creates appropriate tables.
#'
#' @export
VelocityStopTableGenerator <- R6Class(

  classname = "VelocityStopTableGenerator",

  public = list(

    #' @description Create generator
    initialize = function() {
      private$.lmm_analyzer <- LmmAnalyzer$new()
    },

    #' @description Test if load percentage significantly affects velocity-RIR relationship
    #' @param data Data frame with rir, mean_velocity, load_percentage, and id columns
    #' @return LoadImportanceResult object
    test_load_importance = function(data) {
      # Model without load
      base_formula <- mean_velocity ~ rir

      # Model with load
      full_formula <- mean_velocity ~ rir + load_percentage

      # Test importance
      test_result <- private$.lmm_analyzer$test_variable_importance(
        data = data,
        base_formula = base_formula,
        full_formula = full_formula,
        random_formula = ~1 + rir | id
      )

      # Determine recommendation based on multiple criteria
      # Prefer simpler model (global) unless strong evidence for load
      load_significant <- test_result$variable_significant
      bf_favors_complex <- test_result$bayes_factor < 1 / 3
      aic_favors_complex <- test_result$delta_aic > 2

      # Recommendation: use load-specific only if LRT significant AND (BF or AIC support)
      recommendation <- if (load_significant && (bf_favors_complex || aic_favors_complex)) {
        "load_specific"
      } else {
        "global"
      }

      LoadImportanceResult$new(
        load_significant = load_significant,
        lrt_p_value = test_result$lrt_p_value,
        delta_aic = test_result$delta_aic,
        delta_bic = test_result$delta_bic,
        bayes_factor = test_result$bayes_factor,
        bf_interpretation = test_result$bf_interpretation,
        recommendation = recommendation
      )
    },

    #' @description Generate global velocity stop table (load-independent)
    #' @param data Data frame with observations
    #' @param rir_targets Vector of RIR values to generate thresholds for
    #' @param include_random Include random effects in prediction
    #' @return VelocityStopTable object
    generate_global_table = function(data, rir_targets = 0:5, include_random = FALSE) {
      # Fit model without load percentage
      model_result <- private$.lmm_analyzer$fit(
        data = data,
        formula = mean_velocity ~ rir,
        random_formula = ~1 + rir | id,
        model_name = "global_velocity_model"
      )

      # Generate velocity predictions for each RIR
      table <- private$.generate_velocity_table(
        model_result = model_result,
        data = data,
        rir_targets = rir_targets,
        include_load = FALSE,
        include_random = include_random
      )

      VelocityStopTable$new(
        table_type = "global",
        id = NULL,
        table = table,
        model_used = model_result,
        load_included = FALSE
      )
    },

    #' @description Generate load-specific velocity stop table
    #' @param data Data frame with observations
    #' @param rir_targets Vector of RIR values
    #' @param include_random Include random effects
    #' @return VelocityStopTable object
    generate_load_specific_table = function(data, rir_targets = 0:5, include_random = FALSE) {
      # Fit model with load percentage
      model_result <- private$.lmm_analyzer$fit(
        data = data,
        formula = mean_velocity ~ rir + load_percentage,
        random_formula = ~1 + rir | id,
        model_name = "load_specific_velocity_model"
      )

      # Generate velocity predictions for each RIR and load combination
      table <- private$.generate_velocity_table(
        model_result = model_result,
        data = data,
        rir_targets = rir_targets,
        include_load = TRUE,
        include_random = include_random
      )

      VelocityStopTable$new(
        table_type = "load_specific",
        id = NULL,
        table = table,
        model_used = model_result,
        load_included = TRUE
      )
    },

    #' @description Generate individual velocity stop table for a specific participant
    #' @param data Data frame with observations
    #' @param id Participant ID
    #' @param rir_targets Vector of RIR values
    #' @return VelocityStopTable object
    generate_individual_table = function(data, id, rir_targets = 0:5) {
      # Filter to participant data
      participant_data <- data[data$id == id, ]

      if (nrow(participant_data) < 5) {
        warning(paste("Insufficient data for participant", id))
        return(NULL)
      }

      # Fit simple linear model for individual
      model <- stats::lm(mean_velocity ~ rir, data = participant_data)

      # Generate predictions
      newdata <- data.frame(rir = rir_targets)
      predictions <- stats::predict(model, newdata = newdata, se.fit = TRUE)

      table <- data.frame(
        rir = rir_targets,
        velocity = predictions$fit,
        se = predictions$se.fit,
        lower_95 = predictions$fit - 1.96 * predictions$se.fit,
        upper_95 = predictions$fit + 1.96 * predictions$se.fit
      )

      VelocityStopTable$new(
        table_type = "individual",
        id = id,
        table = table,
        model_used = NULL,
        load_included = FALSE
      )
    },

    #' @description Generate individual tables for all participants
    #' @param data Data frame with observations
    #' @param rir_targets Vector of RIR values
    #' @return Named list of VelocityStopTable objects
    generate_all_individual_tables = function(data, rir_targets = 0:5) {
      ids <- unique(data$id)
      tables <- list()

      for (id in ids) {
        table <- self$generate_individual_table(data, id, rir_targets)
        if (!is.null(table)) {
          tables[[as.character(id)]] <- table
        }
      }

      tables
    },

    #' @description Compare general vs individual table accuracy
    #' @param data Data frame with observations
    #' @param test_data Optional separate test data (default: use same data)
    #' @return List with comparison metrics
    compare_general_vs_individual = function(data, test_data = NULL) {
      if (is.null(test_data)) {
        test_data <- data
      }

      # Generate global table
      global_table <- self$generate_global_table(data)
      global_model <- global_table$model_used

      # Predict using global model
      global_predictions <- private$.lmm_analyzer$predict_values(
        global_model,
        test_data,
        include_random = FALSE
      )

      global_mae <- mean(abs(test_data$mean_velocity - global_predictions), na.rm = TRUE)
      global_rmse <- sqrt(mean((test_data$mean_velocity - global_predictions)^2, na.rm = TRUE))

      # Predict using individual models (with random effects)
      individual_predictions <- private$.lmm_analyzer$predict_values(
        global_model,
        test_data,
        include_random = TRUE
      )

      individual_mae <- mean(abs(test_data$mean_velocity - individual_predictions), na.rm = TRUE)
      individual_rmse <- sqrt(mean((test_data$mean_velocity - individual_predictions)^2, na.rm = TRUE))

      # Calculate improvement
      mae_improvement <- (global_mae - individual_mae) / global_mae * 100
      rmse_improvement <- (global_rmse - individual_rmse) / global_rmse * 100

      list(
        global_mae = global_mae,
        global_rmse = global_rmse,
        individual_mae = individual_mae,
        individual_rmse = individual_rmse,
        mae_improvement_pct = mae_improvement,
        rmse_improvement_pct = rmse_improvement,
        recommendation = if (mae_improvement > 10) {
          "Individual calibration recommended (>10% improvement)"
        } else {
          "Global table sufficient (<10% improvement)"
        }
      )
    },

    #' @description Export table to CSV
    #' @param velocity_table VelocityStopTable object
    #' @param path File path for CSV output
    export_to_csv = function(velocity_table, path) {
      utils::write.csv(velocity_table$table, path, row.names = FALSE)
    }
  ),

  private = list(

    .lmm_analyzer = NULL,

    .generate_velocity_table = function(model_result, data, rir_targets, include_load, include_random) {
      if (include_load) {
        # Generate for each load level
        loads <- unique(data$load_percentage)
        tables <- list()

        for (load in loads) {
          newdata <- data.frame(
            rir = rir_targets,
            load_percentage = load,
            id = data$id[1]  # Dummy ID for prediction
          )

          predictions <- private$.lmm_analyzer$predict_values(
            model_result,
            newdata,
            include_random = include_random
          )

          tables[[as.character(load)]] <- data.frame(
            rir = rir_targets,
            load_percentage = load,
            velocity = predictions
          )
        }

        # Combine tables
        table <- do.call(rbind, tables)

        # Add prediction errors (SE from residuals)
        residual_se <- sqrt(mean(stats::residuals(model_result$model)^2))
        table$se <- residual_se
        table$lower_95 <- table$velocity - 1.96 * residual_se
        table$upper_95 <- table$velocity + 1.96 * residual_se

      } else {
        # Global table (no load distinction)
        newdata <- data.frame(
          rir = rir_targets,
          id = data$id[1]  # Dummy ID for prediction
        )

        predictions <- private$.lmm_analyzer$predict_values(
          model_result,
          newdata,
          include_random = include_random
        )

        # Add prediction errors
        residual_se <- sqrt(mean(stats::residuals(model_result$model)^2))

        table <- data.frame(
          rir = rir_targets,
          velocity = predictions,
          se = residual_se,
          lower_95 = predictions - 1.96 * residual_se,
          upper_95 = predictions + 1.96 * residual_se
        )
      }

      table
    }
  )
)
