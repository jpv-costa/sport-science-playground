# R/calculators/marginal_effects_analyzer.R
# Marginal Effects Analysis for Mixed Effects Models
#
# =============================================================================
# EDUCATIONAL OVERVIEW
# =============================================================================
#
# What are Marginal Effects?
# --------------------------
# Marginal effects answer: "What is the expected change in the outcome
# for a unit change in the predictor, averaging across all individuals?"
#
# For mixed models, there are two perspectives:
# 1. **Population average**: Effect for a "typical" individual (re.form = NA)
# 2. **Subject-specific**: Effect accounting for each individual's random effects
#
# Why Marginal Effects Matter:
# - Provide interpretable summaries of model predictions
# - Show the practical magnitude of effects
# - Visualize how effects vary across conditions
#
# =============================================================================

box::use(
  R6[R6Class],
  stats[predict, coef, vcov, qnorm, pnorm, confint, lm, sigma, quantile],
  ggplot2[ggplot, aes, geom_line, geom_ribbon, geom_point, geom_errorbar,
          labs, theme_minimal, theme, element_text, scale_color_manual,
          scale_fill_manual, facet_wrap]
)

#' Marginal Effect Result
#'
#' Immutable value object containing marginal effect calculations
#'
#' @export
MarginalEffectResult <- R6Class(
  classname = "MarginalEffectResult",
  cloneable = FALSE,

  public = list(
    #' @field focal_variable Name of the focal predictor
    focal_variable = NULL,
    #' @field effect_estimate Point estimate of the marginal effect
    effect_estimate = NULL,
    #' @field se Standard error of the effect
    se = NULL,
    #' @field ci_lower Lower bound of confidence interval
    ci_lower = NULL,
    #' @field ci_upper Upper bound of confidence interval
    ci_upper = NULL,
    #' @field prediction_data Data frame with predictions across focal variable range
    prediction_data = NULL,

    #' @description Create a new marginal effect result
    initialize = function(focal_variable, effect_estimate, se, ci_lower, ci_upper,
                          prediction_data) {
      self$focal_variable <- focal_variable
      self$effect_estimate <- effect_estimate
      self$se <- se
      self$ci_lower <- ci_lower
      self$ci_upper <- ci_upper
      self$prediction_data <- prediction_data
    },

    #' @description Interpret effect magnitude
    interpret = function(outcome_unit = "m/s") {
      paste0(
        "For each 1-unit increase in ", self$focal_variable, ", ",
        "velocity changes by ", round(self$effect_estimate, 4), " ", outcome_unit,
        " (95% CI: ", round(self$ci_lower, 4), " to ", round(self$ci_upper, 4), ")"
      )
    },

    #' @description Convert to list for serialization
    to_list = function() {
      list(
        focal_variable = self$focal_variable,
        effect_estimate = self$effect_estimate,
        se = self$se,
        ci_lower = self$ci_lower,
        ci_upper = self$ci_upper
      )
    }
  )
)

#' Conditional Effect Result
#'
#' Marginal effects at specific values of a conditioning variable
#'
#' @export
ConditionalEffectResult <- R6Class(
  classname = "ConditionalEffectResult",
  cloneable = FALSE,

  public = list(
    #' @field focal_variable Name of the focal predictor
    focal_variable = NULL,
    #' @field conditioning_variable Name of the conditioning variable
    conditioning_variable = NULL,
    #' @field effects_at_levels Data frame with effects at different levels
    effects_at_levels = NULL,

    #' @description Create a new conditional effect result
    initialize = function(focal_variable, conditioning_variable, effects_at_levels) {
      self$focal_variable <- focal_variable
      self$conditioning_variable <- conditioning_variable
      self$effects_at_levels <- effects_at_levels
    },

    #' @description Convert to list for serialization
    to_list = function() {
      list(
        focal_variable = self$focal_variable,
        conditioning_variable = self$conditioning_variable,
        effects_at_levels = self$effects_at_levels
      )
    }
  )
)

#' Marginal Effects Analyzer
#'
#' R6 class for calculating and visualizing marginal effects
#' from mixed effects models.
#'
#' @section Key Concepts:
#' - **Average Marginal Effect (AME)**: The effect averaged over all
#'   observed values of other predictors
#' - **Marginal Effect at the Mean (MEM)**: The effect when other
#'   predictors are at their means
#' - **Conditional Marginal Effect**: The effect at specific values
#'   of other predictors
#'
#' @export
MarginalEffectsAnalyzer <- R6Class(
  classname = "MarginalEffectsAnalyzer",
  cloneable = FALSE,

  public = list(
    #' @description Create a new MarginalEffectsAnalyzer instance
    initialize = function() {
      # No initialization needed
    },

    # =========================================================================
    # MARGINAL EFFECTS
    # =========================================================================

    #' @description Calculate Marginal Effect
    #'
    #' Computes the marginal effect of a focal predictor, holding other
    #' variables at specified values or their means.
    #'
    #' @param model Fitted lme4 model
    #' @param focal_var Name of the focal predictor
    #' @param data Data frame used for fitting
    #' @param at Named list of values to hold other variables at (default: means)
    #' @param n_points Number of points for prediction curve
    #' @param level Confidence level
    #' @return MarginalEffectResult object
    calculate_marginal_effect = function(model, focal_var, data, at = NULL,
                                         n_points = 50, level = 0.95) {
      if (!requireNamespace("lme4", quietly = TRUE)) {
        stop("Package 'lme4' required for calculate_marginal_effect()")
      }

      # Get fixed effects
      fixed_effects <- lme4::fixef(model)

      # Extract the coefficient for focal variable
      if (!focal_var %in% names(fixed_effects)) {
        stop("Focal variable '", focal_var, "' not found in model fixed effects")
      }

      effect_estimate <- fixed_effects[focal_var]

      # Get standard error from variance-covariance matrix
      vcov_mat <- as.matrix(vcov(model))
      se <- sqrt(vcov_mat[focal_var, focal_var])

      # Calculate CI
      z_value <- qnorm((1 + level) / 2)
      ci_lower <- effect_estimate - z_value * se
      ci_upper <- effect_estimate + z_value * se

      # Create prediction data across range of focal variable
      focal_range <- range(data[[focal_var]], na.rm = TRUE)
      focal_seq <- seq(focal_range[1], focal_range[2], length.out = n_points)

      # Create new data for predictions
      pred_data <- private$.create_prediction_data(data, focal_var, focal_seq, at)

      # Get predictions (population average)
      pred_data$predicted <- predict(model, newdata = pred_data, re.form = NA,
                                     allow.new.levels = TRUE)

      # Add approximate CIs for predictions
      # This is simplified - full treatment requires delta method or bootstrap
      pred_data$pred_se <- se * abs(pred_data[[focal_var]] - mean(data[[focal_var]], na.rm = TRUE))
      pred_data$pred_lower <- pred_data$predicted - z_value * sqrt(pred_data$pred_se^2 + sigma(model)^2)
      pred_data$pred_upper <- pred_data$predicted + z_value * sqrt(pred_data$pred_se^2 + sigma(model)^2)

      MarginalEffectResult$new(
        focal_variable = focal_var,
        effect_estimate = effect_estimate,
        se = se,
        ci_lower = ci_lower,
        ci_upper = ci_upper,
        prediction_data = pred_data
      )
    },

    # =========================================================================
    # CONDITIONAL EFFECTS
    # =========================================================================

    #' @description Calculate Conditional Effects
    #'
    #' Computes how the effect of the focal predictor varies across
    #' levels of a conditioning variable.
    #'
    #' @param model Fitted lme4 model
    #' @param focal_var Name of the focal predictor
    #' @param by_var Name of the conditioning variable
    #' @param data Data frame used for fitting
    #' @param by_values Specific values of by_var to use (default: quartiles)
    #' @param n_points Number of points for each prediction curve
    #' @param level Confidence level
    #' @return ConditionalEffectResult object
    calculate_conditional_effects = function(model, focal_var, by_var, data,
                                             by_values = NULL, n_points = 50,
                                             level = 0.95) {
      if (!requireNamespace("lme4", quietly = TRUE)) {
        stop("Package 'lme4' required for calculate_conditional_effects()")
      }

      # Determine values of conditioning variable
      if (is.null(by_values)) {
        # Use quartiles by default
        by_values <- quantile(data[[by_var]], probs = c(0.25, 0.5, 0.75), na.rm = TRUE)
      }

      # Range of focal variable
      focal_range <- range(data[[focal_var]], na.rm = TRUE)
      focal_seq <- seq(focal_range[1], focal_range[2], length.out = n_points)

      # Create predictions at each level of by_var
      all_predictions <- lapply(by_values, function(by_val) {
        pred_data <- private$.create_prediction_data(data, focal_var, focal_seq,
                                                     at = setNames(list(by_val), by_var))
        pred_data$predicted <- predict(model, newdata = pred_data, re.form = NA,
                                       allow.new.levels = TRUE)
        pred_data$by_level <- paste0(by_var, " = ", round(by_val, 2))
        pred_data$by_value <- by_val
        pred_data
      })

      combined_predictions <- do.call(rbind, all_predictions)

      # Calculate effect at each level
      effects_at_levels <- data.frame(
        by_value = by_values,
        by_label = paste0(by_var, " = ", round(by_values, 2)),
        stringsAsFactors = FALSE
      )

      # Get fixed effects (effect should be constant if no interaction)
      fixed_effects <- lme4::fixef(model)
      effects_at_levels$effect <- fixed_effects[focal_var]

      # Check for interaction
      interaction_term <- paste0(focal_var, ":", by_var)
      alt_interaction <- paste0(by_var, ":", focal_var)

      if (interaction_term %in% names(fixed_effects)) {
        interaction_coef <- fixed_effects[interaction_term]
        effects_at_levels$effect <- fixed_effects[focal_var] + interaction_coef * by_values
      } else if (alt_interaction %in% names(fixed_effects)) {
        interaction_coef <- fixed_effects[alt_interaction]
        effects_at_levels$effect <- fixed_effects[focal_var] + interaction_coef * by_values
      }

      ConditionalEffectResult$new(
        focal_variable = focal_var,
        conditioning_variable = by_var,
        effects_at_levels = list(
          effect_table = effects_at_levels,
          prediction_curves = combined_predictions
        )
      )
    },

    # =========================================================================
    # VISUALIZATIONS
    # =========================================================================

    #' @description Plot Marginal Effect
    #'
    #' Creates a plot showing predicted values across the range of
    #' the focal predictor with confidence band.
    #'
    #' @param effect_result MarginalEffectResult object
    #' @param data Optional data frame to add observed points
    #' @param outcome_col Name of outcome column in data
    #' @param title Plot title
    #' @return ggplot object
    plot_marginal_effect = function(effect_result, data = NULL, outcome_col = NULL,
                                    title = NULL) {
      pred_data <- effect_result$prediction_data
      focal_var <- effect_result$focal_variable

      if (is.null(title)) {
        title <- paste0("Marginal Effect of ", focal_var)
      }

      p <- ggplot(pred_data, aes_string(x = focal_var)) +
        geom_ribbon(aes(ymin = pred_lower, ymax = pred_upper),
                    fill = "#2E86AB", alpha = 0.2) +
        geom_line(aes(y = predicted), color = "#2E86AB", linewidth = 1.2) +
        labs(
          title = title,
          subtitle = paste0(
            "Effect = ", round(effect_result$effect_estimate, 4),
            " (95% CI: ", round(effect_result$ci_lower, 4),
            " to ", round(effect_result$ci_upper, 4), ")"
          ),
          x = focal_var,
          y = "Predicted Velocity (m/s)",
          caption = "Shaded band = 95% prediction interval"
        ) +
        theme_minimal() +
        theme(
          plot.title = element_text(face = "bold", size = 14),
          plot.subtitle = element_text(size = 10, color = "#666666")
        )

      # Add observed data points if provided
      if (!is.null(data) && !is.null(outcome_col)) {
        p <- p + geom_point(data = data, aes_string(x = focal_var, y = outcome_col),
                            alpha = 0.3, color = "#1D3557")
      }

      p
    },

    #' @description Plot Conditional Effects
    #'
    #' Shows how the effect varies across levels of a conditioning variable.
    #'
    #' @param cond_result ConditionalEffectResult object
    #' @param title Plot title
    #' @return ggplot object
    plot_conditional_effects = function(cond_result, title = NULL) {
      pred_curves <- cond_result$effects_at_levels$prediction_curves
      focal_var <- cond_result$focal_variable
      by_var <- cond_result$conditioning_variable

      if (is.null(title)) {
        title <- paste0("Effect of ", focal_var, " by ", by_var)
      }

      ggplot(pred_curves, aes_string(x = focal_var, y = "predicted",
                                      color = "by_level", fill = "by_level")) +
        geom_line(linewidth = 1.2) +
        labs(
          title = title,
          subtitle = paste0("Predicted velocity across ", focal_var,
                            " at different levels of ", by_var),
          x = focal_var,
          y = "Predicted Velocity (m/s)",
          color = by_var,
          fill = by_var
        ) +
        theme_minimal() +
        theme(
          plot.title = element_text(face = "bold", size = 14),
          plot.subtitle = element_text(size = 10, color = "#666666"),
          legend.position = "bottom"
        )
    },

    #' @description Create Effect Summary Table
    #'
    #' Returns a formatted table of marginal effects for multiple predictors.
    #'
    #' @param model Fitted lme4 model
    #' @param data Data frame used for fitting
    #' @param predictors Vector of predictor names (default: all fixed effects)
    #' @param level Confidence level
    #' @return Data frame with effect summaries
    create_effect_table = function(model, data, predictors = NULL, level = 0.95) {
      if (!requireNamespace("lme4", quietly = TRUE)) {
        stop("Package 'lme4' required for create_effect_table()")
      }

      fixed_effects <- lme4::fixef(model)
      vcov_mat <- as.matrix(vcov(model))

      if (is.null(predictors)) {
        predictors <- names(fixed_effects)[-1]  # Exclude intercept
      }

      z_value <- qnorm((1 + level) / 2)

      effect_table <- data.frame(
        predictor = predictors,
        estimate = fixed_effects[predictors],
        se = sqrt(diag(vcov_mat)[predictors]),
        stringsAsFactors = FALSE
      )

      effect_table$ci_lower <- effect_table$estimate - z_value * effect_table$se
      effect_table$ci_upper <- effect_table$estimate + z_value * effect_table$se
      effect_table$z_value <- effect_table$estimate / effect_table$se
      effect_table$p_value <- 2 * (1 - pnorm(abs(effect_table$z_value)))

      # Round for display
      effect_table$estimate <- round(effect_table$estimate, 4)
      effect_table$se <- round(effect_table$se, 4)
      effect_table$ci_lower <- round(effect_table$ci_lower, 4)
      effect_table$ci_upper <- round(effect_table$ci_upper, 4)
      effect_table$p_value <- round(effect_table$p_value, 4)

      rownames(effect_table) <- NULL
      effect_table
    }
  ),

  private = list(
    #' Create prediction data frame
    .create_prediction_data = function(data, focal_var, focal_seq, at = NULL) {
      # Get numeric columns for means
      numeric_cols <- sapply(data, is.numeric)
      col_means <- colMeans(data[, numeric_cols, drop = FALSE], na.rm = TRUE)

      # Create base prediction data
      n_points <- length(focal_seq)
      pred_data <- data.frame(matrix(NA, nrow = n_points, ncol = ncol(data)))
      names(pred_data) <- names(data)

      # Set focal variable
      pred_data[[focal_var]] <- focal_seq

      # Set other numeric variables to means
      for (col in names(col_means)) {
        if (col != focal_var && col %in% names(pred_data)) {
          pred_data[[col]] <- col_means[col]
        }
      }

      # Override with specified values
      if (!is.null(at)) {
        for (var_name in names(at)) {
          if (var_name %in% names(pred_data)) {
            pred_data[[var_name]] <- at[[var_name]]
          }
        }
      }

      # Handle factors - use first level
      factor_cols <- sapply(data, is.factor)
      for (col in names(data)[factor_cols]) {
        if (col %in% names(pred_data)) {
          pred_data[[col]] <- levels(data[[col]])[1]
        }
      }

      pred_data
    }
  )
)

# =============================================================================
# CONVENIENCE FUNCTIONS
# =============================================================================

#' Calculate Marginal Effect
#'
#' @param model Fitted lme4 model
#' @param focal_var Name of focal predictor
#' @param data Data frame
#' @param ... Additional arguments
#' @return MarginalEffectResult object
#' @export
calculate_marginal_effect <- function(model, focal_var, data, ...) {
  analyzer <- MarginalEffectsAnalyzer$new()
  analyzer$calculate_marginal_effect(model, focal_var, data, ...)
}

#' Plot Marginal Effect
#'
#' @param effect_result MarginalEffectResult object
#' @param ... Additional arguments
#' @return ggplot object
#' @export
plot_marginal_effect <- function(effect_result, ...) {
  analyzer <- MarginalEffectsAnalyzer$new()
  analyzer$plot_marginal_effect(effect_result, ...)
}
