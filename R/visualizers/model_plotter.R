# R/visualizers/model_plotter.R
# Service: Model Visualization for Mixed Effects Models
#
# =============================================================================
# PRINCIPLES APPLIED (from CLAUDE.md)
# =============================================================================
#
# NAMING PRINCIPLES:
# - Understandability: Problem domain terms (model_plot, caterpillar, prediction)
# - Consistency: All public methods start with "plot_", extractors with "extract_"
# - Distinguishability: Specific names (plot_model vs plot_predictions)
# - Conciseness: Short but meaningful (BLUP = Best Linear Unbiased Predictor)
#
# SOLID PRINCIPLES:
# - SRP: Single responsibility - only creates model visualizations
# - OCP: Open for extension via new plot methods
# - DIP: Depends on abstractions (data frames), not lme4 internals
#
# CUPID PRINCIPLES:
# - Composable: Each plot method is independent
# - Unix: Each method does one thing well
# - Predictable: Same inputs -> same outputs
# - Idiomatic: Follows ggplot2/R conventions
# - Domain-based: Names reflect statistical concepts
#
# SCIENTIFIC VALIDITY:
# - Random effects extraction uses lme4::ranef with condVar
# - CIs use 1.96 * SE for 95% intervals
# - Model predictions separate fixed and random components
# =============================================================================

box::use(
  R6[R6Class],
  ggplot2[
    ggplot, aes, geom_point, geom_line, geom_histogram, geom_density,
    geom_segment, geom_vline, geom_hline,
    labs, theme, theme_minimal, element_text, element_blank,
    scale_fill_manual, scale_color_manual,
    after_stat, coord_flip
  ],
  dplyr[mutate, arrange, .data],
  stats[predict],
  patchwork[wrap_plots, plot_annotation]
)

#' Color Palette for Model Visualizations
#'
#' Consistent colorblind-friendly palette.
#'
#' @export
MODEL_COLORS <- list(
  primary = "#2E86AB",
  secondary = "#E63946",
  tertiary = "#F77F00",
  individual = "gray60",
  population = "#E63946",
  gray = "#999999",
  dark_gray = "#333333"
)

#' Model Plotter
#'
#' R6 class for creating mixed model visualizations.
#' Separates data extraction from plotting for testability.
#'
#' @section Principles:
#' - SRP: Only renders model visualizations
#' - Testable: Extraction methods return data, plot methods render
#' - Composable: Individual plots can be combined
#'
#' @export
ModelPlotter <- R6Class(

  classname = "ModelPlotter",
  cloneable = FALSE,

  public = list(
    #' @description Create a new ModelPlotter instance
    #' @param base_size Base font size for plots (default: 14)
    #' @param colors Named list of colors (optional, uses defaults)
    initialize = function(base_size = 14, colors = NULL) {
      private$.base_size <- base_size
      private$.colors <- colors %||% MODEL_COLORS
    },

    #' @description Extract caterpillar plot data from lmer model
    #' @param model An lme4 lmer model object
    #' @param effect_type "intercept" or "slope" (default: both)
    #' @return Data frame with id, effect_type, estimate, lower, upper
    extract_caterpillar_data = function(model, effect_type = c("intercept", "slope")) {
      stopifnot(
        "model must be an lmerMod object" = inherits(model, "lmerMod")
      )

      re <- lme4::ranef(model, condVar = TRUE)
      re_df <- as.data.frame(re[[1]])
      re_df$id <- rownames(re_df)

      # Get conditional variances for CIs
      pv <- attr(re[[1]], "postVar")

      results <- list()

      if ("intercept" %in% effect_type && ncol(re_df) >= 2) {
        intercept_se <- if (!is.null(pv)) sqrt(pv[1, 1, ]) else NA
        results$intercept <- data.frame(
          id = re_df$id,
          effect_type = "intercept",
          estimate = re_df[[1]],
          se = intercept_se,
          lower = re_df[[1]] - 1.96 * intercept_se,
          upper = re_df[[1]] + 1.96 * intercept_se,
          stringsAsFactors = FALSE
        )
      }

      if ("slope" %in% effect_type && ncol(re_df) >= 3) {
        slope_se <- if (!is.null(pv) && dim(pv)[1] >= 2) sqrt(pv[2, 2, ]) else NA
        results$slope <- data.frame(
          id = re_df$id,
          effect_type = "slope",
          estimate = re_df[[2]],
          se = slope_se,
          lower = re_df[[2]] - 1.96 * slope_se,
          upper = re_df[[2]] + 1.96 * slope_se,
          stringsAsFactors = FALSE
        )
      }

      do.call(rbind, results)
    },

    #' @description Create dual caterpillar plot (intercepts and slopes)
    #' @param model An lme4 lmer model object
    #' @param title Plot title
    #' @return ggplot object
    plot_caterpillar_dual = function(model, title = "Individual Random Effects") {
      cat_data <- self$extract_caterpillar_data(model, c("intercept", "slope"))

      # Separate and order by estimate
      intercept_data <- cat_data[cat_data$effect_type == "intercept", ]
      slope_data <- cat_data[cat_data$effect_type == "slope", ]

      intercept_data <- intercept_data[order(intercept_data$estimate), ]
      intercept_data$order <- seq_len(nrow(intercept_data))

      slope_data <- slope_data[order(slope_data$estimate), ]
      slope_data$order <- seq_len(nrow(slope_data))

      p1 <- private$.create_caterpillar_panel(
        intercept_data,
        "Individual Baseline Velocity",
        "Random Intercept (m/s)"
      )

      p2 <- private$.create_caterpillar_panel(
        slope_data,
        "Individual RIR Sensitivity",
        "Random Slope (m/s per RIR)"
      )

      # Combine with patchwork
      wrap_plots(p1, p2, ncol = 2) +
        plot_annotation(
          title = title,
          theme = theme(plot.title = element_text(face = "bold", size = 16))
        )
    },

    #' @description Create model plot with raw data + fitted lines
    #' @param data Data frame with id, rir, mean_velocity
    #' @param model An lme4 lmer model object
    #' @param title Plot title
    #' @return ggplot object
    plot_model = function(data, model, title = "Model Plot: Raw Data + Fitted Lines") {
      stopifnot(
        "data must have 'id' column" = "id" %in% names(data),
        "data must have 'rir' column" = "rir" %in% names(data),
        "data must have 'mean_velocity' column" = "mean_velocity" %in% names(data),
        "model must be an lmerMod object" = inherits(model, "lmerMod")
      )

      # Get predictions
      data$fitted_pop <- predict(model, re.form = NA)
      data$fitted_ind <- predict(model)

      ggplot(data, aes(x = .data$rir, y = .data$mean_velocity)) +
        # Raw data points (colored by individual but no legend)
        geom_point(aes(color = .data$id), alpha = 0.4, size = 2, show.legend = FALSE) +
        # Individual fitted lines
        geom_line(
          aes(y = .data$fitted_ind, group = .data$id),
          color = private$.colors$individual,
          alpha = 0.6,
          linewidth = 0.5
        ) +
        # Population model
        geom_line(
          aes(y = .data$fitted_pop),
          color = private$.colors$population,
          linewidth = 1.5
        ) +
        labs(
          title = title,
          subtitle = "Points = raw data | Gray lines = individual predictions | Red line = population average",
          x = "Repetitions in Reserve (RIR)",
          y = "Mean Velocity (m/s)"
        ) +
        private$.get_base_theme()
    },

    #' @description Create prediction error comparison plot
    #' @param data Data frame with mean_velocity column
    #' @param model An lme4 lmer model object
    #' @param title Plot title
    #' @return ggplot object
    plot_prediction_comparison = function(data, model,
                                           title = "Prediction Error Comparison") {
      stopifnot(
        "data must have 'mean_velocity' column" = "mean_velocity" %in% names(data),
        "model must be an lmerMod object" = inherits(model, "lmerMod")
      )

      # Calculate prediction errors
      pop_pred <- predict(model, re.form = NA)
      ind_pred <- predict(model)

      pop_error <- (data$mean_velocity - pop_pred) * 1000  # Convert to mm/s
      ind_error <- (data$mean_velocity - ind_pred) * 1000

      error_df <- data.frame(
        error = c(pop_error, ind_error),
        method = factor(
          rep(c("Population Average", "Individual Calibrated"), each = nrow(data)),
          levels = c("Population Average", "Individual Calibrated")
        )
      )

      # Calculate MAE for subtitle
      pop_mae <- mean(abs(pop_error))
      ind_mae <- mean(abs(ind_error))
      improvement <- (1 - ind_mae / pop_mae) * 100

      ggplot(error_df, aes(x = .data$error, fill = .data$method)) +
        geom_histogram(
          aes(y = after_stat(density)),
          bins = 30,
          alpha = 0.6,
          position = "identity"
        ) +
        geom_density(aes(color = .data$method), linewidth = 1, fill = NA) +
        geom_vline(xintercept = 0, linetype = "dashed", color = "gray30") +
        scale_fill_manual(values = c(
          "Population Average" = private$.colors$primary,
          "Individual Calibrated" = private$.colors$secondary
        )) +
        scale_color_manual(values = c(
          "Population Average" = private$.colors$primary,
          "Individual Calibrated" = private$.colors$secondary
        )) +
        labs(
          title = title,
          subtitle = sprintf(
            "MAE: Population = %.1f mm/s | Individual = %.1f mm/s (%.0f%% improvement)",
            pop_mae, ind_mae, improvement
          ),
          x = "Prediction Error (mm/s)",
          y = "Density",
          fill = "Approach",
          color = "Approach"
        ) +
        private$.get_base_theme() +
        theme(legend.position = "bottom")
    },

    #' @description Get prediction error statistics
    #' @param data Data frame with mean_velocity column
    #' @param model An lme4 lmer model object
    #' @return List with pop_mae, ind_mae, improvement_pct
    calculate_prediction_errors = function(data, model) {
      pop_pred <- predict(model, re.form = NA)
      ind_pred <- predict(model)

      pop_error <- data$mean_velocity - pop_pred
      ind_error <- data$mean_velocity - ind_pred

      pop_mae <- mean(abs(pop_error)) * 1000
      ind_mae <- mean(abs(ind_error)) * 1000

      list(
        pop_mae_mm = pop_mae,
        ind_mae_mm = ind_mae,
        pop_rmse_mm = sqrt(mean(pop_error^2)) * 1000,
        ind_rmse_mm = sqrt(mean(ind_error^2)) * 1000,
        improvement_pct = (1 - ind_mae / pop_mae) * 100
      )
    }
  ),

  private = list(
    .base_size = NULL,
    .colors = NULL,

    .get_base_theme = function() {
      theme_minimal(base_size = private$.base_size) +
        theme(
          plot.title = element_text(face = "bold", size = 14),
          plot.subtitle = element_text(color = private$.colors$gray, size = 10),
          panel.grid.minor = element_blank()
        )
    },

    .create_caterpillar_panel = function(data, title, x_label) {
      ggplot(data, aes(x = .data$estimate, y = .data$order)) +
        geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
        geom_segment(
          aes(x = .data$lower, xend = .data$upper, yend = .data$order),
          color = private$.colors$primary,
          alpha = 0.5
        ) +
        geom_point(color = private$.colors$primary, size = 2) +
        labs(
          title = title,
          x = x_label,
          y = "Participant (ordered)"
        ) +
        private$.get_base_theme() +
        theme(
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank()
        )
    }
  )
)

#' Convenience function to create model plot
#'
#' @param data Data frame with id, rir, mean_velocity
#' @param model An lme4 lmer model object
#' @param ... Additional arguments passed to ModelPlotter$plot_model
#' @return ggplot object
#' @export
plot_model <- function(data, model, ...) {
  plotter <- ModelPlotter$new()
  plotter$plot_model(data, model, ...)
}

#' Convenience function to create caterpillar plot
#'
#' @param model An lme4 lmer model object
#' @param ... Additional arguments passed to ModelPlotter$plot_caterpillar_dual
#' @return ggplot object
#' @export
plot_caterpillar <- function(model, ...) {
  plotter <- ModelPlotter$new()
  plotter$plot_caterpillar_dual(model, ...)
}

#' Convenience function to create prediction comparison
#'
#' @param data Data frame with mean_velocity
#' @param model An lme4 lmer model object
#' @param ... Additional arguments
#' @return ggplot object
#' @export
plot_prediction_comparison <- function(data, model, ...) {
  plotter <- ModelPlotter$new()
  plotter$plot_prediction_comparison(data, model, ...)
}
