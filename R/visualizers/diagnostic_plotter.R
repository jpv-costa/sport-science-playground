# R/visualizers/diagnostic_plotter.R
# Service: Statistical Diagnostic Visualizations for VBT Research
#
# =============================================================================
# PRINCIPLES APPLIED (from CLAUDE.md)
# =============================================================================
#
# NAMING PRINCIPLES:
# - Understandability: Problem domain terms (bland_altman, caterpillar, qq_plot)
# - Consistency: All public methods start with "plot_", private with "."
# - Distinguishability: Specific names (plot_bland_altman vs plot_scatter)
# - Conciseness: Short but meaningful (LOA = limits of agreement)
#
# SOLID PRINCIPLES:
# - SRP: Single responsibility - only creates diagnostic visualizations
# - OCP: Open for extension via inheritance, closed for modification
# - DIP: Depends on abstractions (data frames), not concrete implementations
#
# CUPID PRINCIPLES:
# - Composable: Each plot method is independent, can be combined
# - Unix: Each method does one thing well (one plot type)
# - Predictable: Same inputs -> same outputs (pure rendering functions)
# - Idiomatic: Follows ggplot2/R conventions, proper roxygen docs
# - Domain-based: Names reflect statistical concepts, not implementation
#
# SCIENTIFIC VALIDITY:
# - Bland-Altman: Uses validated 1.96 SD for 95% limits of agreement
# - Q-Q plots: Standard approach for normality assessment
# - Bootstrap: Proper percentile method for CI visualization
#
# TESTABILITY:
# - All calculations extracted to pure private methods
# - Statistics methods return data structures, not plots
# - Plots can be tested via vdiffr or manual inspection
# =============================================================================

box::use(
  R6[R6Class],
  ggplot2[
    ggplot, aes, geom_point, geom_line, geom_hline, geom_vline,
    geom_smooth, geom_histogram, geom_density, geom_errorbar,
    geom_col, geom_text, geom_abline, geom_ribbon,
    stat_qq, stat_qq_line,
    labs, theme, theme_minimal, element_text, element_blank,
    annotate, coord_flip, coord_cartesian,
    scale_x_continuous, scale_fill_manual, scale_fill_identity, scale_color_manual,
    after_stat
  ],
  patchwork[...],
  dplyr[group_by, summarize, mutate, arrange, n_distinct, .data],
  tidyr[crossing],
  stats[lm, coef, predict, sd]
)

#' Color Palette for Visualizations
#'
#' Colorblind-friendly palette used across all plots.
#' Based on ColorBrewer recommendations.
#'
#' @export
PLOT_COLORS <- list(
  primary = "#2E86AB",
  secondary = "#E63946",
  tertiary = "#F77F00",
  gray = "#999999",
  dark_gray = "#333333",
  light_gray = "#CCCCCC"
)

#' Bland-Altman Statistics Result
#'
#' Value object containing Bland-Altman analysis results.
#' Testable: Can verify calculations without rendering plots.
#'
#' @export
BlandAltmanResult <- R6Class(
  classname = "BlandAltmanResult",
  cloneable = FALSE,

  public = list(
    #' @field mean_diff Mean difference (bias)
    mean_diff = NULL,

    #' @field sd_diff Standard deviation of differences
    sd_diff = NULL,

    #' @field upper_loa Upper limit of agreement (+1.96 SD)
    upper_loa = NULL,

    #' @field lower_loa Lower limit of agreement (-1.96 SD)
    lower_loa = NULL,

    #' @field n Sample size
    n = NULL,

    #' @description Create Bland-Altman result
    #' @param day1 Vector of Day 1 measurements
    #' @param day2 Vector of Day 2 measurements
    initialize = function(day1, day2) {
      stopifnot(
        "day1 must be numeric" = is.numeric(day1),
        "day2 must be numeric" = is.numeric(day2),
        "day1 and day2 must have same length" = length(day1) == length(day2),
        "Need at least 3 observations" = length(day1) >= 3
      )

      differences <- day1 - day2
      self$mean_diff <- mean(differences)
      self$sd_diff <- sd(differences)
      self$upper_loa <- self$mean_diff + 1.96 * self$sd_diff
      self$lower_loa <- self$mean_diff - 1.96 * self$sd_diff
      self$n <- length(day1)
    },

    #' @description Check if bias is clinically significant
    #' @param threshold Maximum acceptable bias
    #' @return Logical
    is_bias_acceptable = function(threshold) {
      abs(self$mean_diff) < threshold
    },

    #' @description Get limits of agreement width
    #' @return Numeric
    get_loa_width = function() {
      self$upper_loa - self$lower_loa
    },

    #' @description Convert to list for serialization
    #' @return List
    to_list = function() {
      list(
        mean_diff = self$mean_diff,
        sd_diff = self$sd_diff,
        upper_loa = self$upper_loa,
        lower_loa = self$lower_loa,
        n = self$n,
        loa_width = self$get_loa_width()
      )
    }
  )
)

#' Spaghetti Plot Data Generator
#'
#' Prepares data for individual trajectory visualization.
#' Separated from plotting for testability.
#'
#' @export
SpaghettiDataGenerator <- R6Class(
  classname = "SpaghettiDataGenerator",
  cloneable = FALSE,

  public = list(
    #' @description Generate individual regression line data
    #' @param data Data frame with id, rir, mean_velocity columns
    #' @param n_points Number of points per line (default: 50)
    #' @return List with individual_lines and population_line data frames
    generate = function(data, n_points = 50) {
      stopifnot(
        "data must have 'id' column" = "id" %in% names(data),
        "data must have 'rir' column" = "rir" %in% names(data),
        "data must have 'mean_velocity' column" = "mean_velocity" %in% names(data)
      )

      individual_fits <- private$.calculate_individual_fits(data)
      rir_range <- seq(min(data$rir), max(data$rir), length.out = n_points)

      individual_lines <- private$.generate_individual_lines(
        individual_fits, rir_range
      )
      population_line <- private$.generate_population_line(data, rir_range)

      list(
        individual_lines = individual_lines,
        population_line = population_line,
        n_individuals = n_distinct(data$id),
        rir_range = range(data$rir)
      )
    }
  ),

  private = list(
    .calculate_individual_fits = function(data) {
      data |>
        group_by(.data$id) |>
        summarize(
          intercept = coef(lm(mean_velocity ~ rir))[1],
          slope = coef(lm(mean_velocity ~ rir))[2],
          .groups = "drop"
        )
    },

    .generate_individual_lines = function(fits, rir_range) {
      fits |>
        crossing(rir = rir_range) |>
        mutate(predicted_velocity = .data$intercept + .data$slope * .data$rir)
    },

    .generate_population_line = function(data, rir_range) {
      pop_model <- lm(mean_velocity ~ rir, data = data)
      data.frame(
        rir = rir_range,
        predicted_velocity = predict(
          pop_model,
          newdata = data.frame(rir = rir_range)
        )
      )
    }
  )
)

#' Diagnostic Plotter
#'
#' R6 class for creating statistical diagnostic visualizations.
#' All statistical calculations are delegated to testable helper classes.
#'
#' @section Principles:
#' - SRP: Only renders visualizations, calculations done by helper classes
#' - Testable: Plotting logic separated from statistical calculations
#' - Composable: Individual plots can be combined with patchwork
#'
#' @export
DiagnosticPlotter <- R6Class(
  classname = "DiagnosticPlotter",
  cloneable = FALSE,

  public = list(
    #' @description Create a new DiagnosticPlotter instance
    #' @param base_size Base font size for plots (default: 14)
    #' @param colors Named list of colors (optional, uses defaults)
    initialize = function(base_size = 14, colors = NULL) {
      private$.base_size <- base_size
      private$.colors <- colors %||% PLOT_COLORS
      private$.spaghetti_generator <- SpaghettiDataGenerator$new()
    },

    #' @description Create spaghetti plot of individual trajectories
    #' @param data Data frame with id, rir, mean_velocity columns
    #' @param highlight_population Show population average line
    #' @param title Plot title
    #' @return ggplot object
    plot_spaghetti = function(data,
                              highlight_population = TRUE,
                              title = "Individual Velocity-RIR Trajectories") {
      plot_data <- private$.spaghetti_generator$generate(data)

      p <- ggplot() +
        geom_line(
          data = plot_data$individual_lines,
          aes(x = .data$rir, y = .data$predicted_velocity, group = .data$id),
          color = private$.colors$gray,
          alpha = 0.4,
          linewidth = 0.5
        ) +
        geom_point(
          data = data,
          aes(x = .data$rir, y = .data$mean_velocity),
          color = private$.colors$primary,
          alpha = 0.2,
          size = 1
        )

      if (highlight_population) {
        p <- p +
          geom_line(
            data = plot_data$population_line,
            aes(x = .data$rir, y = .data$predicted_velocity),
            color = private$.colors$secondary,
            linewidth = 2
          ) +
          annotate(
            "text",
            x = max(plot_data$population_line$rir) - 0.5,
            y = max(plot_data$population_line$predicted_velocity) + 0.02,
            label = "Population Average",
            color = private$.colors$secondary,
            fontface = "bold",
            hjust = 1,
            size = 4
          )
      }

      p +
        labs(
          title = title,
          subtitle = sprintf("n = %d individuals", plot_data$n_individuals),
          x = "Repetitions in Reserve (RIR)",
          y = "Mean Velocity (m/s)"
        ) +
        scale_x_continuous(breaks = seq(
          plot_data$rir_range[1],
          plot_data$rir_range[2]
        )) +
        private$.get_base_theme()
    },

    #' @description Create Q-Q plot and residuals vs fitted panel
    #' @param residuals Vector of model residuals
    #' @param fitted Vector of fitted values
    #' @return Combined ggplot object
    plot_diagnostics_panel = function(residuals, fitted) {
      stopifnot(
        "residuals must be numeric" = is.numeric(residuals),
        "fitted must be numeric" = is.numeric(fitted),
        "residuals and fitted must have same length" =
          length(residuals) == length(fitted)
      )

      diag_data <- data.frame(residuals = residuals, fitted = fitted)

      p_qq <- private$.create_qq_plot(diag_data)
      p_resid <- private$.create_residual_plot(diag_data)

      p_qq + p_resid +
        patchwork::plot_annotation(
          title = "Model Diagnostic Plots",
          theme = theme(plot.title = element_text(face = "bold", size = 16))
        )
    },

    #' @description Create Bland-Altman plot for reliability assessment
    #' @param day1 Vector of Day 1 measurements
    #' @param day2 Vector of Day 2 measurements
    #' @param title Plot title
    #' @param y_label Label for measurement type
    #' @return ggplot object
    plot_bland_altman = function(day1, day2,
                                  title = "Bland-Altman Plot",
                                  y_label = "Measurement") {
      ba_result <- BlandAltmanResult$new(day1, day2)

      ba_data <- data.frame(
        mean_value = (day1 + day2) / 2,
        difference = day1 - day2
      )

      private$.build_bland_altman_plot(ba_data, ba_result, title, y_label)
    },

    #' @description Create caterpillar plot for random effects
    #' @param random_effects Data frame with id, estimate, lower, upper
    #' @param title Plot title
    #' @return ggplot object
    plot_caterpillar = function(random_effects,
                                title = "Individual Random Effects") {
      stopifnot(
        "random_effects must have 'id' column" =
          "id" %in% names(random_effects),
        "random_effects must have 'estimate' column" =
          "estimate" %in% names(random_effects),
        "random_effects must have 'lower' column" =
          "lower" %in% names(random_effects),
        "random_effects must have 'upper' column" =
          "upper" %in% names(random_effects)
      )

      ordered_data <- random_effects |>
        arrange(.data$estimate) |>
        mutate(id = factor(.data$id, levels = .data$id))

      ggplot(ordered_data, aes(x = .data$id, y = .data$estimate)) +
        geom_hline(
          yintercept = 0,
          linetype = "dashed",
          color = private$.colors$gray
        ) +
        geom_errorbar(
          aes(ymin = .data$lower, ymax = .data$upper),
          width = 0,
          color = private$.colors$primary,
          alpha = 0.6
        ) +
        geom_point(color = private$.colors$secondary, size = 2) +
        coord_flip() +
        labs(
          title = title,
          subtitle = "Deviation from population average (with 95% CI)",
          x = "Participant",
          y = "Random Effect Estimate"
        ) +
        private$.get_base_theme()
    },

    #' @description Create bootstrap distribution plot
    #' @param bootstrap_samples Vector of bootstrap samples
    #' @param observed_value Observed estimate
    #' @param ci_lower Lower CI bound
    #' @param ci_upper Upper CI bound
    #' @param title Plot title
    #' @param x_label X-axis label
    #' @return ggplot object
    plot_bootstrap_distribution = function(bootstrap_samples,
                                            observed_value,
                                            ci_lower,
                                            ci_upper,
                                            title = "Bootstrap Distribution",
                                            x_label = "Estimate") {
      stopifnot(
        "bootstrap_samples must be numeric" = is.numeric(bootstrap_samples),
        "Need at least 100 bootstrap samples" = length(bootstrap_samples) >= 100,
        "observed_value must be numeric" = is.numeric(observed_value),
        "ci_lower must be numeric" = is.numeric(ci_lower),
        "ci_upper must be numeric" = is.numeric(ci_upper),
        "ci_lower must be less than ci_upper" = ci_lower < ci_upper
      )

      boot_data <- data.frame(value = bootstrap_samples)

      ggplot(boot_data, aes(x = .data$value)) +
        annotate(
          "rect",
          xmin = ci_lower, xmax = ci_upper,
          ymin = -Inf, ymax = Inf,
          fill = private$.colors$tertiary, alpha = 0.2
        ) +
        geom_histogram(
          aes(y = after_stat(density)),
          bins = 50,
          fill = private$.colors$primary,
          alpha = 0.7,
          color = "white"
        ) +
        geom_density(color = private$.colors$secondary, linewidth = 1) +
        geom_vline(
          xintercept = c(ci_lower, ci_upper),
          linetype = "dashed",
          color = private$.colors$tertiary,
          linewidth = 1
        ) +
        geom_vline(
          xintercept = observed_value,
          color = private$.colors$secondary,
          linewidth = 1.5
        ) +
        labs(
          title = title,
          subtitle = sprintf("95%% CI: [%.4f, %.4f]", ci_lower, ci_upper),
          x = x_label,
          y = "Density"
        ) +
        private$.get_base_theme()
    },

    #' @description Create velocity zones chart for practitioners
    #' @param velocity_thresholds Data frame with rir, velocity, lower, upper
    #' @param title Plot title
    #' @return ggplot object
    plot_velocity_zones = function(velocity_thresholds,
                                   title = "Velocity Zones by RIR") {
      stopifnot(
        "velocity_thresholds must have 'rir' column" =
          "rir" %in% names(velocity_thresholds),
        "velocity_thresholds must have 'velocity' column" =
          "velocity" %in% names(velocity_thresholds)
      )

      # Create zone data with colors
      zone_data <- velocity_thresholds |>
        mutate(
          zone_label = paste0("RIR ", .data$rir),
          zone_color = private$.get_zone_color(.data$rir)
        )

      # Bar chart with CI whiskers if available
      p <- ggplot(zone_data, aes(
        x = factor(.data$rir),
        y = .data$velocity,
        fill = .data$zone_color
      )) +
        geom_col(alpha = 0.8, width = 0.7)

      # Add CI if present
      if ("lower" %in% names(zone_data) && "upper" %in% names(zone_data)) {
        p <- p + geom_errorbar(
          aes(ymin = .data$lower, ymax = .data$upper),
          width = 0.2,
          color = private$.colors$dark_gray
        )
      }

      p +
        scale_fill_identity() +
        labs(
          title = title,
          subtitle = "Mean velocity thresholds for each RIR level",
          x = "Repetitions in Reserve (RIR)",
          y = "Mean Velocity (m/s)"
        ) +
        private$.get_base_theme() +
        theme(legend.position = "none")
    },

    #' @description Create coverage calibration plot
    #' @param calibration_data Data frame with target_coverage, empirical_coverage
    #' @param title Plot title
    #' @return ggplot object
    plot_coverage_calibration = function(calibration_data,
                                          title = "Conformal Prediction Calibration") {
      stopifnot(
        "calibration_data must have 'target_coverage' column" =
          "target_coverage" %in% names(calibration_data),
        "calibration_data must have 'empirical_coverage' column" =
          "empirical_coverage" %in% names(calibration_data)
      )

      ggplot(calibration_data, aes(
        x = .data$target_coverage,
        y = .data$empirical_coverage
      )) +
        geom_abline(
          slope = 1, intercept = 0,
          linetype = "dashed",
          color = private$.colors$gray
        ) +
        geom_point(
          color = private$.colors$primary,
          size = 3
        ) +
        geom_line(
          color = private$.colors$primary,
          linewidth = 1
        ) +
        annotate(
          "text",
          x = 0.95, y = 0.5,
          label = "Perfect\nCalibration",
          color = private$.colors$gray,
          fontface = "italic",
          hjust = 1
        ) +
        coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
        labs(
          title = title,
          subtitle = "Points on diagonal indicate well-calibrated prediction intervals",
          x = "Target Coverage",
          y = "Empirical Coverage"
        ) +
        private$.get_base_theme()
    },

    #' @description Create exercise comparison plot (e.g., deadlift vs squat)
    #' @param comparison_data Data frame with exercise, slope, slope_se
    #' @param title Plot title
    #' @return ggplot object
    plot_exercise_comparison = function(comparison_data,
                                         title = "Velocity-RIR Slope by Exercise") {
      stopifnot(
        "comparison_data must have 'exercise' column" =
          "exercise" %in% names(comparison_data),
        "comparison_data must have 'slope' column" =
          "slope" %in% names(comparison_data),
        "comparison_data must have 'slope_se' column" =
          "slope_se" %in% names(comparison_data)
      )

      plot_data <- comparison_data |>
        mutate(
          lower = .data$slope - 1.96 * .data$slope_se,
          upper = .data$slope + 1.96 * .data$slope_se
        )

      ggplot(plot_data, aes(
        x = .data$exercise,
        y = .data$slope,
        fill = .data$exercise
      )) +
        geom_col(alpha = 0.8, width = 0.6) +
        geom_errorbar(
          aes(ymin = .data$lower, ymax = .data$upper),
          width = 0.15,
          color = private$.colors$dark_gray
        ) +
        scale_fill_manual(values = c(
          private$.colors$primary,
          private$.colors$secondary
        )) +
        labs(
          title = title,
          subtitle = "Slope of velocity vs RIR relationship (with 95% CI)",
          x = "Exercise",
          y = "Slope (m/s per RIR)"
        ) +
        private$.get_base_theme() +
        theme(legend.position = "none")
    }
  ),

  private = list(
    .base_size = NULL,
    .colors = NULL,
    .spaghetti_generator = NULL,

    .get_base_theme = function() {
      theme_minimal(base_size = private$.base_size) +
        theme(
          plot.title = element_text(face = "bold", size = 16),
          plot.subtitle = element_text(color = private$.colors$gray),
          panel.grid.minor = element_blank()
        )
    },

    .create_qq_plot = function(data) {
      ggplot(data, aes(sample = .data$residuals)) +
        stat_qq(color = private$.colors$primary, alpha = 0.6) +
        stat_qq_line(color = private$.colors$secondary, linewidth = 1) +
        labs(
          title = "Q-Q Plot of Residuals",
          subtitle = "Points should follow the line for normality",
          x = "Theoretical Quantiles",
          y = "Sample Quantiles"
        ) +
        private$.get_base_theme()
    },

    .create_residual_plot = function(data) {
      ggplot(data, aes(x = .data$fitted, y = .data$residuals)) +
        geom_hline(
          yintercept = 0,
          linetype = "dashed",
          color = private$.colors$gray
        ) +
        geom_point(color = private$.colors$primary, alpha = 0.5) +
        geom_smooth(
          method = "loess",
          formula = y ~ x,
          color = private$.colors$secondary,
          se = FALSE,
          linewidth = 1
        ) +
        labs(
          title = "Residuals vs Fitted Values",
          subtitle = "Look for patterns (heteroscedasticity)",
          x = "Fitted Values",
          y = "Residuals"
        ) +
        private$.get_base_theme()
    },

    .build_bland_altman_plot = function(data, ba_result, title, y_label) {
      ggplot(data, aes(x = .data$mean_value, y = .data$difference)) +
        geom_hline(
          yintercept = ba_result$upper_loa,
          linetype = "dashed",
          color = private$.colors$secondary
        ) +
        geom_hline(
          yintercept = ba_result$lower_loa,
          linetype = "dashed",
          color = private$.colors$secondary
        ) +
        geom_hline(
          yintercept = ba_result$mean_diff,
          color = private$.colors$primary,
          linewidth = 1
        ) +
        geom_point(color = private$.colors$primary, alpha = 0.7, size = 3) +
        annotate(
          "text",
          x = max(data$mean_value),
          y = ba_result$upper_loa,
          label = sprintf("+1.96 SD (%.3f)", ba_result$upper_loa),
          hjust = 1, vjust = -0.5,
          color = private$.colors$secondary, size = 3.5
        ) +
        annotate(
          "text",
          x = max(data$mean_value),
          y = ba_result$lower_loa,
          label = sprintf("-1.96 SD (%.3f)", ba_result$lower_loa),
          hjust = 1, vjust = 1.5,
          color = private$.colors$secondary, size = 3.5
        ) +
        annotate(
          "text",
          x = max(data$mean_value),
          y = ba_result$mean_diff,
          label = sprintf("Bias: %.3f", ba_result$mean_diff),
          hjust = 1, vjust = -0.5,
          color = private$.colors$primary, fontface = "bold", size = 3.5
        ) +
        labs(
          title = title,
          subtitle = sprintf(
            "Limits of Agreement: [%.3f, %.3f] | Bias: %.3f",
            ba_result$lower_loa, ba_result$upper_loa, ba_result$mean_diff
          ),
          x = paste("Mean of Day 1 and Day 2", y_label),
          y = paste("Difference (Day 1 - Day 2)", y_label)
        ) +
        private$.get_base_theme()
    },

    .get_zone_color = function(rir) {
      # Color gradient from red (RIR 0) to green (RIR 7+)
      # Based on traffic light metaphor for effort zones
      colors <- c(
        "#D32F2F",  # 0 - Red (failure)
        "#E64A19",  # 1 - Deep orange
        "#F57C00",  # 2 - Orange
        "#FFA000",  # 3 - Amber
        "#FFC107",  # 4 - Yellow
        "#CDDC39",  # 5 - Lime
        "#8BC34A",  # 6 - Light green
        "#4CAF50"   # 7+ - Green (easy)
      )
      sapply(rir, function(r) colors[min(r + 1, length(colors))])
    }
  )
)
