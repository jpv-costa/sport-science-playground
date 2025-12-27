# _common.R
# Shared utilities for all Quarto reports in the Deadlift Study project
# Source this file in setup chunks: source("_common.R")

# ==============================================================================
# PACKAGES
# ==============================================================================
suppressPackageStartupMessages({
  library(ggplot2)
  library(knitr)
  library(dplyr)
  library(scales)
})

# ==============================================================================
# GLOBAL OPTIONS
# ==============================================================================
knitr::opts_chunk$set(

  echo = TRUE,
  message = FALSE,
  warning = FALSE,
  fig.width = 10,
  fig.height = 6,
  fig.retina = 2,
  out.width = "100%"
)

options(
  digits = 4,
  scipen = 999,
  knitr.kable.NA = ""
)

# ==============================================================================
# COLOR PALETTE (Colorblind-Friendly)
# ==============================================================================
COLORS <- list(
  # Primary accent colors
  primary = "#2E86AB",
  secondary = "#E63946",
  tertiary = "#F77F00",

 # Load percentage colors (from ColorBrewer, colorblind-safe)
  load_80 = "#E69F00",
  load_90 = "#009E73",

  # Outcome colors
  strength = "#377EB8",
  hypertrophy = "#E41A1C",

  # Significance colors
  significant = "#E41A1C",
  not_significant = "#999999",

  # Neutral
 gray = "#999999",
  dark_gray = "#333333",
  light_gray = "#CCCCCC"
)

# ==============================================================================
# THEME
# ==============================================================================

#' Create consistent theme for all plots
#'
#' @param base_size Base font size (default: 14)
#' @return ggplot2 theme object
theme_study <- function(base_size = 14) {
  theme_minimal(base_size = base_size) +
    theme(
      # Text
      plot.title = element_text(face = "bold", size = rel(1.1)),
      plot.subtitle = element_text(color = COLORS$gray, size = rel(0.9)),
      plot.caption = element_text(color = COLORS$gray, size = rel(0.8)),

      # Axes
      axis.title = element_text(size = rel(0.9)),
      axis.text = element_text(size = rel(0.85)),

      # Legend
      legend.position = "bottom",
      legend.title = element_text(size = rel(0.9)),
      legend.text = element_text(size = rel(0.85)),

      # Panel
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),

      # Strip (facets)
      strip.text = element_text(face = "bold", size = rel(0.95))
    )
}

# Set as default theme
theme_set(theme_study())

# ==============================================================================
# TABLE UTILITIES
# ==============================================================================

#' Create formatted kable table with consistent styling
#'
#' @param df Data frame to display
#' @param caption Table caption (optional)
#' @param digits Number of decimal places for numeric columns (default: 3)
#' @param col.names Optional column name overrides
#' @return kable object
format_table <- function(df, caption = NULL, digits = 3, col.names = NULL) {
  df <- df |>
    dplyr::mutate(dplyr::across(
      where(is.numeric),
      ~ round(.x, digits)
    ))

  knitr::kable(
    df,
    caption = caption,
    col.names = col.names,
    format.args = list(big.mark = ",")
  )
}

#' Create simple summary statistics table
#'
#' @param metrics Character vector of metric names
#' @param values Vector of values (will be converted to character)
#' @param caption Table caption (optional)
#' @return kable object
summary_table <- function(metrics, values, caption = NULL) {
  df <- data.frame(
    Metric = metrics,
    Value = as.character(values),
    check.names = FALSE
  )
  knitr::kable(df, caption = caption, col.names = c("", ""))
}

#' Create model comparison table
#'
#' @param models Character vector of model names
#' @param r2_values Numeric vector of R-squared values
#' @param interpretations Character vector of interpretations
#' @param caption Table caption (optional)
#' @return kable object
model_comparison_table <- function(models, r2_values, interpretations = NULL, caption = NULL) {
  df <- data.frame(
    `Model Type` = models,
    `R²` = round(r2_values, 3),
    check.names = FALSE
  )

  if (!is.null(interpretations)) {
    df$Interpretation <- interpretations
  }

  format_table(df, caption = caption)
}

# ==============================================================================
# FORMATTING HELPERS
# ==============================================================================

#' Format p-value for display
#'
#' @param p P-value
#' @param digits Decimal places for non-tiny values (default: 3)
#' @return Formatted string
format_p <- function(p, digits = 3) {
  if (is.na(p)) return("NA")
  if (p < 0.001) return("<0.001")
  sprintf(paste0("%.", digits, "f"), p)
}

#' Format estimate with confidence interval
#'
#' @param estimate Point estimate
#' @param lower Lower bound
#' @param upper Upper bound
#' @param digits Decimal places (default: 3)
#' @return Formatted string like "0.23 [0.15, 0.31]"
format_ci <- function(estimate, lower, upper, digits = 3) {
  sprintf(
    "%.*f [%.*f, %.*f]",
    digits, estimate,
    digits, lower,
    digits, upper
  )
}

#' Format number for inline reporting
#'
#' @param x Numeric value
#' @param digits Decimal places (default: 2)
#' @return Formatted string with thousand separators
inline_num <- function(x, digits = 2) {
  format(round(x, digits), nsmall = digits, big.mark = ",")
}

#' Format percentage for inline reporting
#'
#' @param x Numeric value (0-1 scale or 0-100 scale)
#' @param digits Decimal places (default: 1)
#' @param scale If TRUE, multiply by 100 (default: TRUE for 0-1 scale)
#' @return Formatted string with % suffix
inline_pct <- function(x, digits = 1, scale = TRUE) {
  if (scale && max(x, na.rm = TRUE) <= 1) {
    x <- x * 100
  }
  paste0(format(round(x, digits), nsmall = digits), "%")
}

# ==============================================================================
# COLOR SCALES
# ==============================================================================

#' Create color scale for load percentages
#' @return scale_color_manual object
scale_color_load <- function() {
  scale_color_manual(
    values = c("80%" = COLORS$load_80, "90%" = COLORS$load_90),
    name = "Load"
  )
}

#' Create fill scale for load percentages
#' @return scale_fill_manual object
scale_fill_load <- function() {
  scale_fill_manual(
    values = c("80%" = COLORS$load_80, "90%" = COLORS$load_90),
    name = "Load"
  )
}

#' Create color scale for significance
#' @return scale_color_manual object
scale_color_significance <- function() {
  scale_color_manual(
    values = c("Yes" = COLORS$significant, "No" = COLORS$not_significant),
    name = "Significant"
  )
}

#' Create color scale for outcomes (strength vs hypertrophy)
#' @return scale_color_manual object
scale_color_outcome <- function() {
  scale_color_manual(
    values = c("Strength" = COLORS$strength, "Hypertrophy" = COLORS$hypertrophy),
    name = "Outcome"
  )
}

# ==============================================================================
# DATA LOADING UTILITIES
# ==============================================================================

#' Safely load RDS file with error handling
#'
#' @param path Path to RDS file
#' @param fallback Value to return if file not found (default: NULL)
#' @return Loaded object or fallback
safe_load_rds <- function(path, fallback = NULL) {
  tryCatch(
    readRDS(path),
    error = function(e) {
      warning(sprintf("Could not load: %s. Error: %s", path, e$message))
      fallback
    }
  )
}

#' Load results from data/processed/ directory
#'
#' @param name Name of results file (without path or extension)
#' @return Loaded results object
load_results <- function(name) {
  path <- file.path("..", "data", "processed", paste0(name, ".rds"))
  safe_load_rds(path)
}

# ==============================================================================
# COMMON PLOT FUNCTIONS
# ==============================================================================

#' Create velocity-RIR scatter plot
#'
#' @param data Data frame with mean_velocity and rir columns
#' @param color_var Variable for color aesthetic (optional)
#' @param facet_var Variable for faceting (optional)
#' @param title Plot title
#' @return ggplot object
plot_velocity_rir <- function(data, color_var = NULL, facet_var = NULL, title = "Velocity-RIR Relationship") {
  p <- ggplot(data, aes(x = mean_velocity, y = rir))

  if (!is.null(color_var)) {
    p <- p + aes(color = .data[[color_var]])
  }

  p <- p +
    geom_point(alpha = 0.4, size = 2) +
    geom_smooth(method = "lm", formula = y ~ poly(x, 2), se = TRUE, linewidth = 1) +
    labs(
      title = title,
      x = "Mean Velocity (m/s)",
      y = "Repetitions in Reserve"
    )

  if (!is.null(facet_var)) {
    p <- p + facet_wrap(as.formula(paste("~", facet_var)))
  }

  p
}

#' Create coefficient comparison plot with confidence intervals
#'
#' @param data Data frame with columns: name, estimate, lower, upper, significant
#' @param title Plot title
#' @return ggplot object
plot_coefficients <- function(data, title = "Model Coefficients") {
  ggplot(data, aes(x = name, y = estimate, color = significant)) +
    geom_point(size = 4) +
    geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2, linewidth = 1) +
    geom_hline(yintercept = 0, linetype = "dashed", color = COLORS$gray) +
    scale_color_significance() +
    coord_flip() +
    labs(title = title, x = NULL, y = "Estimate")
}
