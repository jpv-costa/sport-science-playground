# R/calculators/influence_diagnostics.R
# Influence Diagnostics for Mixed Effects Models
#
# =============================================================================
# PRINCIPLES APPLIED (from CLAUDE.md)
# =============================================================================
#
# NAMING PRINCIPLES:
# - Understandability: Statistical terms (leverage, cooks_distance)
# - Consistency: calculate_* for computations, plot_* for visualizations
# - Distinguishability: Separate observation vs participant level
# - Conciseness: Clear abbreviations where standard (DFBETAS)
#
# SOLID PRINCIPLES:
# - SRP: Single responsibility - only influence diagnostics
# - OCP: Open for extension via new diagnostic measures
#
# EDUCATIONAL FOCUS:
# - Plain language explanations in method documentation
# - Practical interpretation guidance
# =============================================================================

box::use(
  R6[R6Class],
  stats[hatvalues, cooks.distance, residuals, coef, lm, predict, sd, quantile, sigma],
  ggplot2[ggplot, aes, geom_point, geom_hline, geom_vline, geom_text,
          geom_segment, labs, theme_minimal, theme, element_text,
          scale_color_gradient2, scale_size_continuous, coord_cartesian]
)

#' Influence Diagnostics Result
#'
#' Immutable value object containing influence diagnostic measures
#'
#' @export
InfluenceDiagnosticsResult <- R6Class(
 classname = "InfluenceDiagnosticsResult",
  cloneable = FALSE,

  public = list(
    #' @field leverage Vector of leverage values (hat values)
    leverage = NULL,
    #' @field cooks_d Vector of Cook's distance values
    cooks_d = NULL,
    #' @field standardized_residuals Vector of standardized residuals
    standardized_residuals = NULL,
    #' @field n_observations Number of observations
    n_observations = NULL,
    #' @field influential_threshold Threshold for influential points
    influential_threshold = NULL,
    #' @field diagnostics_data Data frame with all measures
    diagnostics_data = NULL,

    #' @description Create a new diagnostics result
    initialize = function(leverage, cooks_d, standardized_residuals,
                          n_observations, influential_threshold, diagnostics_data) {
      self$leverage <- leverage
      self$cooks_d <- cooks_d
      self$standardized_residuals <- standardized_residuals
      self$n_observations <- n_observations
      self$influential_threshold <- influential_threshold
      self$diagnostics_data <- diagnostics_data
    },

    #' @description Count influential observations
    count_influential = function() {
      sum(self$cooks_d > self$influential_threshold, na.rm = TRUE)
    },

    #' @description Get indices of influential observations
    get_influential_indices = function() {
      which(self$cooks_d > self$influential_threshold)
    },

    #' @description Convert to list for serialization
    to_list = function() {
      list(
        n_observations = self$n_observations,
        n_influential = self$count_influential(),
        influential_threshold = self$influential_threshold,
        max_cooks_d = max(self$cooks_d, na.rm = TRUE),
        max_leverage = max(self$leverage, na.rm = TRUE),
        mean_leverage = mean(self$leverage, na.rm = TRUE)
      )
    }
  )
)

#' Participant Influence Result
#'
#' Immutable value object for participant-level influence analysis
#'
#' @export
ParticipantInfluenceResult <- R6Class(
  classname = "ParticipantInfluenceResult",
  cloneable = FALSE,

  public = list(
    #' @field participant_effects Data frame with per-participant influence
    participant_effects = NULL,
    #' @field coefficient_name Name of the coefficient being analyzed
    coefficient_name = NULL,
    #' @field original_estimate Original coefficient estimate
    original_estimate = NULL,
    #' @field max_change Maximum change when removing a participant
    max_change = NULL,
    #' @field most_influential_id ID of most influential participant
    most_influential_id = NULL,

    #' @description Create a new participant influence result
    initialize = function(participant_effects, coefficient_name, original_estimate,
                          max_change, most_influential_id) {
      self$participant_effects <- participant_effects
      self$coefficient_name <- coefficient_name
      self$original_estimate <- original_estimate
      self$max_change <- max_change
      self$most_influential_id <- most_influential_id
    },

    #' @description Interpret influence magnitude
    interpret = function() {
      pct_change <- abs(self$max_change / self$original_estimate) * 100
      if (pct_change < 5) {
        "minimal"
      } else if (pct_change < 10) {
        "small"
      } else if (pct_change < 20) {
        "moderate"
      } else {
        "substantial"
      }
    },

    #' @description Convert to list for serialization
    to_list = function() {
      list(
        coefficient_name = self$coefficient_name,
        original_estimate = self$original_estimate,
        max_change = self$max_change,
        pct_change = abs(self$max_change / self$original_estimate) * 100,
        most_influential_id = self$most_influential_id,
        interpretation = self$interpret()
      )
    }
  )
)

#' Influence Diagnostics Analyzer
#'
#' R6 class for calculating and visualizing influence diagnostics
#' for mixed effects models.
#'
#' @section What is Influence?
#' Influence measures how much each observation (or participant) affects
#' the model's conclusions. High-influence points could be:
#' - Outliers with unusual Y values
#' - High-leverage points with unusual X values
#' - Or both
#'
#' @section Why It Matters:
#' If your conclusions depend heavily on a few unusual observations,
#' you should investigate those cases and consider sensitivity analyses.
#'
#' @export
InfluenceDiagnostics <- R6Class(
  classname = "InfluenceDiagnostics",
  cloneable = FALSE,

  public = list(
    #' @description Create a new InfluenceDiagnostics instance
    initialize = function() {
      # No initialization needed
    },

    # =========================================================================
    # OBSERVATION-LEVEL DIAGNOSTICS
    # =========================================================================

    #' @description Calculate Observation-Level Influence Measures
    #'
    #' Computes leverage, Cook's distance, and standardized residuals
    #' for each observation.
    #'
    #' @section Interpretation Guide:
    #' - **Leverage**: How unusual is this observation's predictor values?
    #'   High leverage means the point has "pull" on the regression line.
    #'   Threshold: > 2*(p+1)/n where p = number of predictors
    #'
    #' - **Cook's Distance**: How much would results change without this point?
    #'   Combines leverage and residual size.
    #'   Threshold: > 4/n (common rule of thumb)
    #'
    #' - **Standardized Residuals**: How unusual is the Y value?
    #'   Values > |2| are worth investigating.
    #'
    #' @param model Fitted lme4 model
    #' @return InfluenceDiagnosticsResult object
    calculate_observation_influence = function(model) {
      if (!requireNamespace("lme4", quietly = TRUE)) {
        stop("Package 'lme4' required for calculate_observation_influence()")
      }

      # Get hat values (leverage)
      leverage <- hatvalues(model)

      # Get residuals
      resid_raw <- residuals(model)
      sigma_est <- sigma(model)
      std_resid <- resid_raw / sigma_est

      # Calculate Cook's distance approximation for mixed models
      # Cook's D = (residual^2 / p*MSE) * (leverage / (1-leverage)^2)
      n <- length(resid_raw)
      p <- length(lme4::fixef(model))  # Number of fixed effects

      # Simplified Cook's D calculation
      cooks_d <- private$.calculate_cooks_d(model, leverage, resid_raw, p)

      # Standard threshold: 4/n
      threshold <- 4 / n

      # Create diagnostics data frame
      diagnostics_data <- data.frame(
        observation = seq_along(leverage),
        leverage = leverage,
        cooks_d = cooks_d,
        std_residual = std_resid,
        is_influential = cooks_d > threshold
      )

      InfluenceDiagnosticsResult$new(
        leverage = leverage,
        cooks_d = cooks_d,
        standardized_residuals = std_resid,
        n_observations = n,
        influential_threshold = threshold,
        diagnostics_data = diagnostics_data
      )
    },

    # =========================================================================
    # PARTICIPANT-LEVEL INFLUENCE
    # =========================================================================

    #' @description Calculate Participant-Level Influence
    #'
    #' Leave-one-participant-out analysis: How much does removing each
    #' participant change the key coefficient estimate?
    #'
    #' @section Why Participant-Level?
    #' In clustered data, removing individual observations doesn't capture
    #' the full influence of a participant. This method removes ALL
    #' observations from each participant and refits the model.
    #'
    #' @param model Fitted lme4 model
    #' @param data Data frame used for fitting
    #' @param id_col Name of participant ID column
    #' @param coef_name Name of coefficient to track (default: second fixed effect)
    #' @return ParticipantInfluenceResult object
    calculate_participant_influence = function(model, data, id_col = "id",
                                                coef_name = NULL) {
      if (!requireNamespace("lme4", quietly = TRUE)) {
        stop("Package 'lme4' required for calculate_participant_influence()")
      }

      # Get original estimate
      fixed_effects <- lme4::fixef(model)
      if (is.null(coef_name)) {
        coef_name <- names(fixed_effects)[2]  # Usually the main predictor
      }
      original_estimate <- fixed_effects[coef_name]

      # Get formula
      formula_obj <- stats::formula(model)

      # Get unique participants
      participants <- unique(data[[id_col]])
      n_participants <- length(participants)

      # Store results
      loo_estimates <- numeric(n_participants)
      loo_change <- numeric(n_participants)

      for (i in seq_along(participants)) {
        # Remove this participant
        loo_data <- data[data[[id_col]] != participants[i], ]

        # Refit model
        loo_model <- tryCatch({
          lme4::lmer(formula_obj, data = loo_data, REML = FALSE)
        }, error = function(e) NULL)

        if (is.null(loo_model)) {
          loo_estimates[i] <- NA
          loo_change[i] <- NA
        } else {
          loo_fe <- lme4::fixef(loo_model)
          loo_estimates[i] <- loo_fe[coef_name]
          loo_change[i] <- loo_estimates[i] - original_estimate
        }
      }

      # Create participant effects table
      participant_effects <- data.frame(
        participant_id = participants,
        estimate_without = loo_estimates,
        change = loo_change,
        pct_change = (loo_change / original_estimate) * 100,
        stringsAsFactors = FALSE
      )

      # Order by absolute change
      participant_effects <- participant_effects[order(-abs(participant_effects$change)), ]

      # Find most influential
      max_idx <- which.max(abs(loo_change))

      ParticipantInfluenceResult$new(
        participant_effects = participant_effects,
        coefficient_name = coef_name,
        original_estimate = original_estimate,
        max_change = loo_change[max_idx],
        most_influential_id = as.character(participants[max_idx])
      )
    },

    # =========================================================================
    # VISUALIZATIONS
    # =========================================================================

    #' @description Plot Leverage vs Residuals
    #'
    #' Classic diagnostic plot showing leverage on x-axis, standardized
    #' residuals on y-axis, with point size proportional to Cook's D.
    #'
    #' @param diagnostics_result InfluenceDiagnosticsResult object
    #' @param title Plot title
    #' @param label_threshold Only label points with Cook's D above this
    #' @return ggplot object
    plot_leverage_residual = function(diagnostics_result,
                                      title = "Influence Diagnostics: Leverage vs Residuals",
                                      label_threshold = NULL) {
      diag_data <- diagnostics_result$diagnostics_data

      if (is.null(label_threshold)) {
        label_threshold <- diagnostics_result$influential_threshold
      }

      # Add labels for influential points
      diag_data$label <- ifelse(
        diag_data$cooks_d > label_threshold,
        diag_data$observation,
        NA
      )

      # Reference lines
      n <- diagnostics_result$n_observations

      ggplot(diag_data, aes(x = leverage, y = std_residual)) +
        # Reference lines
        geom_hline(yintercept = c(-2, 0, 2), linetype = c("dashed", "solid", "dashed"),
                   color = c("#E63946", "#666666", "#E63946"), alpha = 0.7) +
        # Points sized by Cook's D
        geom_point(aes(size = cooks_d, color = cooks_d), alpha = 0.7) +
        # Labels for influential points
        geom_text(aes(label = label), hjust = -0.2, vjust = 0.5, size = 3,
                  na.rm = TRUE, color = "#1D3557") +
        # Color scale
        scale_color_gradient2(
          low = "#2E86AB", mid = "#F77F00", high = "#E63946",
          midpoint = diagnostics_result$influential_threshold,
          name = "Cook's D"
        ) +
        scale_size_continuous(name = "Cook's D", range = c(1, 6)) +
        labs(
          title = title,
          subtitle = paste0(
            diagnostics_result$count_influential(), " influential observations ",
            "(Cook's D > ", round(diagnostics_result$influential_threshold, 4), ")"
          ),
          x = "Leverage (hat value)",
          y = "Standardized Residual",
          caption = "Dashed lines at ±2 standard deviations"
        ) +
        theme_minimal() +
        theme(
          plot.title = element_text(face = "bold", size = 14),
          plot.subtitle = element_text(size = 10, color = "#666666"),
          legend.position = "right"
        )
    },

    #' @description Plot Participant Influence
    #'
    #' Shows how much each participant's removal changes the key coefficient.
    #'
    #' @param influence_result ParticipantInfluenceResult object
    #' @param title Plot title
    #' @param top_n Number of most influential participants to highlight
    #' @return ggplot object
    plot_participant_influence = function(influence_result,
                                          title = "Participant Influence Analysis",
                                          top_n = 5) {
      effect_data <- influence_result$participant_effects

      # Mark top N
      effect_data$is_top <- seq_len(nrow(effect_data)) <= top_n

      # Add color based on direction
      effect_data$direction <- ifelse(effect_data$change > 0, "increases", "decreases")

      ggplot(effect_data, aes(x = reorder(participant_id, abs(change)),
                               y = pct_change)) +
        geom_hline(yintercept = 0, linetype = "solid", color = "#666666") +
        geom_hline(yintercept = c(-10, 10), linetype = "dashed",
                   color = "#E63946", alpha = 0.5) +
        geom_segment(aes(xend = participant_id, y = 0, yend = pct_change,
                         color = is_top),
                     linewidth = 1) +
        geom_point(aes(color = is_top), size = 3) +
        scale_color_manual(values = c("TRUE" = "#E63946", "FALSE" = "#2E86AB"),
                           guide = "none") +
        coord_flip() +
        labs(
          title = title,
          subtitle = paste0(
            "Effect on '", influence_result$coefficient_name,
            "' coefficient (original = ",
            round(influence_result$original_estimate, 4), ")"
          ),
          x = "Participant",
          y = "% Change in Coefficient When Removed",
          caption = paste0(
            "Most influential: ", influence_result$most_influential_id,
            " (", round(abs(influence_result$max_change / influence_result$original_estimate) * 100, 1),
            "% change)"
          )
        ) +
        theme_minimal() +
        theme(
          plot.title = element_text(face = "bold", size = 14),
          plot.subtitle = element_text(size = 10, color = "#666666"),
          axis.text.y = element_text(size = 8)
        )
    },

    #' @description Create Summary Table of Influential Cases
    #'
    #' Returns a formatted table of the most influential observations
    #' and participants for reporting.
    #'
    #' @param obs_diagnostics InfluenceDiagnosticsResult object
    #' @param participant_influence ParticipantInfluenceResult object (optional)
    #' @param top_n Number of cases to include
    #' @return Data frame suitable for reporting
    create_influence_summary = function(obs_diagnostics, participant_influence = NULL,
                                        top_n = 10) {
      # Top observations by Cook's D
      obs_data <- obs_diagnostics$diagnostics_data
      top_obs <- obs_data[order(-obs_data$cooks_d), ][1:min(top_n, nrow(obs_data)), ]

      summary_list <- list(
        observation_summary = data.frame(
          observation = top_obs$observation,
          leverage = round(top_obs$leverage, 4),
          cooks_d = round(top_obs$cooks_d, 4),
          std_residual = round(top_obs$std_residual, 2),
          influential = top_obs$is_influential
        )
      )

      if (!is.null(participant_influence)) {
        summary_list$participant_summary <- participant_influence$participant_effects[
          1:min(top_n, nrow(participant_influence$participant_effects)),
          c("participant_id", "change", "pct_change")
        ]
        summary_list$participant_summary$pct_change <- round(
          summary_list$participant_summary$pct_change, 2
        )
      }

      summary_list
    }
  ),

  private = list(
    #' Calculate Cook's distance for mixed model
    .calculate_cooks_d = function(model, leverage, residuals, p) {
      sigma_est <- sigma(model)
      n <- length(residuals)

      # Standard Cook's D formula adapted for mixed models
      # This is an approximation since true Cook's D for LMM is complex
      cooks_d <- (residuals^2 / (p * sigma_est^2)) *
                 (leverage / (1 - leverage)^2)

      # Handle edge cases
      cooks_d[leverage >= 1] <- NA
      cooks_d[is.infinite(cooks_d)] <- NA

      cooks_d
    }
  )
)

# =============================================================================
# CONVENIENCE FUNCTIONS
# =============================================================================

#' Calculate Observation-Level Influence
#'
#' @param model Fitted lme4 model
#' @return InfluenceDiagnosticsResult object
#' @export
calculate_influence <- function(model) {
  diagnostics <- InfluenceDiagnostics$new()
  diagnostics$calculate_observation_influence(model)
}

#' Calculate Participant-Level Influence
#'
#' @param model Fitted lme4 model
#' @param data Data frame
#' @param ... Additional arguments
#' @return ParticipantInfluenceResult object
#' @export
calculate_participant_influence <- function(model, data, ...) {
  diagnostics <- InfluenceDiagnostics$new()
  diagnostics$calculate_participant_influence(model, data, ...)
}
