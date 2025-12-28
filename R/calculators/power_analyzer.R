# R/calculators/power_analyzer.R
# Power Analysis for Mixed Effects Models
#
# =============================================================================
# EDUCATIONAL OVERVIEW
# =============================================================================
#
# What is Statistical Power?
# --------------------------
# Power is the probability of detecting an effect IF IT EXISTS.
# - Power = 1 - β (where β is Type II error rate)
# - Conventional target: 80% power
#
# What Affects Power in Mixed Models?
# -----------------------------------
# 1. Effect size (larger effects are easier to detect)
# 2. Number of participants (more clusters = more power)
# 3. Observations per participant (more obs = more precise estimates)
# 4. ICC (higher ICC = less independent information)
# 5. Variance explained by random effects
#
# Important Caveats:
# ------------------
# - **Retrospective power** (post-hoc) is controversial because it's
#   mathematically related to the p-value - low power for non-significant
#   results is circular
# - More useful: "What effect size COULD we have detected?"
# - Best use: Planning FUTURE studies
#
# =============================================================================

box::use(
  R6[R6Class],
  stats[pnorm, qnorm, pt, qt, var, sd, coef, vcov, power.t.test, sigma]
)

#' Power Analysis Result
#'
#' Immutable value object containing power analysis results
#'
#' @export
PowerAnalysisResult <- R6Class(
  classname = "PowerAnalysisResult",
  cloneable = FALSE,

  public = list(
    #' @field observed_effect Observed effect size from the model
    observed_effect = NULL,
    #' @field standardized_effect Cohen's d or similar
    standardized_effect = NULL,
    #' @field observed_power Post-hoc power (use with caution!)
    observed_power = NULL,
    #' @field detectable_effect Minimum detectable effect at 80% power
    detectable_effect = NULL,
    #' @field n_participants Number of participants in study
    n_participants = NULL,
    #' @field n_observations Total observations
    n_observations = NULL,
    #' @field effective_n Effective sample size after accounting for clustering
    effective_n = NULL,
    #' @field icc Intraclass correlation
    icc = NULL,

    #' @description Create a new power analysis result
    initialize = function(observed_effect, standardized_effect, observed_power,
                          detectable_effect, n_participants, n_observations,
                          effective_n, icc) {
      self$observed_effect <- observed_effect
      self$standardized_effect <- standardized_effect
      self$observed_power <- observed_power
      self$detectable_effect <- detectable_effect
      self$n_participants <- n_participants
      self$n_observations <- n_observations
      self$effective_n <- effective_n
      self$icc <- icc
    },

    #' @description Interpret power adequacy
    interpret = function() {
      if (is.na(self$observed_power)) {
        return("Power calculation not available")
      }

      power_level <- if (self$observed_power >= 0.90) {
        "excellent"
      } else if (self$observed_power >= 0.80) {
        "adequate"
      } else if (self$observed_power >= 0.60) {
        "moderate"
      } else {
        "low"
      }

      paste0(
        "Observed power: ", round(self$observed_power * 100, 1), "% (",
        power_level, ")\n",
        "With this design, we could reliably detect effects ≥ ",
        round(self$detectable_effect, 4), " (at 80% power)"
      )
    },

    #' @description Convert to list for serialization
    to_list = function() {
      list(
        observed_effect = self$observed_effect,
        standardized_effect = self$standardized_effect,
        observed_power = self$observed_power,
        detectable_effect = self$detectable_effect,
        n_participants = self$n_participants,
        n_observations = self$n_observations,
        effective_n = self$effective_n,
        icc = self$icc,
        interpretation = self$interpret()
      )
    }
  )
)

#' Sample Size Result
#'
#' Result from sample size calculation
#'
#' @export
SampleSizeResult <- R6Class(
  classname = "SampleSizeResult",
  cloneable = FALSE,

  public = list(
    #' @field required_n_participants Required number of participants
    required_n_participants = NULL,
    #' @field obs_per_participant Assumed observations per participant
    obs_per_participant = NULL,
    #' @field target_effect Target effect size
    target_effect = NULL,
    #' @field target_power Target power level
    target_power = NULL,
    #' @field alpha Significance level
    alpha = NULL,
    #' @field assumed_icc Assumed ICC for calculation
    assumed_icc = NULL,

    #' @description Create a new sample size result
    initialize = function(required_n_participants, obs_per_participant,
                          target_effect, target_power, alpha, assumed_icc) {
      self$required_n_participants <- required_n_participants
      self$obs_per_participant <- obs_per_participant
      self$target_effect <- target_effect
      self$target_power <- target_power
      self$alpha <- alpha
      self$assumed_icc <- assumed_icc
    },

    #' @description Convert to list for serialization
    to_list = function() {
      list(
        required_n_participants = self$required_n_participants,
        total_observations = self$required_n_participants * self$obs_per_participant,
        target_effect = self$target_effect,
        target_power = self$target_power,
        alpha = self$alpha,
        assumed_icc = self$assumed_icc
      )
    }
  )
)

#' Power Analyzer
#'
#' R6 class for power analysis and sample size calculations
#' for mixed effects models.
#'
#' @section Caveats:
#' - Retrospective power is controversial - interpret with caution
#' - Better question: "What effect could we have detected?"
#' - These calculations are approximations for mixed models
#'
#' @export
PowerAnalyzer <- R6Class(
  classname = "PowerAnalyzer",
  cloneable = FALSE,

  public = list(
    #' @description Create a new PowerAnalyzer instance
    initialize = function() {
      # No initialization needed
    },

    # =========================================================================
    # POWER CALCULATIONS
    # =========================================================================

    #' @description Calculate Power Analysis for Mixed Model
    #'
    #' Computes observed power and minimum detectable effect size.
    #'
    #' @section Important Note:
    #' Retrospective power is mathematically related to the p-value.
    #' For non-significant results, power will necessarily be low.
    #' The more useful metric is the "minimum detectable effect."
    #'
    #' @param model Fitted lme4 model
    #' @param data Data frame used for fitting
    #' @param id_col Name of participant ID column
    #' @param coef_name Name of coefficient to analyze (default: main predictor)
    #' @param alpha Significance level (default 0.05)
    #' @param target_power Target power for detectable effect (default 0.80)
    #' @return PowerAnalysisResult object
    calculate_power = function(model, data, id_col = "id", coef_name = NULL,
                               alpha = 0.05, target_power = 0.80) {
      if (!requireNamespace("lme4", quietly = TRUE)) {
        stop("Package 'lme4' required for calculate_power()")
      }

      # Get fixed effects and SEs
      fixed_effects <- lme4::fixef(model)
      vcov_mat <- as.matrix(vcov(model))

      if (is.null(coef_name)) {
        coef_name <- names(fixed_effects)[2]  # Usually main predictor
      }

      observed_effect <- fixed_effects[coef_name]
      se <- sqrt(vcov_mat[coef_name, coef_name])

      # Sample sizes
      n_observations <- nrow(data)
      n_participants <- length(unique(data[[id_col]]))
      obs_per_participant <- n_observations / n_participants

      # Calculate ICC
      vc <- lme4::VarCorr(model)
      sigma_between <- sqrt(sum(sapply(vc, function(x) sum(diag(x)))))
      sigma_within <- sigma(model)
      icc <- sigma_between^2 / (sigma_between^2 + sigma_within^2)

      # Design effect and effective sample size
      design_effect <- 1 + (obs_per_participant - 1) * icc
      effective_n <- n_observations / design_effect

      # Standardized effect size (approximate Cohen's d)
      total_sd <- sqrt(sigma_between^2 + sigma_within^2)
      standardized_effect <- observed_effect / total_sd

      # Calculate observed power using z-test approximation
      z_crit <- qnorm(1 - alpha / 2)
      ncp <- abs(observed_effect) / se  # Non-centrality parameter
      observed_power <- pnorm(ncp - z_crit) + pnorm(-ncp - z_crit)

      # Calculate minimum detectable effect at target power
      # Solve for effect size: power = Φ(|effect|/se - z_crit)
      z_power <- qnorm(target_power)
      detectable_effect <- (z_crit + z_power) * se

      PowerAnalysisResult$new(
        observed_effect = observed_effect,
        standardized_effect = standardized_effect,
        observed_power = observed_power,
        detectable_effect = detectable_effect,
        n_participants = n_participants,
        n_observations = n_observations,
        effective_n = effective_n,
        icc = icc
      )
    },

    # =========================================================================
    # SAMPLE SIZE CALCULATIONS
    # =========================================================================

    #' @description Calculate Required Sample Size
    #'
    #' Estimates the number of participants needed to detect a
    #' specified effect size with target power.
    #'
    #' @param target_effect Target effect size (in original units)
    #' @param residual_sd Residual standard deviation
    #' @param icc Intraclass correlation coefficient
    #' @param obs_per_participant Expected observations per participant
    #' @param power Target power (default 0.80)
    #' @param alpha Significance level (default 0.05)
    #' @return SampleSizeResult object
    calculate_required_n = function(target_effect, residual_sd, icc,
                                    obs_per_participant, power = 0.80,
                                    alpha = 0.05) {
      # Design effect
      design_effect <- 1 + (obs_per_participant - 1) * icc

      # Standardized effect
      d <- abs(target_effect) / residual_sd

      # Required n for simple design (no clustering)
      z_alpha <- qnorm(1 - alpha / 2)
      z_power <- qnorm(power)
      n_simple <- 2 * ((z_alpha + z_power) / d)^2

      # Adjust for clustering
      n_per_group <- ceiling(n_simple * design_effect / obs_per_participant)

      SampleSizeResult$new(
        required_n_participants = n_per_group,
        obs_per_participant = obs_per_participant,
        target_effect = target_effect,
        target_power = power,
        alpha = alpha,
        assumed_icc = icc
      )
    },

    #' @description Calculate Sample Size from Existing Model
    #'
    #' Uses variance estimates from a pilot/existing study to
    #' calculate required sample size for a future study.
    #'
    #' @param model Fitted lme4 model (used for variance estimates)
    #' @param data Data frame used for fitting
    #' @param id_col Name of participant ID column
    #' @param target_effect Target effect size (in original units)
    #' @param power Target power (default 0.80)
    #' @param alpha Significance level (default 0.05)
    #' @return SampleSizeResult object
    calculate_required_n_from_model = function(model, data, id_col = "id",
                                                target_effect, power = 0.80,
                                                alpha = 0.05) {
      if (!requireNamespace("lme4", quietly = TRUE)) {
        stop("Package 'lme4' required for calculate_required_n_from_model()")
      }

      # Extract variance components
      vc <- lme4::VarCorr(model)
      sigma_between <- sqrt(sum(sapply(vc, function(x) sum(diag(x)))))
      sigma_within <- sigma(model)

      # Calculate ICC
      icc <- sigma_between^2 / (sigma_between^2 + sigma_within^2)

      # Current observations per participant
      n_observations <- nrow(data)
      n_participants <- length(unique(data[[id_col]]))
      obs_per_participant <- n_observations / n_participants

      self$calculate_required_n(
        target_effect = target_effect,
        residual_sd = sigma_within,
        icc = icc,
        obs_per_participant = obs_per_participant,
        power = power,
        alpha = alpha
      )
    },

    # =========================================================================
    # POWER CURVES
    # =========================================================================

    #' @description Generate Power Curve
    #'
    #' Shows how power varies with sample size for a given effect.
    #'
    #' @param effect_size Effect size to detect
    #' @param se Standard error at reference sample size
    #' @param n_range Vector of sample sizes to evaluate
    #' @param n_reference Reference sample size (for scaling SE)
    #' @param alpha Significance level
    #' @return Data frame with power at each sample size
    generate_power_curve = function(effect_size, se, n_range,
                                    n_reference, alpha = 0.05) {
      z_crit <- qnorm(1 - alpha / 2)

      # SE scales with sqrt(n)
      se_scaled <- se * sqrt(n_reference / n_range)

      # Non-centrality parameter
      ncp <- abs(effect_size) / se_scaled

      # Power
      power <- pnorm(ncp - z_crit) + pnorm(-ncp - z_crit)

      data.frame(
        n = n_range,
        se = se_scaled,
        power = power
      )
    },

    # =========================================================================
    # SUMMARY AND REPORTING
    # =========================================================================

    #' @description Create Power Analysis Summary
    #'
    #' Generates a comprehensive summary table for reporting.
    #'
    #' @param model Fitted lme4 model
    #' @param data Data frame used for fitting
    #' @param id_col Name of participant ID column
    #' @param coef_name Name of coefficient (optional)
    #' @return List with summary tables and interpretation
    create_power_summary = function(model, data, id_col = "id", coef_name = NULL) {
      # Calculate power analysis
      power_result <- self$calculate_power(model, data, id_col, coef_name)

      # Sample size for replication
      required_n <- self$calculate_required_n_from_model(
        model, data, id_col,
        target_effect = power_result$observed_effect,
        power = 0.80
      )

      # Sample size for smaller effect
      small_effect <- power_result$observed_effect * 0.5
      required_n_small <- self$calculate_required_n_from_model(
        model, data, id_col,
        target_effect = small_effect,
        power = 0.80
      )

      list(
        current_study = data.frame(
          metric = c("Participants", "Total observations", "Effective N",
                     "ICC", "Observed effect", "Standardized effect (d)",
                     "Standard error", "Observed power", "Minimum detectable effect"),
          value = c(
            power_result$n_participants,
            power_result$n_observations,
            round(power_result$effective_n, 1),
            round(power_result$icc, 3),
            round(power_result$observed_effect, 4),
            round(power_result$standardized_effect, 2),
            round(sqrt(as.matrix(vcov(model))[2, 2]), 4),
            round(power_result$observed_power, 3),
            round(power_result$detectable_effect, 4)
          ),
          stringsAsFactors = FALSE
        ),
        future_study = data.frame(
          scenario = c("Replicate current effect", "Detect half the effect"),
          target_effect = c(round(power_result$observed_effect, 4),
                            round(small_effect, 4)),
          required_n = c(required_n$required_n_participants,
                         required_n_small$required_n_participants),
          power = c(0.80, 0.80),
          stringsAsFactors = FALSE
        ),
        interpretation = power_result$interpret(),
        caveats = c(
          "Retrospective power should be interpreted cautiously.",
          "Post-hoc power is mathematically related to p-values.",
          "More useful: What effect size could we have detected?",
          "For planning, use prospective power calculations."
        )
      )
    }
  )
)

# =============================================================================
# CONVENIENCE FUNCTIONS
# =============================================================================

#' Calculate Power Analysis
#'
#' @param model Fitted lme4 model
#' @param data Data frame
#' @param ... Additional arguments
#' @return PowerAnalysisResult object
#' @export
calculate_power <- function(model, data, ...) {
  analyzer <- PowerAnalyzer$new()
  analyzer$calculate_power(model, data, ...)
}

#' Calculate Required Sample Size
#'
#' @param model Fitted lme4 model
#' @param data Data frame
#' @param target_effect Target effect size
#' @param ... Additional arguments
#' @return SampleSizeResult object
#' @export
calculate_required_n <- function(model, data, target_effect, ...) {
  analyzer <- PowerAnalyzer$new()
  analyzer$calculate_required_n_from_model(model, data, target_effect = target_effect, ...)
}
