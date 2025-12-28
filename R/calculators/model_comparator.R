# R/calculators/model_comparator.R
# Service: Model Comparison and Selection for Mixed Effects Models
#
# =============================================================================
# PRINCIPLES APPLIED (from CLAUDE.md)
# =============================================================================
#
# NAMING PRINCIPLES:
# - Understandability: Domain terms (bayes_factor, aic_comparison, nested_test)
# - Consistency: All public methods follow verb_noun pattern
# - Distinguishability: Clear method names (compare_aic vs compare_bayes)
# - Conciseness: BF = Bayes Factor, LRT = Likelihood Ratio Test
#
# SOLID PRINCIPLES:
# - SRP: Single responsibility - only model comparison logic
# - OCP: Open for extension via new comparison methods
# - DIP: Depends on model objects, not specific implementations
#
# CUPID PRINCIPLES:
# - Composable: Each comparison method is independent
# - Unix: Each method does one comparison type
# - Predictable: Same models -> same comparison results
# - Idiomatic: Follows R/lme4 conventions
# - Domain-based: Names reflect statistical concepts
#
# SCIENTIFIC VALIDITY:
# - Bayes factors use BIC approximation (Wagenmakers 2007)
# - LRT uses chi-square distribution for nested models
# - Effect sizes use standardized and prediction difference methods
# =============================================================================

box::use(
  R6[R6Class],
  stats[AIC, BIC, anova, logLik, predict, sd, pchisq]
)

#' Model Comparison Result
#'
#' Value object containing comparison metrics between two models.
#'
#' @export
ComparisonResult <- R6Class(
  classname = "ComparisonResult",
  cloneable = FALSE,

 public = list(
    #' @field model1_name Name of first model
    model1_name = NULL,

    #' @field model2_name Name of second model
    model2_name = NULL,

    #' @field aic_diff AIC difference (model2 - model1)
    aic_diff = NULL,

    #' @field bic_diff BIC difference (model2 - model1)
    bic_diff = NULL,

    #' @field bayes_factor Approximate Bayes factor favoring model1
    bayes_factor = NULL,

    #' @field lrt_chisq Likelihood ratio test chi-square
    lrt_chisq = NULL,

    #' @field lrt_df Degrees of freedom for LRT
    lrt_df = NULL,

    #' @field lrt_pvalue P-value for LRT
    lrt_pvalue = NULL,

    #' @field preferred_model Which model is preferred
    preferred_model = NULL,

    #' @field evidence_strength Strength of evidence (weak/moderate/strong)
    evidence_strength = NULL,

    #' @description Create comparison result
    initialize = function(model1_name, model2_name, aic_diff, bic_diff,
                          bayes_factor = NULL, lrt_chisq = NULL,
                          lrt_df = NULL, lrt_pvalue = NULL) {
      self$model1_name <- model1_name
      self$model2_name <- model2_name
      self$aic_diff <- aic_diff
      self$bic_diff <- bic_diff
      self$bayes_factor <- bayes_factor
      self$lrt_chisq <- lrt_chisq
      self$lrt_df <- lrt_df
      self$lrt_pvalue <- lrt_pvalue

      # Determine preferred model and evidence strength
      private$.determine_preference()
    },

    #' @description Get interpretation text
    #' @return Character string with interpretation
    interpret = function() {
      bf_text <- ""
      if (!is.null(self$bayes_factor)) {
        bf_text <- sprintf(
          "Bayes factor = %.2f (%s evidence for %s)",
          self$bayes_factor,
          self$evidence_strength,
          self$preferred_model
        )
      }

      aic_text <- sprintf("ΔAIC = %.1f", self$aic_diff)
      bic_text <- sprintf("ΔBIC = %.1f", self$bic_diff)

      lrt_text <- ""
      if (!is.null(self$lrt_pvalue)) {
        lrt_text <- sprintf(
          "LRT: χ²(%.0f) = %.2f, p = %s",
          self$lrt_df, self$lrt_chisq,
          ifelse(self$lrt_pvalue < 0.001, "<0.001", sprintf("%.3f", self$lrt_pvalue))
        )
      }

      paste(c(bf_text, aic_text, bic_text, lrt_text), collapse = " | ")
    },

    #' @description Convert to list
    #' @return List representation
    to_list = function() {
      list(
        model1 = self$model1_name,
        model2 = self$model2_name,
        aic_diff = self$aic_diff,
        bic_diff = self$bic_diff,
        bayes_factor = self$bayes_factor,
        lrt_chisq = self$lrt_chisq,
        lrt_df = self$lrt_df,
        lrt_pvalue = self$lrt_pvalue,
        preferred_model = self$preferred_model,
        evidence_strength = self$evidence_strength
      )
    }
  ),

  private = list(
    .determine_preference = function() {
      # Use Bayes factor if available, otherwise AIC
      if (!is.null(self$bayes_factor)) {
        if (self$bayes_factor > 1) {
          self$preferred_model <- self$model1_name
          bf <- self$bayes_factor
        } else {
          self$preferred_model <- self$model2_name
          bf <- 1 / self$bayes_factor
        }

        # Interpret Bayes factor strength (Jeffreys scale)
        self$evidence_strength <- if (bf < 3) {
          "weak/anecdotal"
        } else if (bf < 10) {
          "moderate"
        } else if (bf < 30) {
          "strong"
        } else if (bf < 100) {
          "very strong"
        } else {
          "decisive"
        }
      } else {
        # Fall back to AIC
        self$preferred_model <- if (self$aic_diff > 0) {
          self$model1_name
        } else {
          self$model2_name
        }
        self$evidence_strength <- if (abs(self$aic_diff) < 2) {
          "negligible"
        } else if (abs(self$aic_diff) < 7) {
          "moderate"
        } else {
          "strong"
        }
      }
    }
  )
)

#' Model Comparator
#'
#' R6 class for comparing statistical models using multiple metrics.
#' Implements Bayes factors, AIC/BIC, and likelihood ratio tests.
#'
#' @section Statistical Background:
#' - Bayes Factor approximation uses BIC (Wagenmakers, 2007)
#' - BF = exp((BIC_0 - BIC_1) / 2)
#' - Interpretation follows Jeffreys (1961) scale
#'
#' @export
ModelComparator <- R6Class(
  classname = "ModelComparator",
  cloneable = FALSE,

  public = list(
    #' @description Create a new ModelComparator
    initialize = function() {
      # No initialization needed
    },

    #' @description Compare two models using all available metrics
    #' @param model1 First model (typically simpler/null)
    #' @param model2 Second model (typically more complex/alternative)
    #' @param model1_name Name for first model
    #' @param model2_name Name for second model
    #' @param nested Are models nested? (for LRT)
    #' @return ComparisonResult object
    compare = function(model1, model2,
                       model1_name = "Model 1",
                       model2_name = "Model 2",
                       nested = TRUE) {
      # Calculate information criteria differences
      aic1 <- AIC(model1)
      aic2 <- AIC(model2)
      bic1 <- BIC(model1)
      bic2 <- BIC(model2)

      aic_diff <- aic2 - aic1  # Positive means model1 is better
      bic_diff <- bic2 - bic1

      # Calculate Bayes factor approximation from BIC
      # BF_10 = exp((BIC_0 - BIC_1) / 2) where 0 is null, 1 is alternative
      bayes_factor <- exp(bic_diff / 2)

      # Likelihood ratio test for nested models
      lrt_chisq <- NULL
      lrt_df <- NULL
      lrt_pvalue <- NULL

      if (nested) {
        tryCatch({
          lrt_result <- anova(model1, model2)
          if ("Chisq" %in% names(lrt_result)) {
            lrt_chisq <- lrt_result$Chisq[2]
            lrt_df <- lrt_result$Df[2]
            lrt_pvalue <- lrt_result$`Pr(>Chisq)`[2]
          }
        }, error = function(e) {
          # LRT not available for these models
        })
      }

      ComparisonResult$new(
        model1_name = model1_name,
        model2_name = model2_name,
        aic_diff = aic_diff,
        bic_diff = bic_diff,
        bayes_factor = bayes_factor,
        lrt_chisq = lrt_chisq,
        lrt_df = lrt_df,
        lrt_pvalue = lrt_pvalue
      )
    },

    #' @description Compare multiple models
    #' @param models Named list of models
    #' @param reference_model Name of reference model (default: first)
    #' @return Data frame with comparison metrics
    compare_multiple = function(models, reference_model = NULL) {
      if (is.null(reference_model)) {
        reference_model <- names(models)[1]
      }

      ref_mod <- models[[reference_model]]

      results <- lapply(names(models), function(name) {
        if (name == reference_model) {
          return(data.frame(
            model = name,
            aic = AIC(models[[name]]),
            bic = BIC(models[[name]]),
            delta_aic = 0,
            delta_bic = 0,
            bayes_factor = 1,
            evidence = "reference",
            stringsAsFactors = FALSE
          ))
        }

        comp <- self$compare(
          ref_mod, models[[name]],
          reference_model, name,
          nested = FALSE
        )

        data.frame(
          model = name,
          aic = AIC(models[[name]]),
          bic = BIC(models[[name]]),
          delta_aic = -comp$aic_diff,  # Relative to reference
          delta_bic = -comp$bic_diff,
          bayes_factor = comp$bayes_factor,
          evidence = comp$evidence_strength,
          stringsAsFactors = FALSE
        )
      })

      do.call(rbind, results)
    },

    #' @description Calculate prediction difference effect size
    #' @param model Fitted model
    #' @param data Original data
    #' @param predictor Name of predictor variable
    #' @param sd_range Number of SDs above/below mean (default: 1)
    #' @return List with prediction difference metrics
    calculate_prediction_difference = function(model, data, predictor,
                                                sd_range = 1) {
      stopifnot(
        "predictor must be in data" = predictor %in% names(data)
      )

      pred_vals <- data[[predictor]]
      pred_mean <- mean(pred_vals, na.rm = TRUE)
      pred_sd <- sd(pred_vals, na.rm = TRUE)

      # Create prediction data at +/- 1 SD
      low_val <- pred_mean - sd_range * pred_sd
      high_val <- pred_mean + sd_range * pred_sd

      # Create new data for prediction
      new_data_low <- data[1, , drop = FALSE]
      new_data_high <- data[1, , drop = FALSE]

      # Set predictor to low and high values
      new_data_low[[predictor]] <- low_val
      new_data_high[[predictor]] <- high_val

      # Get predictions (population level for LMM)
      pred_low <- predict(model, newdata = new_data_low, re.form = NA)
      pred_high <- predict(model, newdata = new_data_high, re.form = NA)

      prediction_diff <- pred_high - pred_low

      list(
        predictor = predictor,
        low_value = low_val,
        high_value = high_val,
        pred_at_low = pred_low,
        pred_at_high = pred_high,
        prediction_difference = prediction_diff,
        interpretation = sprintf(
          "A %d SD increase in %s is associated with a %.4f unit change in the outcome",
          2 * sd_range, predictor, prediction_diff
        )
      )
    },

    #' @description Generate model comparison summary table
    #' @param models Named list of models
    #' @return Formatted data frame for display
    create_comparison_table = function(models) {
      comparison <- self$compare_multiple(models)

      # Add R² if available (for lmer models)
      comparison$r2_marginal <- NA
      comparison$r2_conditional <- NA

      for (i in seq_len(nrow(comparison))) {
        mod <- models[[comparison$model[i]]]
        if (inherits(mod, "lmerMod")) {
          tryCatch({
            r2 <- MuMIn::r.squaredGLMM(mod)
            comparison$r2_marginal[i] <- r2[1, "R2m"]
            comparison$r2_conditional[i] <- r2[1, "R2c"]
          }, error = function(e) {
            # R² calculation failed
          })
        }
      }

      comparison
    },

    #' @description Interpret Bayes factor with user-friendly text
    #' @param bf Bayes factor value
    #' @return Character string interpretation
    interpret_bayes_factor = function(bf) {
      if (bf < 1/100) {
        "Decisive evidence against"
      } else if (bf < 1/30) {
        "Very strong evidence against"
      } else if (bf < 1/10) {
        "Strong evidence against"
      } else if (bf < 1/3) {
        "Moderate evidence against"
      } else if (bf < 1) {
        "Weak evidence against"
      } else if (bf < 3) {
        "Weak evidence for"
      } else if (bf < 10) {
        "Moderate evidence for"
      } else if (bf < 30) {
        "Strong evidence for"
      } else if (bf < 100) {
        "Very strong evidence for"
      } else {
        "Decisive evidence for"
      }
    }
  )
)

#' Convenience function to compare two models
#'
#' @param model1 First model
#' @param model2 Second model
#' @param ... Additional arguments passed to ModelComparator$compare
#' @return ComparisonResult object
#' @export
compare_models <- function(model1, model2, ...) {
  comparator <- ModelComparator$new()
  comparator$compare(model1, model2, ...)
}
