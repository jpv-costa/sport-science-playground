# R/models/meta_regression.R
# Model: Multi-Level Meta-Regression
#
# SOLID Principles Applied:
# - SRP: Wraps meta-regression fitting and inference (actor: statistician)
# - OCP: Open for extension with different model specifications
# - ISP: Minimal interface for model fitting and results extraction

box::use(
 R6[R6Class]
)

#' Meta-Regression Result
#'
#' Value object containing meta-regression results.
#'
#' @export
MetaRegressionResult <- R6Class(

  classname = "MetaRegressionResult",
  cloneable = FALSE,

  public = list(

    #' @description Create result object
    #' @param model The fitted rma.mv model object
    #' @param outcome Outcome type (Strength/Hypertrophy)
    #' @param n_effects Number of effect sizes
    #' @param n_studies Number of studies
    initialize = function(model, outcome, n_effects, n_studies) {
      private$.model <- model
      private$.outcome <- outcome
      private$.n_effects <- n_effects
      private$.n_studies <- n_studies
    },

    # =========================================================================
    # Query Methods
    # =========================================================================

    #' @description Get coefficient for a predictor
    #' @param predictor Name of predictor variable
    #' @return Named list with estimate, se, pvalue, ci_lower, ci_upper
    get_coefficient = function(predictor) {
      model <- private$.model
      # Handle different ways the predictor might be named
      row_names <- rownames(model$b)
      idx <- which(row_names == predictor)

      # Try partial match if exact match fails
      if (length(idx) == 0) {
        idx <- grep(paste0("^", predictor), row_names)
      }

      if (length(idx) == 0) {
        return(list(
          estimate = NA_real_,
          standard_error = NA_real_,
          p_value = NA_real_,
          ci_lower = NA_real_,
          ci_upper = NA_real_
        ))
      }

      # Take first match if multiple
      idx <- idx[1]

      list(
        estimate = model$b[idx, 1],
        standard_error = model$se[idx],
        p_value = model$pval[idx],
        ci_lower = model$ci.lb[idx],
        ci_upper = model$ci.ub[idx]
      )
    },

    #' @description Check if RIR effect is statistically significant
    #' @param alpha Significance level (default: 0.05)
    #' @return Logical
    is_rir_significant = function(alpha = 0.05) {
      rir_coef <- self$get_coefficient("avg.rir")
      rir_coef$p_value < alpha
    },

    #' @description Get all coefficients as data frame
    #' @return Data frame with all model coefficients
    get_coefficients_table = function() {
      model <- private$.model

      data.frame(
        predictor = rownames(model$b),
        estimate = model$b[, 1],
        se = model$se,
        z_value = model$zval,
        p_value = model$pval,
        ci_lower = model$ci.lb,
        ci_upper = model$ci.ub,
        row.names = NULL
      )
    },

    #' @description Get variance components
    #' @return Named list with sigma2 values
    get_variance_components = function() {
      model <- private$.model
      list(
        sigma2_study = model$sigma2[1],
        sigma2_group = model$sigma2[2],
        sigma2_obs = model$sigma2[3],
        total_heterogeneity = sum(model$sigma2)
      )
    },

    #' @description Get model fit statistics
    #' @return Named list with fit indices
    get_fit_statistics = function() {
      model <- private$.model
      list(
        log_likelihood = model$fit.stats$ML[1],
        aic = model$fit.stats$ML[3],
        bic = model$fit.stats$ML[4],
        aicc = model$fit.stats$ML[5]
      )
    },

    #' @description Create summary for reporting
    #' @return Named list with key results
    summarize = function() {
      rir_coef <- tryCatch(
        self$get_coefficient("avg.rir"),
        error = function(e) NULL
      )

      list(
        outcome = private$.outcome,
        n_effects = private$.n_effects,
        n_studies = private$.n_studies,
        rir_effect = if (!is.null(rir_coef)) rir_coef$estimate else NA,
        rir_p_value = if (!is.null(rir_coef)) rir_coef$p_value else NA,
        rir_significant = if (!is.null(rir_coef)) rir_coef$p_value < 0.05 else NA,
        variance_components = self$get_variance_components()
      )
    }
  ),

  active = list(
    model = function() private$.model,
    outcome = function() private$.outcome,
    n_effects = function() private$.n_effects,
    n_studies = function() private$.n_studies
  ),

  private = list(
    .model = NULL,
    .outcome = NULL,
    .n_effects = NULL,
    .n_studies = NULL
  )
)


#' Multi-Level Meta-Regression Model
#'
#' Fits multi-level meta-regression models for RIR dose-response analysis.
#' Implements the Pelland et al. methodology.
#'
#' @description
#' Service class for fitting and interpreting multi-level meta-regression
#' models with nested random effects structure.
#'
#' @export
MetaRegressionModel <- R6Class(

  classname = "MetaRegressionModel",

  public = list(

    # =========================================================================
    # Constructor
    # =========================================================================

    #' @description Create model specification
    #' @param moderators Character vector of moderator variable names
    #' @param random_effects List of random effects formulas
    #' @param method Estimation method (default: "REML")
    #' @param test Test type (default: "t")
    initialize = function(moderators = c("avg.rir", "load.set", "set.rep.equated",
                                          "weeks", "train.status"),
                          random_effects = list(~1|study, ~1|group, ~1|obs),
                          method = "REML",
                          test = "t") {

      private$.moderators <- moderators
      private$.random_effects <- random_effects
      private$.method <- method
      private$.test <- test
    },

    # =========================================================================
    # Public Methods
    # =========================================================================

    #' @description Fit meta-regression model
    #' @param data Data frame with yi, vi, and moderator columns
    #' @param outcome Outcome type for labeling (Strength/Hypertrophy)
    #' @return MetaRegressionResult object
    fit = function(data, outcome = "Unknown") {
      if (!requireNamespace("metafor", quietly = TRUE)) {
        stop("Package 'metafor' is required for meta-regression")
      }

      private$.validate_data(data)

      # Build formula
      formula <- private$.build_formula()

      # Count studies
      n_studies <- length(unique(data$study))
      n_effects <- nrow(data)

      # Fit model
      model <- metafor::rma.mv(
        yi = yi,
        V = vi,
        mods = formula,
        data = data,
        random = private$.random_effects,
        method = private$.method,
        test = private$.test,
        dfs = "contain"
      )

      MetaRegressionResult$new(
        model = model,
        outcome = outcome,
        n_effects = n_effects,
        n_studies = n_studies
      )
    },

    #' @description Fit models for both outcomes
    #' @param strength_data Data for strength outcome
    #' @param hypertrophy_data Data for hypertrophy outcome
    #' @return Named list with strength and hypertrophy results
    fit_both_outcomes = function(strength_data, hypertrophy_data) {
      list(
        strength = self$fit(strength_data, "Strength"),
        hypertrophy = self$fit(hypertrophy_data, "Hypertrophy")
      )
    }
  ),

  private = list(

    .moderators = NULL,
    .random_effects = NULL,
    .method = NULL,
    .test = NULL,

    .validate_data = function(data) {
      required_cols <- c("yi", "vi", "study", "group", "obs", private$.moderators)
      missing_cols <- setdiff(required_cols, names(data))

      if (length(missing_cols) > 0) {
        stop(paste("Missing required columns:", paste(missing_cols, collapse = ", ")))
      }

      if (any(is.na(data$yi)) || any(is.na(data$vi))) {
        stop("Effect sizes (yi) and variances (vi) cannot contain NA values")
      }
    },

    .build_formula = function() {
      stats::as.formula(paste("~", paste(private$.moderators, collapse = " + ")))
    }
  )
)
