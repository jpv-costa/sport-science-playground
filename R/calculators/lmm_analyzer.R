# R/calculators/lmm_analyzer.R
# Service: Linear Mixed Effects Model Analysis
#
# SOLID Principles Applied:
# - SRP: Single responsibility - fit and evaluate LMM models
# - OCP: Open for extension with different random effect structures

box::use(
  R6[R6Class]
)

#' LMM Model Result
#'
#' Value object containing LMM model fit results.
#'
#' @export
LmmModelResult <- R6Class(

  classname = "LmmModelResult",
  cloneable = FALSE,

  public = list(

    #' @description Create LMM model result
    #' @param model The fitted lmer model object
    #' @param model_name Name identifier for the model
    #' @param formula_string String representation of the formula
    #' @param aic AIC value
    #' @param bic BIC value
    #' @param log_likelihood Log-likelihood value
    #' @param r2_marginal Marginal R-squared (fixed effects only)
    #' @param r2_conditional Conditional R-squared (fixed + random)
    #' @param n_observations Number of observations
    #' @param n_groups Number of grouping levels
    initialize = function(model,
                          model_name,
                          formula_string,
                          aic,
                          bic,
                          log_likelihood,
                          r2_marginal,
                          r2_conditional,
                          n_observations,
                          n_groups) {
      private$.model <- model
      private$.model_name <- model_name
      private$.formula_string <- formula_string
      private$.aic <- aic
      private$.bic <- bic
      private$.log_likelihood <- log_likelihood
      private$.r2_marginal <- r2_marginal
      private$.r2_conditional <- r2_conditional
      private$.n_observations <- n_observations
      private$.n_groups <- n_groups
    },

    #' @description Convert to list
    to_list = function() {
      list(
        model_name = private$.model_name,
        formula = private$.formula_string,
        aic = private$.aic,
        bic = private$.bic,
        log_likelihood = private$.log_likelihood,
        r2_marginal = private$.r2_marginal,
        r2_conditional = private$.r2_conditional,
        n_observations = private$.n_observations,
        n_groups = private$.n_groups
      )
    }
  ),

  active = list(
    model = function() private$.model,
    model_name = function() private$.model_name,
    formula_string = function() private$.formula_string,
    aic = function() private$.aic,
    bic = function() private$.bic,
    log_likelihood = function() private$.log_likelihood,
    r2_marginal = function() private$.r2_marginal,
    r2_conditional = function() private$.r2_conditional,
    n_observations = function() private$.n_observations,
    n_groups = function() private$.n_groups
  ),

  private = list(
    .model = NULL,
    .model_name = NULL,
    .formula_string = NULL,
    .aic = NULL,
    .bic = NULL,
    .log_likelihood = NULL,
    .r2_marginal = NULL,
    .r2_conditional = NULL,
    .n_observations = NULL,
    .n_groups = NULL
  )
)


#' Model Comparison Result
#'
#' Value object containing results from comparing multiple LMM models.
#'
#' @export
LmmModelComparison <- R6Class(

  classname = "LmmModelComparison",
  cloneable = FALSE,

  public = list(

    #' @description Create comparison result
    #' @param comparison_table Data frame with model comparison metrics
    #' @param best_model_name Name of the best model by selected criterion
    #' @param best_model The best LmmModelResult object
    #' @param anova_results ANOVA comparison results (if applicable)
    initialize = function(comparison_table, best_model_name, best_model, anova_results = NULL) {
      private$.comparison_table <- comparison_table
      private$.best_model_name <- best_model_name
      private$.best_model <- best_model
      private$.anova_results <- anova_results
    },

    #' @description Convert to list
    to_list = function() {
      list(
        comparison_table = private$.comparison_table,
        best_model_name = private$.best_model_name,
        anova_results = private$.anova_results
      )
    }
  ),

  active = list(
    comparison_table = function() private$.comparison_table,
    best_model_name = function() private$.best_model_name,
    best_model = function() private$.best_model,
    anova_results = function() private$.anova_results
  ),

  private = list(
    .comparison_table = NULL,
    .best_model_name = NULL,
    .best_model = NULL,
    .anova_results = NULL
  )
)


#' Model Diagnostics Result
#'
#' Value object containing model diagnostic test results.
#'
#' @export
LmmDiagnostics <- R6Class(

  classname = "LmmDiagnostics",
  cloneable = FALSE,

  public = list(

    #' @description Create diagnostics result
    #' @param normality_test List with normality test results
    #' @param homoscedasticity_test List with homoscedasticity assessment
    #' @param residuals Vector of residuals
    #' @param fitted Vector of fitted values
    #' @param random_effects Data frame of random effects
    initialize = function(normality_test,
                          homoscedasticity_test,
                          residuals,
                          fitted,
                          random_effects) {
      private$.normality_test <- normality_test
      private$.homoscedasticity_test <- homoscedasticity_test
      private$.residuals <- residuals
      private$.fitted <- fitted
      private$.random_effects <- random_effects
    },

    #' @description Convert to list
    to_list = function() {
      list(
        normality_test = private$.normality_test,
        homoscedasticity_test = private$.homoscedasticity_test,
        n_residuals = length(private$.residuals)
      )
    }
  ),

  active = list(
    normality_test = function() private$.normality_test,
    homoscedasticity_test = function() private$.homoscedasticity_test,
    residuals = function() private$.residuals,
    fitted = function() private$.fitted,
    random_effects = function() private$.random_effects
  ),

  private = list(
    .normality_test = NULL,
    .homoscedasticity_test = NULL,
    .residuals = NULL,
    .fitted = NULL,
    .random_effects = NULL
  )
)


#' LMM Analyzer
#'
#' Fits and evaluates Linear Mixed Effects Models for nested/hierarchical data.
#' Supports random intercepts and slopes, interaction terms, and model diagnostics.
#'
#' @export
LmmAnalyzer <- R6Class(

  classname = "LmmAnalyzer",

  public = list(

    #' @description Create analyzer
    initialize = function() {
      private$.validate_dependencies()
    },

    #' @description Fit LMM with specified formula and random effects
    #' @param data Data frame with observations
    #' @param formula Formula for fixed effects (e.g., mean_velocity ~ rir + load_percentage)
    #' @param random_formula Random effects formula (e.g., ~1 + rir | id)
    #' @param model_name Optional name for the model
    #' @param REML Use REML estimation (TRUE) or ML (FALSE)
    #' @return LmmModelResult object
    fit = function(data, formula, random_formula = ~1 | id, model_name = NULL, REML = TRUE) {
      if (!requireNamespace("lme4", quietly = TRUE)) {
        stop("Package 'lme4' is required for LMM fitting")
      }

      # Build full formula with random effects
      full_formula <- private$.build_full_formula(formula, random_formula)
      formula_string <- deparse(full_formula, width.cutoff = 500)

      if (is.null(model_name)) {
        model_name <- formula_string
      }

      # Fit model using lmerTest for p-values
      model <- lme4::lmer(full_formula, data = data, REML = REML)

      # Calculate R-squared values
      r2_values <- private$.calculate_r_squared(model)

      LmmModelResult$new(
        model = model,
        model_name = model_name,
        formula_string = formula_string,
        aic = stats::AIC(model),
        bic = stats::BIC(model),
        log_likelihood = as.numeric(stats::logLik(model)),
        r2_marginal = r2_values$marginal,
        r2_conditional = r2_values$conditional,
        n_observations = nrow(data),
        n_groups = length(unique(lme4::getME(model, "flist")[[1]]))
      )
    },

    #' @description Fit LMM with interaction terms
    #' @param data Data frame with observations
    #' @param base_formula Base formula without interactions
    #' @param interactions Character vector of interaction terms (e.g., c("rir:load_percentage"))
    #' @param random_formula Random effects formula
    #' @param model_name Optional name for the model
    #' @return LmmModelResult object
    fit_with_interactions = function(data,
                                     base_formula = mean_velocity ~ rir + load_percentage + day,
                                     interactions = c("rir:load_percentage"),
                                     random_formula = ~1 | id,
                                     model_name = NULL) {
      # Parse base formula
      response <- as.character(base_formula)[2]
      predictors <- as.character(base_formula)[3]

      # Add interaction terms
      full_predictors <- paste(c(predictors, interactions), collapse = " + ")
      full_formula <- stats::as.formula(paste(response, "~", full_predictors))

      if (is.null(model_name)) {
        model_name <- paste("Interactions:", paste(interactions, collapse = ", "))
      }

      self$fit(data, full_formula, random_formula, model_name)
    },

    #' @description Compare multiple models
    #' @param model_results Named list of LmmModelResult objects
    #' @param criterion Comparison criterion ("AIC", "BIC", or "logLik")
    #' @return LmmModelComparison object
    compare_models = function(model_results, criterion = "AIC") {
      # Build comparison table
      comparison_df <- data.frame(
        model = character(),
        AIC = numeric(),
        BIC = numeric(),
        logLik = numeric(),
        R2_marginal = numeric(),
        R2_conditional = numeric(),
        n_obs = integer(),
        n_groups = integer(),
        delta_AIC = numeric(),
        delta_BIC = numeric(),
        bayes_factor_approx = numeric(),
        stringsAsFactors = FALSE
      )

      for (name in names(model_results)) {
        result <- model_results[[name]]
        comparison_df <- rbind(comparison_df, data.frame(
          model = name,
          AIC = result$aic,
          BIC = result$bic,
          logLik = result$log_likelihood,
          R2_marginal = result$r2_marginal,
          R2_conditional = result$r2_conditional,
          n_obs = result$n_observations,
          n_groups = result$n_groups,
          delta_AIC = NA_real_,
          delta_BIC = NA_real_,
          bayes_factor_approx = NA_real_,
          stringsAsFactors = FALSE
        ))
      }

      # Calculate delta values relative to best model
      min_aic <- min(comparison_df$AIC, na.rm = TRUE)
      min_bic <- min(comparison_df$BIC, na.rm = TRUE)
      comparison_df$delta_AIC <- comparison_df$AIC - min_aic
      comparison_df$delta_BIC <- comparison_df$BIC - min_bic

      # Approximate Bayes Factor using BIC (Raftery 1995)
      # BF ≈ exp(-0.5 * delta_BIC)
      comparison_df$bayes_factor_approx <- exp(-0.5 * comparison_df$delta_BIC)

      # Determine best model
      if (criterion == "AIC") {
        best_idx <- which.min(comparison_df$AIC)
      } else if (criterion == "BIC") {
        best_idx <- which.min(comparison_df$BIC)
      } else {
        best_idx <- which.max(comparison_df$logLik)
      }

      best_name <- comparison_df$model[best_idx]

      # Perform ANOVA if models are nested
      anova_results <- NULL
      if (length(model_results) >= 2) {
        anova_results <- tryCatch({
          models_list <- lapply(model_results, function(r) r$model)
          do.call(stats::anova, models_list)
        }, error = function(e) NULL)
      }

      LmmModelComparison$new(
        comparison_table = comparison_df,
        best_model_name = best_name,
        best_model = model_results[[best_name]],
        anova_results = anova_results
      )
    },

    #' @description Test if a variable significantly improves model fit
    #' @param data Data frame with observations
    #' @param base_formula Formula without the variable
    #' @param full_formula Formula with the variable
    #' @param random_formula Random effects formula
    #' @return List with test results (LRT, AIC difference, BF approximation)
    test_variable_importance = function(data,
                                        base_formula,
                                        full_formula,
                                        random_formula = ~1 | id) {
      # Fit both models with ML (required for LRT comparison)
      model_without <- self$fit(data, base_formula, random_formula, "without_variable", REML = FALSE)
      model_with <- self$fit(data, full_formula, random_formula, "with_variable", REML = FALSE)

      # Likelihood ratio test
      lrt <- stats::anova(model_without$model, model_with$model)

      # Calculate differences
      delta_aic <- model_without$aic - model_with$aic  # Positive = with is better
      delta_bic <- model_without$bic - model_with$bic
      bayes_factor <- exp(0.5 * delta_bic)  # BF in favor of simpler model

      # Interpret Bayes factor
      bf_interpretation <- private$.interpret_bayes_factor(bayes_factor)

      list(
        lrt_chisq = if (!is.null(lrt$Chisq)) lrt$Chisq[2] else NA,
        lrt_df = if (!is.null(lrt$Df)) lrt$Df[2] else NA,
        lrt_p_value = if (!is.null(lrt$`Pr(>Chisq)`)) lrt$`Pr(>Chisq)`[2] else NA,
        delta_aic = delta_aic,
        delta_bic = delta_bic,
        bayes_factor = bayes_factor,
        bf_interpretation = bf_interpretation,
        variable_significant = !is.null(lrt$`Pr(>Chisq)`) && lrt$`Pr(>Chisq)`[2] < 0.05,
        model_without = model_without,
        model_with = model_with
      )
    },

    #' @description Test model assumptions
    #' @param model_result LmmModelResult object
    #' @return LmmDiagnostics object
    test_assumptions = function(model_result) {
      model <- model_result$model
      resid <- stats::residuals(model)
      fitted <- stats::fitted(model)

      # Normality test on residuals
      normality_test <- private$.test_normality(resid)

      # Homoscedasticity assessment
      homoscedasticity_test <- private$.test_homoscedasticity(resid, fitted)

      # Extract random effects
      random_effects <- lme4::ranef(model)[[1]]

      LmmDiagnostics$new(
        normality_test = normality_test,
        homoscedasticity_test = homoscedasticity_test,
        residuals = resid,
        fitted = fitted,
        random_effects = random_effects
      )
    },

    #' @description Get data for Q-Q plot
    #' @param diagnostics LmmDiagnostics object
    #' @return Data frame with theoretical and sample quantiles
    qq_plot_data = function(diagnostics) {
      resid <- diagnostics$residuals
      n <- length(resid)
      theoretical <- stats::qnorm(stats::ppoints(n))
      sample <- sort(resid)

      data.frame(
        theoretical = theoretical,
        sample = sample
      )
    },

    #' @description Get data for residuals vs fitted plot
    #' @param diagnostics LmmDiagnostics object
    #' @return Data frame with fitted values and residuals
    residual_plot_data = function(diagnostics) {
      data.frame(
        fitted = diagnostics$fitted,
        residuals = diagnostics$residuals
      )
    },

    #' @description Get fixed effects coefficients
    #' @param model_result LmmModelResult object
    #' @return Data frame with estimates, SE, and p-values
    get_fixed_effects = function(model_result) {
      model <- model_result$model

      # Try to use lmerTest for p-values, fall back to lme4 summary
      coef_table <- tryCatch({
        if (requireNamespace("lmerTest", quietly = TRUE)) {
          lmer_model <- lmerTest::as_lmerModLmerTest(model)
          as.data.frame(summary(lmer_model)$coefficients)
        } else {
          as.data.frame(summary(model)$coefficients)
        }
      }, error = function(e) {
        # lmerTest failed, use basic lme4 summary
        as.data.frame(summary(model)$coefficients)
      })

      coef_table$term <- rownames(coef_table)
      rownames(coef_table) <- NULL

      # Handle different column names from lme4 vs lmerTest
      n_cols <- ncol(coef_table) - 1  # Excluding term

      if (n_cols == 5) {
        # lmerTest output: Estimate, Std. Error, df, t value, Pr(>|t|)
        names(coef_table) <- c("estimate", "std_error", "df", "t_value", "p_value", "term")
      } else if (n_cols == 3) {
        # lme4 output: Estimate, Std. Error, t value (no p-value or df)
        names(coef_table) <- c("estimate", "std_error", "t_value", "term")
        coef_table$df <- NA_real_
        coef_table$p_value <- NA_real_
      } else {
        # Fallback
        names(coef_table)[ncol(coef_table)] <- "term"
        if (!"df" %in% names(coef_table)) coef_table$df <- NA_real_
        if (!"p_value" %in% names(coef_table)) coef_table$p_value <- NA_real_
      }

      coef_table <- coef_table[, c("term", "estimate", "std_error", "df", "t_value", "p_value")]
      coef_table
    },

    #' @description Get random effects
    #' @param model_result LmmModelResult object
    #' @return List with random effect estimates and variances
    get_random_effects = function(model_result) {
      model <- model_result$model

      list(
        estimates = lme4::ranef(model),
        variances = lme4::VarCorr(model)
      )
    },

    #' @description Predict values for new data
    #' @param model_result LmmModelResult object
    #' @param newdata Data frame for prediction
    #' @param include_random Include random effects in prediction
    #' @return Vector of predicted values
    predict_values = function(model_result, newdata, include_random = TRUE) {
      model <- model_result$model

      if (include_random) {
        stats::predict(model, newdata = newdata, allow.new.levels = TRUE)
      } else {
        stats::predict(model, newdata = newdata, re.form = NA)
      }
    }
  ),

  private = list(

    .validate_dependencies = function() {
      if (!requireNamespace("lme4", quietly = TRUE)) {
        warning("Package 'lme4' is recommended for LMM analysis")
      }
    },

    .build_full_formula = function(fixed_formula, random_formula) {
      # Convert to characters
      fixed_str <- deparse(fixed_formula, width.cutoff = 500)
      random_str <- deparse(random_formula, width.cutoff = 500)

      # Remove leading "~" from random formula
      random_str <- gsub("^~\\s*", "", random_str)

      # Combine
      full_str <- paste(fixed_str, "+", "(", random_str, ")")
      stats::as.formula(full_str)
    },

    .calculate_r_squared = function(model) {
      # Use Nakagawa & Schielzeth (2013) R-squared
      tryCatch({
        # Fixed effects variance
        fixed_var <- stats::var(stats::predict(model, re.form = NA))

        # Random effects variance
        vc <- lme4::VarCorr(model)
        random_var <- sum(sapply(vc, function(x) sum(diag(as.matrix(x)))))

        # Residual variance
        resid_var <- attr(vc, "sc")^2

        # Total variance
        total_var <- fixed_var + random_var + resid_var

        list(
          marginal = fixed_var / total_var,
          conditional = (fixed_var + random_var) / total_var
        )
      }, error = function(e) {
        list(marginal = NA_real_, conditional = NA_real_)
      })
    },

    .test_normality = function(residuals) {
      # Shapiro-Wilk test (limited to 5000 observations)
      n <- length(residuals)
      if (n > 5000) {
        sample_resid <- sample(residuals, 5000)
      } else {
        sample_resid <- residuals
      }

      test <- stats::shapiro.test(sample_resid)

      list(
        test_name = "Shapiro-Wilk",
        statistic = test$statistic,
        p_value = test$p.value,
        is_normal = test$p.value > 0.05,
        interpretation = if (test$p.value > 0.05) {
          "Residuals appear normally distributed (p > 0.05)"
        } else {
          "Residuals may deviate from normality (p < 0.05)"
        }
      )
    },

    .test_homoscedasticity = function(residuals, fitted) {
      # Breusch-Pagan-like test using correlation
      abs_resid <- abs(residuals)
      correlation <- stats::cor(abs_resid, fitted)
      test <- stats::cor.test(abs_resid, fitted)

      list(
        test_name = "Residual-Fitted Correlation",
        correlation = correlation,
        p_value = test$p.value,
        is_homoscedastic = abs(correlation) < 0.1 && test$p.value > 0.05,
        interpretation = if (abs(correlation) < 0.1) {
          "Residuals show no strong pattern with fitted values (homoscedastic)"
        } else {
          paste0("Residuals show correlation with fitted values (r = ",
                 round(correlation, 3), "), potential heteroscedasticity")
        }
      )
    },

    .interpret_bayes_factor = function(bf) {
      # Interpret Bayes Factor using Jeffreys (1961) scale
      # BF > 1 favors simpler model, BF < 1 favors complex model
      if (is.na(bf)) return("Unable to calculate")

      if (bf > 100) {
        "Decisive evidence for simpler model"
      } else if (bf > 30) {
        "Very strong evidence for simpler model"
      } else if (bf > 10) {
        "Strong evidence for simpler model"
      } else if (bf > 3) {
        "Moderate evidence for simpler model"
      } else if (bf > 1) {
        "Weak evidence for simpler model"
      } else if (bf > 1 / 3) {
        "Weak evidence for complex model"
      } else if (bf > 1 / 10) {
        "Moderate evidence for complex model"
      } else if (bf > 1 / 30) {
        "Strong evidence for complex model"
      } else if (bf > 1 / 100) {
        "Very strong evidence for complex model"
      } else {
        "Decisive evidence for complex model"
      }
    }
  )
)
