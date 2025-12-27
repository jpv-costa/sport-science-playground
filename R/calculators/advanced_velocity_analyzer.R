# R/calculators/advanced_velocity_analyzer.R
# Service: Advanced Velocity-RIR Analyses for Study 6
#
# =============================================================================
# PRINCIPLES APPLIED (from CLAUDE.md)
# =============================================================================
#
# NAMING PRINCIPLES:
# - Understandability: Domain terms (mvt, icc, decay_rate, failure_prediction)
# - Consistency: All public methods verb-based (analyze, calculate, compare, build)
# - Distinguishability: analyze_mvt_variability vs analyze_velocity_decay
# - Conciseness: mvt for Minimum Velocity Threshold, icc for Intraclass Correlation
#
# SOLID PRINCIPLES:
# - SRP: Each public method handles one hypothesis (H2-H6), one actor (researcher)
# - OCP: New analyses can be added without modifying existing code
# - DIP: Depends on data frames (abstractions), not specific data sources
#
# CUPID PRINCIPLES:
# - Composable: Results are R6 objects that can be combined and compared
# - Unix: Each method does one thing (one hypothesis per method)
# - Predictable: Same data -> same results (deterministic statistical tests)
# - Idiomatic: Follows R conventions, proper roxygen docs, R6 classes
# - Domain-based: Names reflect VBT research concepts (MVT, ICC, decay, prediction)
#
# SCIENTIFIC VALIDITY:
# - ICC(2,1): Two-way random effects, single measurement (Shrout & Fleiss, 1979)
# - Interpretation thresholds from Koo & Li (2016): Poor/Moderate/Good/Excellent
# - SEM and MDC95 formulas follow COSMIN guidelines
# - Model comparison uses AIC with >2 difference threshold (Burnham & Anderson)
#
# TESTABILITY:
# - All private methods are pure functions (input -> output)
# - No side effects in statistical calculations
# - Each method can be tested in isolation with synthetic data
#
# =============================================================================
# ANALYSES IMPLEMENTED (Hypotheses H2-H6)
# =============================================================================
# - H2: Minimum Velocity Threshold (MVT) variability at failure
# - H3: Day-to-day reliability of individual velocity profiles
# - H4: Polynomial vs linear model comparison
# - H5: Velocity decay rate within sets
# - H6: Early rep velocity for failure prediction
# =============================================================================

box::use(
  R6[R6Class],
  stats[lm, anova, coef, predict, residuals, sd, var, qt, t.test, wilcox.test],
  stats[AIC, BIC, confint, cor, IQR, quantile, aggregate]
)

# =============================================================================
# Result Data Classes
# =============================================================================

#' MVT Analysis Result
#'
#' Contains minimum velocity threshold statistics at failure (RIR=0).
#' @export
MvtAnalysisResult <- R6Class(
  classname = "MvtAnalysisResult",

  public = list(
    #' @field population_stats Named list with mean, sd, cv, iqr
    population_stats = NULL,

    #' @field individual_stats Data frame with per-participant MVT
    individual_stats = NULL,

    #' @field sex_comparison Named list with test results
    sex_comparison = NULL,

    #' @field load_comparison Named list with test results
    load_comparison = NULL,

    #' @description Create MVT result
    initialize = function(population_stats, individual_stats, sex_comparison, load_comparison) {
      self$population_stats <- population_stats
      self$individual_stats <- individual_stats
      self$sex_comparison <- sex_comparison
      self$load_comparison <- load_comparison
    }
  )
)

#' Day-to-Day Reliability Result
#'
#' Contains ICC and reliability statistics for velocity profiles.
#' @export
ReliabilityResult <- R6Class(
  classname = "ReliabilityResult",

  public = list(
    #' @field slope_icc Named list with ICC, CI, SEM, MDC
    slope_icc = NULL,

    #' @field intercept_icc Named list with ICC, CI, SEM, MDC
    intercept_icc = NULL,

    #' @field mvt_icc Named list with ICC, CI, SEM, MDC
    mvt_icc = NULL,

    #' @field day_parameters Data frame with Day 1 and Day 2 parameters
    day_parameters = NULL,

    #' @description Create reliability result
    initialize = function(slope_icc, intercept_icc, mvt_icc, day_parameters) {
      self$slope_icc <- slope_icc
      self$intercept_icc <- intercept_icc
      self$mvt_icc <- mvt_icc
      self$day_parameters <- day_parameters
    }
  )
)

#' Model Comparison Result
#'
#' Contains comparison of linear vs polynomial models.
#' @export
ModelComparisonResult <- R6Class(
  classname = "ModelComparisonResult",

  public = list(
    #' @field individual_results Data frame with per-participant metrics
    individual_results = NULL,

    #' @field population_comparison Named list with LMM comparison
    population_comparison = NULL,

    #' @field best_model_summary Named list with model selection summary
    best_model_summary = NULL,

    #' @description Create model comparison result
    initialize = function(individual_results, population_comparison, best_model_summary) {
      self$individual_results <- individual_results
      self$population_comparison <- population_comparison
      self$best_model_summary <- best_model_summary
    }
  )
)

#' Velocity Decay Result
#'
#' Contains velocity decay rate analysis.
#' @export
VelocityDecayResult <- R6Class(
  classname = "VelocityDecayResult",

  public = list(
    #' @field decay_summary Named list with average decay statistics
    decay_summary = NULL,

    #' @field decay_acceleration Named list with acceleration test results
    decay_acceleration = NULL,

    #' @field set_trajectories Data frame with per-set velocity trajectories
    set_trajectories = NULL,

    #' @field breakpoint Named list with changepoint analysis
    breakpoint = NULL,

    #' @description Create velocity decay result
    initialize = function(decay_summary, decay_acceleration, set_trajectories, breakpoint) {
      self$decay_summary <- decay_summary
      self$decay_acceleration <- decay_acceleration
      self$set_trajectories <- set_trajectories
      self$breakpoint <- breakpoint
    }
  )
)

#' Failure Prediction Result
#'
#' Contains early rep prediction model for failure.
#' @export
FailurePredictionResult <- R6Class(
  classname = "FailurePredictionResult",

  public = list(
    #' @field prediction_models Named list of fitted prediction models
    prediction_models = NULL,

    #' @field cv_results Named list with cross-validation metrics
    cv_results = NULL,

    #' @field lookup_table Data frame with practical velocity-to-reps table
    lookup_table = NULL,

    #' @description Create failure prediction result
    initialize = function(prediction_models, cv_results, lookup_table) {
      self$prediction_models <- prediction_models
      self$cv_results <- cv_results
      self$lookup_table <- lookup_table
    }
  )
)

# =============================================================================
# Main Analyzer Class
# =============================================================================

#' Advanced Velocity Analyzer
#'
#' R6 class implementing 5 novel velocity-RIR analyses for Study 6.
#'
#' @export
AdvancedVelocityAnalyzer <- R6Class(
  classname = "AdvancedVelocityAnalyzer",

  public = list(

    # =========================================================================
    # H2: Minimum Velocity Threshold (MVT) Variability
    # =========================================================================

    #' @description Analyze MVT variability at failure (RIR=0)
    #'
    #' Research Question: How variable is the velocity at muscular failure?
    #'
    #' @param data Data frame with columns: id, sex, load_percentage, rir, mean_velocity
    #' @return MvtAnalysisResult object
    analyze_mvt_variability = function(data) {
      # Extract failure reps only (RIR = 0)
      failure_data <- data[data$rir == 0, ]

      if (nrow(failure_data) == 0) {
        stop("No failure observations (RIR=0) found in data")
      }

      # Population-level statistics
      population_stats <- private$.calculate_population_mvt(failure_data)

      # Individual-level statistics
      individual_stats <- private$.calculate_individual_mvt(failure_data)

      # Sex comparison
      sex_comparison <- private$.compare_mvt_by_sex(failure_data)

      # Load comparison
      load_comparison <- private$.compare_mvt_by_load(failure_data)

      MvtAnalysisResult$new(
        population_stats = population_stats,
        individual_stats = individual_stats,
        sex_comparison = sex_comparison,
        load_comparison = load_comparison
      )
    },

    # =========================================================================
    # H3: Day-to-Day Reliability
    # =========================================================================

    #' @description Calculate day-to-day reliability of individual velocity profiles
    #'
    #' Research Question: Are individual velocity-RIR relationships stable across days?
    #'
    #' @param data Data frame with columns: id, day, rir, mean_velocity
    #' @param load Optional load to filter by (e.g., "90%")
    #' @return ReliabilityResult object
    calculate_day_reliability = function(data, load = NULL) {
      # Filter by load if specified
      if (!is.null(load)) {
        data <- data[data$load_percentage == load, ]
      }

      # Fit individual linear models for each day
      day_parameters <- private$.fit_individual_day_models(data)

      if (nrow(day_parameters) == 0) {
        stop("No participants with data on both days")
      }

      # Calculate ICC for slopes
      slope_icc <- private$.calculate_icc(
        day_parameters$slope_day1,
        day_parameters$slope_day2,
        "slope"
      )

      # Calculate ICC for intercepts
      intercept_icc <- private$.calculate_icc(
        day_parameters$intercept_day1,
        day_parameters$intercept_day2,
        "intercept"
      )

      # Calculate ICC for predicted MVT (velocity at RIR=0)
      mvt_icc <- private$.calculate_icc(
        day_parameters$mvt_day1,
        day_parameters$mvt_day2,
        "mvt"
      )

      ReliabilityResult$new(
        slope_icc = slope_icc,
        intercept_icc = intercept_icc,
        mvt_icc = mvt_icc,
        day_parameters = day_parameters
      )
    },

    # =========================================================================
    # H4: Polynomial vs Linear Model Comparison
    # =========================================================================

    #' @description Compare polynomial vs linear models for velocity-RIR
    #'
    #' Research Question: Does a quadratic model fit better than linear?
    #'
    #' @param data Data frame with columns: id, rir, mean_velocity
    #' @return ModelComparisonResult object
    compare_polynomial_models = function(data) {
      # Individual-level model comparison
      individual_results <- private$.compare_individual_models(data)

      # Population-level LMM comparison
      population_comparison <- private$.compare_population_models(data)

      # Summarize best model selection
      best_model_summary <- private$.summarize_best_model(
        individual_results,
        population_comparison
      )

      ModelComparisonResult$new(
        individual_results = individual_results,
        population_comparison = population_comparison,
        best_model_summary = best_model_summary
      )
    },

    # =========================================================================
    # H5: Velocity Decay Rate Within Sets
    # =========================================================================

    #' @description Analyze velocity decay rate within sets
    #'
    #' Research Question: How does velocity loss per rep change as a set progresses?
    #'
    #' @param data Data frame with columns: set_id, rep_number, mean_velocity
    #' @return VelocityDecayResult object
    analyze_velocity_decay = function(data) {
      # Validate required columns
      required_cols <- c("set_id", "rep_number", "mean_velocity")
      missing_cols <- setdiff(required_cols, names(data))
      if (length(missing_cols) > 0) {
        stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
      }

      # Calculate per-set trajectories and decay
      set_trajectories <- private$.calculate_set_trajectories(data)

      # Summary of decay rates
      decay_summary <- private$.summarize_decay_rates(set_trajectories)

      # Test if decay accelerates (slope of delta_v vs rep_number)
      decay_acceleration <- private$.test_decay_acceleration(set_trajectories)

      # Changepoint detection
      breakpoint <- private$.detect_decay_breakpoint(set_trajectories)

      VelocityDecayResult$new(
        decay_summary = decay_summary,
        decay_acceleration = decay_acceleration,
        set_trajectories = set_trajectories,
        breakpoint = breakpoint
      )
    },

    # =========================================================================
    # H6: Early Rep Failure Prediction
    # =========================================================================

    #' @description Build prediction model for failure from early rep velocities
    #'
    #' Research Question: Can first 1-3 rep velocities predict when failure occurs?
    #'
    #' @param data Data frame with columns: set_id, rep_number, mean_velocity,
    #'             load_percentage, reps_to_failure
    #' @return FailurePredictionResult object
    build_failure_predictor = function(data) {
      # Validate required columns
      required_cols <- c("set_id", "rep_number", "mean_velocity",
                         "load_percentage", "reps_to_failure")
      missing_cols <- setdiff(required_cols, names(data))
      if (length(missing_cols) > 0) {
        stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
      }

      # Prepare set-level data with early rep velocities
      set_data <- private$.prepare_set_level_data(data)

      # Fit prediction models
      prediction_models <- private$.fit_prediction_models(set_data)

      # Cross-validation
      cv_results <- private$.cross_validate_predictions(set_data)

      # Create practical lookup table
      lookup_table <- private$.create_prediction_lookup(set_data, prediction_models)

      FailurePredictionResult$new(
        prediction_models = prediction_models,
        cv_results = cv_results,
        lookup_table = lookup_table
      )
    }
  ),

  private = list(

    # =========================================================================
    # H2: MVT Private Methods
    # =========================================================================

    .calculate_population_mvt = function(failure_data) {
      velocities <- failure_data$mean_velocity

      mean_mvt <- mean(velocities, na.rm = TRUE)
      sd_mvt <- sd(velocities, na.rm = TRUE)
      cv_mvt <- (sd_mvt / mean_mvt) * 100
      iqr_mvt <- IQR(velocities, na.rm = TRUE)
      q25 <- quantile(velocities, 0.25, na.rm = TRUE)
      q75 <- quantile(velocities, 0.75, na.rm = TRUE)
      min_mvt <- min(velocities, na.rm = TRUE)
      max_mvt <- max(velocities, na.rm = TRUE)
      n_observations <- length(velocities)

      list(
        mean = mean_mvt,
        sd = sd_mvt,
        cv_percent = cv_mvt,
        iqr = iqr_mvt,
        q25 = q25,
        q75 = q75,
        min = min_mvt,
        max = max_mvt,
        n = n_observations
      )
    },

    .calculate_individual_mvt = function(failure_data) {
      individual_mvt <- aggregate(
        mean_velocity ~ id + sex,
        data = failure_data,
        FUN = function(x) c(mean = mean(x), sd = sd(x), n = length(x))
      )

      # Flatten the matrix column
      mvt_values <- individual_mvt$mean_velocity
      individual_mvt$mean_mvt <- mvt_values[, "mean"]
      individual_mvt$sd_mvt <- mvt_values[, "sd"]
      individual_mvt$n_sets <- mvt_values[, "n"]
      individual_mvt$mean_velocity <- NULL

      individual_mvt
    },

    .compare_mvt_by_sex = function(failure_data) {
      male_mvt <- failure_data$mean_velocity[failure_data$sex == "male"]
      female_mvt <- failure_data$mean_velocity[failure_data$sex == "female"]

      # Use Mann-Whitney for robustness
      if (length(female_mvt) >= 3 && length(male_mvt) >= 3) {
        test_result <- wilcox.test(male_mvt, female_mvt)
        t_result <- t.test(male_mvt, female_mvt)

        list(
          male_mean = mean(male_mvt, na.rm = TRUE),
          male_sd = sd(male_mvt, na.rm = TRUE),
          female_mean = mean(female_mvt, na.rm = TRUE),
          female_sd = sd(female_mvt, na.rm = TRUE),
          difference = mean(male_mvt, na.rm = TRUE) - mean(female_mvt, na.rm = TRUE),
          wilcox_p = test_result$p.value,
          t_test_p = t_result$p.value,
          significant = test_result$p.value < 0.05
        )
      } else {
        list(
          male_mean = mean(male_mvt, na.rm = TRUE),
          female_mean = mean(female_mvt, na.rm = TRUE),
          note = "Insufficient female observations for formal testing"
        )
      }
    },

    .compare_mvt_by_load = function(failure_data) {
      load_80 <- failure_data$mean_velocity[failure_data$load_percentage == "80%"]
      load_90 <- failure_data$mean_velocity[failure_data$load_percentage == "90%"]

      test_result <- wilcox.test(load_80, load_90, paired = FALSE)
      t_result <- t.test(load_80, load_90)

      list(
        load_80_mean = mean(load_80, na.rm = TRUE),
        load_80_sd = sd(load_80, na.rm = TRUE),
        load_90_mean = mean(load_90, na.rm = TRUE),
        load_90_sd = sd(load_90, na.rm = TRUE),
        difference = mean(load_80, na.rm = TRUE) - mean(load_90, na.rm = TRUE),
        wilcox_p = test_result$p.value,
        t_test_p = t_result$p.value,
        significant = test_result$p.value < 0.05
      )
    },

    # =========================================================================
    # H3: Reliability Private Methods
    # =========================================================================

    .fit_individual_day_models = function(data) {
      participants <- unique(data$id)
      results <- list()

      for (pid in participants) {
        pid_data <- data[data$id == pid, ]

        day1_data <- pid_data[pid_data$day == "Day 1", ]
        day2_data <- pid_data[pid_data$day == "Day 2", ]

        if (nrow(day1_data) < 3 || nrow(day2_data) < 3) {
          next
        }

        # Fit linear models for each day
        model_day1 <- lm(mean_velocity ~ rir, data = day1_data)
        model_day2 <- lm(mean_velocity ~ rir, data = day2_data)

        coef1 <- coef(model_day1)
        coef2 <- coef(model_day2)

        results[[pid]] <- data.frame(
          id = pid,
          slope_day1 = coef1["rir"],
          intercept_day1 = coef1["(Intercept)"],
          mvt_day1 = predict(model_day1, newdata = data.frame(rir = 0)),
          slope_day2 = coef2["rir"],
          intercept_day2 = coef2["(Intercept)"],
          mvt_day2 = predict(model_day2, newdata = data.frame(rir = 0)),
          stringsAsFactors = FALSE
        )
      }

      if (length(results) == 0) {
        return(data.frame())
      }

      do.call(rbind, results)
    },

    #' Calculate ICC(2,1) - Two-way random effects, single measurement
    #' Orchestrates ICC calculation using smaller focused methods
    #' Based on Shrout & Fleiss (1979)
    .calculate_icc = function(values1, values2, name) {
      n <- length(values1)

      if (!private$.is_sufficient_for_icc(n)) {
        return(private$.create_insufficient_icc_result())
      }

      ratings <- cbind(values1, values2)
      variance <- private$.calculate_icc_variance_components(ratings)
      icc <- private$.calculate_icc_value(variance)
      reliability <- private$.calculate_reliability_measures(icc, variance$ms_within)
      ci <- private$.calculate_icc_confidence_interval(icc, n)
      interpretation <- private$.interpret_icc(icc)

      list(
        icc = icc,
        ci_lower = ci$lower,
        ci_upper = ci$upper,
        sem = reliability$sem,
        mdc95 = reliability$mdc95,
        interpretation = interpretation,
        n = n
      )
    },

    #' Check if sample size is sufficient for ICC calculation
    .is_sufficient_for_icc = function(n) {
      n >= 3
    },

    #' Create result object for insufficient data
    .create_insufficient_icc_result = function() {
      list(
        icc = NA,
        ci_lower = NA,
        ci_upper = NA,
        sem = NA,
        mdc95 = NA,
        interpretation = "Insufficient data"
      )
    },

    #' Calculate variance components for ICC
    .calculate_icc_variance_components = function(ratings) {
      subject_means <- rowMeans(ratings)
      ms_between <- var(subject_means) * 2  # k=2 raters
      ms_within <- mean(apply(ratings, 1, var))

      list(
        ms_between = ms_between,
        ms_within = ms_within
      )
    },

    #' Calculate ICC value from variance components
    .calculate_icc_value = function(variance) {
      icc <- (variance$ms_between - variance$ms_within) /
        (variance$ms_between + variance$ms_within)
      max(0, min(1, icc))  # Bound between 0 and 1
    },

    #' Calculate SEM and MDC95 (COSMIN guidelines)
    .calculate_reliability_measures = function(icc, ms_within) {
      pooled_sd <- sqrt(ms_within)
      sem <- pooled_sd * sqrt(1 - icc)
      mdc95 <- sem * 1.96 * sqrt(2)

      list(sem = sem, mdc95 = mdc95)
    },

    #' Calculate confidence interval for ICC
    .calculate_icc_confidence_interval = function(icc, n) {
      se_icc <- sqrt(2 * (1 - icc)^2 * (1 + icc)^2 / (n * (n - 1)))
      list(
        lower = max(0, icc - 1.96 * se_icc),
        upper = min(1, icc + 1.96 * se_icc)
      )
    },

    #' Interpret ICC based on Koo & Li (2016) thresholds
    .interpret_icc = function(icc) {
      if (icc < 0.50) "Poor"
      else if (icc < 0.75) "Moderate"
      else if (icc < 0.90) "Good"
      else "Excellent"
    },

    # =========================================================================
    # H4: Model Comparison Private Methods
    # =========================================================================

    #' Compare linear vs polynomial models for each participant
    #' Orchestrates individual model comparison using smaller focused methods
    .compare_individual_models = function(data) {
      participants <- unique(data$id)
      results <- list()

      for (pid in participants) {
        pid_data <- data[data$id == pid, ]

        if (!private$.has_sufficient_observations(pid_data, min_obs = 4)) {
          next
        }

        models <- private$.fit_velocity_models(pid_data)
        metrics <- private$.calculate_model_comparison_metrics(models)
        best_model <- private$.select_best_model(metrics)

        results[[pid]] <- private$.build_model_comparison_result(
          pid, nrow(pid_data), metrics, best_model
        )
      }

      do.call(rbind, results)
    },

    #' Check if participant has sufficient observations
    .has_sufficient_observations = function(data, min_obs) {
      nrow(data) >= min_obs
    },

    #' Fit linear and quadratic velocity-RIR models
    .fit_velocity_models = function(data) {
      list(
        linear = lm(mean_velocity ~ rir, data = data),
        quadratic = lm(mean_velocity ~ rir + I(rir^2), data = data)
      )
    },

    #' Calculate comparison metrics for fitted models
    .calculate_model_comparison_metrics = function(models) {
      linear <- models$linear
      quad <- models$quadratic

      list(
        rmse_linear = sqrt(mean(residuals(linear)^2)),
        rmse_quad = sqrt(mean(residuals(quad)^2)),
        r2_adj_linear = summary(linear)$adj.r.squared,
        r2_adj_quad = summary(quad)$adj.r.squared,
        aic_linear = AIC(linear),
        aic_quad = AIC(quad),
        bic_linear = BIC(linear),
        bic_quad = BIC(quad)
      )
    },

    #' Select best model based on AIC (>2 threshold per Burnham & Anderson)
    .select_best_model = function(metrics) {
      if (metrics$aic_quad < metrics$aic_linear - 2) "quadratic" else "linear"
    },

    #' Build result data frame for model comparison
    .build_model_comparison_result = function(pid, n_obs, metrics, best_model) {
      data.frame(
        id = pid,
        n_obs = n_obs,
        r2_adj_linear = metrics$r2_adj_linear,
        r2_adj_quad = metrics$r2_adj_quad,
        rmse_linear = metrics$rmse_linear,
        rmse_quad = metrics$rmse_quad,
        aic_linear = metrics$aic_linear,
        aic_quad = metrics$aic_quad,
        bic_linear = metrics$bic_linear,
        bic_quad = metrics$bic_quad,
        delta_aic = metrics$aic_quad - metrics$aic_linear,
        best_model = best_model,
        stringsAsFactors = FALSE
      )
    },

    .compare_population_models = function(data) {
      # Check if lme4 is available
      if (!requireNamespace("lme4", quietly = TRUE)) {
        return(list(note = "lme4 package required for population-level comparison"))
      }

      # Fit population LMM models
      model_linear <- lme4::lmer(
        mean_velocity ~ rir + (1 + rir | id),
        data = data,
        REML = FALSE
      )

      model_quad <- lme4::lmer(
        mean_velocity ~ rir + I(rir^2) + (1 + rir | id),
        data = data,
        REML = FALSE
      )

      # Likelihood ratio test
      lrt <- anova(model_linear, model_quad)

      list(
        aic_linear = AIC(model_linear),
        aic_quad = AIC(model_quad),
        bic_linear = BIC(model_linear),
        bic_quad = BIC(model_quad),
        lrt_chisq = lrt$Chisq[2],
        lrt_df = lrt$Df[2],
        lrt_p = lrt$`Pr(>Chisq)`[2],
        quad_significant = lrt$`Pr(>Chisq)`[2] < 0.05
      )
    },

    .summarize_best_model = function(individual_results, population_comparison) {
      if (is.null(individual_results) || nrow(individual_results) == 0) {
        return(list(note = "No individual results available"))
      }

      n_linear <- sum(individual_results$best_model == "linear")
      n_quad <- sum(individual_results$best_model == "quadratic")
      n_total <- nrow(individual_results)

      # Average improvement from quadratic
      avg_r2_improvement <- mean(
        individual_results$r2_adj_quad - individual_results$r2_adj_linear,
        na.rm = TRUE
      )

      list(
        n_participants = n_total,
        n_linear_best = n_linear,
        n_quad_best = n_quad,
        pct_quad_best = (n_quad / n_total) * 100,
        avg_r2_improvement = avg_r2_improvement,
        population_quad_significant = if (!is.null(population_comparison$quad_significant))
          population_comparison$quad_significant else NA,
        recommendation = if (n_quad / n_total > 0.5 || avg_r2_improvement > 0.05)
          "Consider quadratic model" else "Linear model adequate"
      )
    },

    # =========================================================================
    # H5: Velocity Decay Private Methods
    # =========================================================================

    .calculate_set_trajectories = function(data) {
      sets <- unique(data$set_id)
      trajectories <- list()

      for (sid in sets) {
        set_data <- data[data$set_id == sid, ]
        set_data <- set_data[order(set_data$rep_number), ]

        n_reps <- nrow(set_data)
        if (n_reps < 2) {
          next
        }

        # Calculate rep-to-rep velocity loss
        velocities <- set_data$mean_velocity
        delta_v <- diff(velocities)  # v(n) - v(n-1)
        rep_numbers <- set_data$rep_number[-1]  # Rep number for each delta

        # First rep velocity
        v1 <- velocities[1]

        # Cumulative velocity loss
        cumulative_loss <- v1 - velocities

        trajectories[[sid]] <- data.frame(
          set_id = sid,
          rep_number = rep_numbers,
          velocity = velocities[-1],
          delta_v = delta_v,
          cumulative_loss = cumulative_loss[-1],
          first_rep_velocity = v1,
          n_reps = n_reps,
          stringsAsFactors = FALSE
        )
      }

      do.call(rbind, trajectories)
    },

    .summarize_decay_rates = function(set_trajectories) {
      if (is.null(set_trajectories) || nrow(set_trajectories) == 0) {
        return(list(note = "No trajectory data available"))
      }

      # Average decay per rep
      avg_decay <- mean(set_trajectories$delta_v, na.rm = TRUE)
      sd_decay <- sd(set_trajectories$delta_v, na.rm = TRUE)

      # Decay in early vs late reps
      early_reps <- set_trajectories[set_trajectories$rep_number <= 3, ]
      late_reps <- set_trajectories[set_trajectories$rep_number >= 4, ]

      early_decay <- mean(early_reps$delta_v, na.rm = TRUE)
      late_decay <- mean(late_reps$delta_v, na.rm = TRUE)

      list(
        avg_decay_per_rep = avg_decay,
        sd_decay = sd_decay,
        early_reps_decay = early_decay,
        late_reps_decay = late_decay,
        decay_acceleration = late_decay - early_decay,
        n_observations = nrow(set_trajectories)
      )
    },

    .test_decay_acceleration = function(set_trajectories) {
      if (is.null(set_trajectories) || nrow(set_trajectories) < 10) {
        return(list(note = "Insufficient data for acceleration test"))
      }

      # Model: delta_v ~ rep_number
      # If slope is negative, decay accelerates (velocity drops faster as set progresses)
      model <- lm(delta_v ~ rep_number, data = set_trajectories)
      coefs <- summary(model)$coefficients

      slope <- coefs["rep_number", "Estimate"]
      se <- coefs["rep_number", "Std. Error"]
      t_value <- coefs["rep_number", "t value"]
      p_value <- coefs["rep_number", "Pr(>|t|)"]

      list(
        slope = slope,
        se = se,
        t_value = t_value,
        p_value = p_value,
        accelerating = slope < 0 && p_value < 0.05,
        interpretation = if (slope < 0 && p_value < 0.05) {
          "Decay accelerates toward failure"
        } else if (slope > 0 && p_value < 0.05) {
          "Decay decelerates (unexpected)"
        } else {
          "Decay rate is constant"
        }
      )
    },

    .detect_decay_breakpoint = function(set_trajectories) {
      if (is.null(set_trajectories) || nrow(set_trajectories) < 10) {
        return(list(note = "Insufficient data for breakpoint detection"))
      }

      # Aggregate delta_v by rep_number across all sets
      agg <- aggregate(delta_v ~ rep_number, data = set_trajectories, FUN = mean)
      agg <- agg[order(agg$rep_number), ]

      if (nrow(agg) < 3) {
        return(list(note = "Too few rep positions for breakpoint detection"))
      }

      # Simple breakpoint: where decay becomes notably larger
      # Find rep where decay drops below average by >1 SD
      avg_decay <- mean(agg$delta_v)
      sd_decay <- sd(agg$delta_v)
      threshold <- avg_decay - sd_decay

      breakpoint_rep <- NA
      for (i in seq_len(nrow(agg))) {
        if (agg$delta_v[i] < threshold) {
          breakpoint_rep <- agg$rep_number[i]
          break
        }
      }

      list(
        breakpoint_rep = breakpoint_rep,
        avg_decay = avg_decay,
        threshold = threshold,
        decay_by_rep = agg,
        interpretation = if (!is.na(breakpoint_rep)) {
          paste0("Accelerated decay starts at rep ", breakpoint_rep)
        } else {
          "No clear breakpoint detected"
        }
      )
    },

    # =========================================================================
    # H6: Prediction Private Methods
    # =========================================================================

    .prepare_set_level_data = function(data) {
      sets <- unique(data$set_id)
      set_data <- list()

      for (sid in sets) {
        s <- data[data$set_id == sid, ]
        s <- s[order(s$rep_number), ]

        if (nrow(s) < 3) {
          next
        }

        # Extract first 3 rep velocities
        v1 <- s$mean_velocity[1]
        v2 <- if (nrow(s) >= 2) s$mean_velocity[2] else NA
        v3 <- if (nrow(s) >= 3) s$mean_velocity[3] else NA

        set_data[[sid]] <- data.frame(
          set_id = sid,
          id = s$id[1],
          load_percentage = s$load_percentage[1],
          v1 = v1,
          v2 = v2,
          v3 = v3,
          reps_to_failure = s$reps_to_failure[1],
          stringsAsFactors = FALSE
        )
      }

      do.call(rbind, set_data)
    },

    .fit_prediction_models = function(set_data) {
      models <- list()

      # Model 1: n_reps ~ v1
      models$v1_only <- lm(reps_to_failure ~ v1, data = set_data)

      # Model 2: n_reps ~ v1 + v2
      complete_v2 <- set_data[!is.na(set_data$v2), ]
      if (nrow(complete_v2) > 5) {
        models$v1_v2 <- lm(reps_to_failure ~ v1 + v2, data = complete_v2)
      }

      # Model 3: n_reps ~ v1 + v2 + v3
      complete_v3 <- set_data[!is.na(set_data$v3), ]
      if (nrow(complete_v3) > 5) {
        models$v1_v2_v3 <- lm(reps_to_failure ~ v1 + v2 + v3, data = complete_v3)
      }

      # Model 4: n_reps ~ v1 + load
      # Convert load to numeric
      set_data$load_numeric <- ifelse(set_data$load_percentage == "90%", 90, 80)
      models$v1_load <- lm(reps_to_failure ~ v1 + load_numeric, data = set_data)

      models
    },

    .cross_validate_predictions = function(set_data) {
      # Leave-one-set-out cross-validation
      n <- nrow(set_data)
      predictions <- numeric(n)
      actuals <- set_data$reps_to_failure

      for (i in seq_len(n)) {
        train <- set_data[-i, ]
        test <- set_data[i, ]

        model <- lm(reps_to_failure ~ v1, data = train)
        predictions[i] <- predict(model, newdata = test)
      }

      errors <- actuals - predictions
      mae <- mean(abs(errors))
      rmse <- sqrt(mean(errors^2))
      r2 <- cor(actuals, predictions)^2

      # Within 1 rep accuracy
      within_1_rep <- mean(abs(errors) <= 1) * 100
      within_2_reps <- mean(abs(errors) <= 2) * 100

      list(
        mae = mae,
        rmse = rmse,
        r2 = r2,
        within_1_rep_pct = within_1_rep,
        within_2_reps_pct = within_2_reps,
        n_sets = n
      )
    },

    .create_prediction_lookup = function(set_data, prediction_models) {
      # Create lookup table: v1 -> expected reps
      v1_range <- seq(0.15, 0.55, by = 0.05)

      model <- prediction_models$v1_only
      if (is.null(model)) {
        return(data.frame(note = "No prediction model available"))
      }

      lookup <- data.frame(
        first_rep_velocity = v1_range,
        predicted_reps = predict(model, newdata = data.frame(v1 = v1_range))
      )

      # Add confidence intervals
      ci <- predict(model, newdata = data.frame(v1 = v1_range),
                    interval = "prediction", level = 0.95)
      lookup$lower_95 <- ci[, "lwr"]
      lookup$upper_95 <- ci[, "upr"]

      # Round to practical values
      lookup$predicted_reps <- pmax(1, round(lookup$predicted_reps, 1))
      lookup$lower_95 <- pmax(1, round(lookup$lower_95, 1))
      lookup$upper_95 <- round(lookup$upper_95, 1)

      lookup
    }
  )
)
