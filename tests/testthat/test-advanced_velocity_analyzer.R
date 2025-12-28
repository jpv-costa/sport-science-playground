# tests/testthat/test-advanced_velocity_analyzer.R
# Unit tests for AdvancedVelocityAnalyzer (Study 6: H2-H6)

box::use(
  testthat[...],
  ../../R/calculators/advanced_velocity_analyzer[
    AdvancedVelocityAnalyzer,
    MvtAnalysisResult,
    ReliabilityResult,
    ModelComparisonResult,
    VelocityDecayResult,
    FailurePredictionResult
  ]
)

# ==============================================================================
# Test Fixtures
# ==============================================================================

#' Create MVT test data with failure observations (RIR=0)
create_mvt_data <- function(n_participants = 10, n_sets_per = 4, seed = 42) {
  set.seed(seed)

  data <- expand.grid(
    participant_num = 1:n_participants,
    set_num = 1:n_sets_per
  )

  data$id <- paste0("P", data$participant_num)
  data$sex <- ifelse(data$participant_num <= n_participants / 2, "male", "female")
  data$load_percentage <- ifelse(data$set_num <= 2, "80%", "90%")
  data$rir <- 0  # All failure observations

  # MVT varies by sex and load
  base_mvt <- 0.20
  sex_effect <- ifelse(data$sex == "male", 0.02, 0)
  load_effect <- ifelse(data$load_percentage == "90%", -0.03, 0)
  participant_effect <- rnorm(nrow(data), 0, 0.02)

  data$mean_velocity <- base_mvt + sex_effect + load_effect + participant_effect
  data$mean_velocity <- pmax(0.05, data$mean_velocity)  # Floor at 0.05 m/s

  data[, c("id", "sex", "load_percentage", "rir", "mean_velocity")]
}

#' Create reliability test data with Day 1 and Day 2 observations
create_reliability_data <- function(n_participants = 8, n_obs_per_day = 8, seed = 42) {
  set.seed(seed)

  data <- expand.grid(
    participant_num = 1:n_participants,
    day = c("Day 1", "Day 2"),
    obs = 1:n_obs_per_day
  )

  data$id <- paste0("P", data$participant_num)
  data$rir <- rep(0:7, length.out = nrow(data))
  data$load_percentage <- "90%"

  # Individual slope and intercept
  participant_effects <- data.frame(
    id = paste0("P", 1:n_participants),
    intercept = 0.20 + rnorm(n_participants, 0, 0.02),
    slope = 0.04 + rnorm(n_participants, 0, 0.005)
  )

  data <- merge(data, participant_effects, by = "id")

  # Velocity = intercept + slope * RIR + day effect + noise
  day_effect <- ifelse(data$day == "Day 2", rnorm(nrow(data), 0, 0.01), 0)
  data$mean_velocity <- data$intercept + data$slope * data$rir + day_effect +
    rnorm(nrow(data), 0, 0.015)

  data[, c("id", "day", "rir", "load_percentage", "mean_velocity")]
}

#' Create polynomial model comparison data
create_polynomial_data <- function(n_participants = 10, n_obs_per = 10, seed = 42) {
  set.seed(seed)

  data <- expand.grid(
    participant_num = 1:n_participants,
    obs = 1:n_obs_per
  )

  data$id <- paste0("P", data$participant_num)
  data$rir <- sample(0:7, nrow(data), replace = TRUE)

  # Generate velocity with slight curvature
  # v = 0.2 + 0.04*rir + 0.002*rir^2 + noise
  participant_effects <- data.frame(
    id = paste0("P", 1:n_participants),
    intercept_effect = rnorm(n_participants, 0, 0.02)
  )

  data <- merge(data, participant_effects, by = "id")

  data$mean_velocity <- 0.20 +
    0.04 * data$rir +
    0.002 * data$rir^2 +
    data$intercept_effect +
    rnorm(nrow(data), 0, 0.015)

  data[, c("id", "rir", "mean_velocity")]
}

#' Create velocity decay test data with set trajectories
create_decay_data <- function(n_sets = 20, max_reps = 10, seed = 42) {
  set.seed(seed)

  sets <- list()

  for (i in 1:n_sets) {
    n_reps <- sample(5:max_reps, 1)
    v1 <- runif(1, 0.35, 0.50)

    # Decay accelerates toward failure
    rep_numbers <- 1:n_reps
    decay_rate <- 0.025 + 0.003 * rep_numbers
    velocities <- v1 - cumsum(c(0, decay_rate[-n_reps])) + rnorm(n_reps, 0, 0.01)
    velocities <- pmax(0.10, velocities)

    sets[[i]] <- data.frame(
      set_id = paste0("S", i),
      rep_number = rep_numbers,
      mean_velocity = velocities
    )
  }

  do.call(rbind, sets)
}

#' Create failure prediction test data
create_prediction_data <- function(n_sets = 30, seed = 42) {
  set.seed(seed)

  sets <- list()

  for (i in 1:n_sets) {
    # Higher v1 -> more reps to failure
    v1 <- runif(1, 0.25, 0.55)
    reps_to_failure <- round(5 + 15 * (v1 - 0.25) / 0.30 + rnorm(1, 0, 1))
    reps_to_failure <- max(3, min(20, reps_to_failure))

    n_reps <- reps_to_failure
    velocities <- seq(v1, 0.15, length.out = n_reps) + rnorm(n_reps, 0, 0.01)

    sets[[i]] <- data.frame(
      set_id = paste0("S", i),
      id = paste0("P", sample(1:5, 1)),
      rep_number = 1:n_reps,
      mean_velocity = velocities,
      load_percentage = sample(c("80%", "90%"), 1),
      reps_to_failure = reps_to_failure
    )
  }

  do.call(rbind, sets)
}

# ==============================================================================
# Result Class Tests
# ==============================================================================

describe("MvtAnalysisResult", {

  describe("initialization", {
    it("stores all fields correctly", {
      result <- MvtAnalysisResult$new(
        population_stats = list(mean = 0.20, sd = 0.03),
        individual_stats = data.frame(id = "P1", mean_mvt = 0.20),
        sex_comparison = list(male_mean = 0.22, female_mean = 0.18),
        load_comparison = list(load_80_mean = 0.22, load_90_mean = 0.18)
      )

      expect_s3_class(result, "MvtAnalysisResult")
      expect_equal(result$population_stats$mean, 0.20)
      expect_equal(result$sex_comparison$male_mean, 0.22)
    })
  })
})

describe("ReliabilityResult", {

  describe("initialization", {
    it("stores ICC results correctly", {
      result <- ReliabilityResult$new(
        slope_icc = list(icc = 0.85, interpretation = "Good"),
        intercept_icc = list(icc = 0.90, interpretation = "Excellent"),
        mvt_icc = list(icc = 0.80, interpretation = "Good"),
        day_parameters = data.frame(id = "P1", slope_day1 = 0.04)
      )

      expect_s3_class(result, "ReliabilityResult")
      expect_equal(result$slope_icc$icc, 0.85)
      expect_equal(result$intercept_icc$interpretation, "Excellent")
    })
  })
})

describe("ModelComparisonResult", {

  describe("initialization", {
    it("stores model comparison fields", {
      result <- ModelComparisonResult$new(
        individual_results = data.frame(id = "P1", best_model = "linear"),
        population_comparison = list(aic_linear = 100, aic_quad = 102),
        best_model_summary = list(recommendation = "Linear model adequate")
      )

      expect_s3_class(result, "ModelComparisonResult")
      expect_equal(result$best_model_summary$recommendation, "Linear model adequate")
    })
  })
})

describe("VelocityDecayResult", {

  describe("initialization", {
    it("stores decay analysis fields", {
      result <- VelocityDecayResult$new(
        decay_summary = list(avg_decay_per_rep = -0.03),
        decay_acceleration = list(accelerating = TRUE),
        set_trajectories = data.frame(set_id = "S1", rep_number = 2),
        breakpoint = list(breakpoint_rep = 5)
      )

      expect_s3_class(result, "VelocityDecayResult")
      expect_true(result$decay_acceleration$accelerating)
    })
  })
})

describe("FailurePredictionResult", {

  describe("initialization", {
    it("stores prediction fields", {
      result <- FailurePredictionResult$new(
        prediction_models = list(v1_only = NULL),
        cv_results = list(mae = 1.5, r2 = 0.75),
        lookup_table = data.frame(first_rep_velocity = 0.30, predicted_reps = 8)
      )

      expect_s3_class(result, "FailurePredictionResult")
      expect_equal(result$cv_results$mae, 1.5)
    })
  })
})

# ==============================================================================
# AdvancedVelocityAnalyzer Tests
# ==============================================================================

describe("AdvancedVelocityAnalyzer", {

  describe("initialization", {
    it("creates analyzer without error", {
      analyzer <- AdvancedVelocityAnalyzer$new()
      expect_s3_class(analyzer, "AdvancedVelocityAnalyzer")
    })
  })

  # ============================================================================
  # H2: MVT Variability Tests
  # ============================================================================

  describe("analyze_mvt_variability (H2)", {

    it("returns MvtAnalysisResult object", {
      analyzer <- AdvancedVelocityAnalyzer$new()
      data <- create_mvt_data()

      result <- analyzer$analyze_mvt_variability(data)

      expect_s3_class(result, "MvtAnalysisResult")
    })

    it("calculates population statistics correctly", {
      analyzer <- AdvancedVelocityAnalyzer$new()
      data <- create_mvt_data()

      result <- analyzer$analyze_mvt_variability(data)
      pop_stats <- result$population_stats

      expect_true("mean" %in% names(pop_stats))
      expect_true("sd" %in% names(pop_stats))
      expect_true("cv_percent" %in% names(pop_stats))
      expect_true("iqr" %in% names(pop_stats))
      expect_true("n" %in% names(pop_stats))

      # Mean MVT should be around 0.20 m/s
      expect_true(pop_stats$mean > 0.10 && pop_stats$mean < 0.35)
    })

    it("calculates individual participant MVT", {
      analyzer <- AdvancedVelocityAnalyzer$new()
      data <- create_mvt_data(n_participants = 5)

      result <- analyzer$analyze_mvt_variability(data)
      ind_stats <- result$individual_stats

      expect_true("id" %in% names(ind_stats))
      expect_true("mean_mvt" %in% names(ind_stats))
      expect_equal(length(unique(ind_stats$id)), 5)
    })

    it("compares MVT by sex", {
      analyzer <- AdvancedVelocityAnalyzer$new()
      data <- create_mvt_data()

      result <- analyzer$analyze_mvt_variability(data)
      sex_comp <- result$sex_comparison

      expect_true("male_mean" %in% names(sex_comp))
      expect_true("female_mean" %in% names(sex_comp))
    })

    it("compares MVT by load", {
      analyzer <- AdvancedVelocityAnalyzer$new()
      data <- create_mvt_data()

      result <- analyzer$analyze_mvt_variability(data)
      load_comp <- result$load_comparison

      expect_true("load_80_mean" %in% names(load_comp))
      expect_true("load_90_mean" %in% names(load_comp))
      expect_true("wilcox_p" %in% names(load_comp))
    })

    it("throws error when no failure observations exist", {
      analyzer <- AdvancedVelocityAnalyzer$new()
      data <- create_mvt_data()
      data$rir <- 1  # No RIR=0 observations

      expect_error(
        analyzer$analyze_mvt_variability(data),
        "No failure observations"
      )
    })
  })

  # ============================================================================
  # H3: Day-to-Day Reliability Tests
  # ============================================================================

  describe("calculate_day_reliability (H3)", {

    it("returns ReliabilityResult object", {
      analyzer <- AdvancedVelocityAnalyzer$new()
      data <- create_reliability_data()

      result <- analyzer$calculate_day_reliability(data)

      expect_s3_class(result, "ReliabilityResult")
    })

    it("calculates ICC for slopes", {
      analyzer <- AdvancedVelocityAnalyzer$new()
      data <- create_reliability_data()

      result <- analyzer$calculate_day_reliability(data)
      slope_icc <- result$slope_icc

      expect_true("icc" %in% names(slope_icc))
      expect_true("ci_lower" %in% names(slope_icc))
      expect_true("ci_upper" %in% names(slope_icc))
      expect_true("sem" %in% names(slope_icc))
      expect_true("mdc95" %in% names(slope_icc))
      expect_true("interpretation" %in% names(slope_icc))

      # ICC should be between 0 and 1
      expect_true(slope_icc$icc >= 0 && slope_icc$icc <= 1)
    })

    it("calculates ICC for intercepts", {
      analyzer <- AdvancedVelocityAnalyzer$new()
      data <- create_reliability_data()

      result <- analyzer$calculate_day_reliability(data)

      expect_true(!is.na(result$intercept_icc$icc))
    })

    it("calculates ICC for MVT predictions", {
      analyzer <- AdvancedVelocityAnalyzer$new()
      data <- create_reliability_data()

      result <- analyzer$calculate_day_reliability(data)

      expect_true(!is.na(result$mvt_icc$icc))
    })

    it("returns day parameter estimates", {
      analyzer <- AdvancedVelocityAnalyzer$new()
      data <- create_reliability_data(n_participants = 5)

      result <- analyzer$calculate_day_reliability(data)
      day_params <- result$day_parameters

      expect_true("slope_day1" %in% names(day_params))
      expect_true("slope_day2" %in% names(day_params))
      expect_true("intercept_day1" %in% names(day_params))
      expect_true("intercept_day2" %in% names(day_params))
    })

    it("interprets ICC using Koo & Li thresholds", {
      analyzer <- AdvancedVelocityAnalyzer$new()
      data <- create_reliability_data()

      result <- analyzer$calculate_day_reliability(data)

      valid_interpretations <- c("Poor", "Moderate", "Good", "Excellent",
                                 "Insufficient data")
      expect_true(result$slope_icc$interpretation %in% valid_interpretations)
    })

    it("filters by load when specified", {
      analyzer <- AdvancedVelocityAnalyzer$new()
      data <- create_reliability_data()
      data$load_percentage <- sample(c("80%", "90%"), nrow(data), replace = TRUE)

      result <- analyzer$calculate_day_reliability(data, load = "90%")

      expect_s3_class(result, "ReliabilityResult")
    })

    it("throws error when no participants have both days", {
      analyzer <- AdvancedVelocityAnalyzer$new()
      data <- create_reliability_data()
      data$day <- "Day 1"  # Only Day 1 data

      expect_error(
        analyzer$calculate_day_reliability(data),
        "No participants with data on both days"
      )
    })
  })

  # ============================================================================
  # H4: Polynomial vs Linear Model Comparison Tests
  # ============================================================================

  describe("compare_polynomial_models (H4)", {

    it("returns ModelComparisonResult object", {
      skip_if_not_installed("lme4")

      analyzer <- AdvancedVelocityAnalyzer$new()
      data <- create_polynomial_data()

      result <- analyzer$compare_polynomial_models(data)

      expect_s3_class(result, "ModelComparisonResult")
    })

    it("compares individual participant models", {
      skip_if_not_installed("lme4")

      analyzer <- AdvancedVelocityAnalyzer$new()
      data <- create_polynomial_data(n_participants = 5)

      result <- analyzer$compare_polynomial_models(data)
      ind_results <- result$individual_results

      expect_true("id" %in% names(ind_results))
      expect_true("r2_adj_linear" %in% names(ind_results))
      expect_true("r2_adj_quad" %in% names(ind_results))
      expect_true("aic_linear" %in% names(ind_results))
      expect_true("aic_quad" %in% names(ind_results))
      expect_true("best_model" %in% names(ind_results))
    })

    it("computes delta AIC for model selection", {
      skip_if_not_installed("lme4")

      analyzer <- AdvancedVelocityAnalyzer$new()
      data <- create_polynomial_data()

      result <- analyzer$compare_polynomial_models(data)
      ind_results <- result$individual_results

      expect_true("delta_aic" %in% names(ind_results))
      # delta_aic = aic_quad - aic_linear
      expected_delta <- ind_results$aic_quad[1] - ind_results$aic_linear[1]
      expect_equal(ind_results$delta_aic[1], expected_delta, tolerance = 0.01)
    })

    it("includes population-level LMM comparison", {
      skip_if_not_installed("lme4")

      analyzer <- AdvancedVelocityAnalyzer$new()
      data <- create_polynomial_data()

      result <- analyzer$compare_polynomial_models(data)
      pop_comp <- result$population_comparison

      expect_true("aic_linear" %in% names(pop_comp))
      expect_true("aic_quad" %in% names(pop_comp))
      expect_true("lrt_p" %in% names(pop_comp))
    })

    it("provides model selection summary", {
      skip_if_not_installed("lme4")

      analyzer <- AdvancedVelocityAnalyzer$new()
      data <- create_polynomial_data()

      result <- analyzer$compare_polynomial_models(data)
      summary <- result$best_model_summary

      expect_true("n_participants" %in% names(summary))
      expect_true("n_linear_best" %in% names(summary))
      expect_true("n_quad_best" %in% names(summary))
      expect_true("recommendation" %in% names(summary))
    })
  })

  # ============================================================================
  # H5: Velocity Decay Tests
  # ============================================================================

  describe("analyze_velocity_decay (H5)", {

    it("returns VelocityDecayResult object", {
      analyzer <- AdvancedVelocityAnalyzer$new()
      data <- create_decay_data()

      result <- analyzer$analyze_velocity_decay(data)

      expect_s3_class(result, "VelocityDecayResult")
    })

    it("calculates set trajectories", {
      analyzer <- AdvancedVelocityAnalyzer$new()
      data <- create_decay_data(n_sets = 10)

      result <- analyzer$analyze_velocity_decay(data)
      trajectories <- result$set_trajectories

      expect_true("set_id" %in% names(trajectories))
      expect_true("rep_number" %in% names(trajectories))
      expect_true("delta_v" %in% names(trajectories))
      expect_true("cumulative_loss" %in% names(trajectories))
    })

    it("summarizes decay rates", {
      analyzer <- AdvancedVelocityAnalyzer$new()
      data <- create_decay_data()

      result <- analyzer$analyze_velocity_decay(data)
      decay_summary <- result$decay_summary

      expect_true("avg_decay_per_rep" %in% names(decay_summary))
      expect_true("sd_decay" %in% names(decay_summary))
      expect_true("early_reps_decay" %in% names(decay_summary))
      expect_true("late_reps_decay" %in% names(decay_summary))

      # Decay should be negative (velocity decreases)
      expect_true(decay_summary$avg_decay_per_rep < 0)
    })

    it("tests for decay acceleration", {
      analyzer <- AdvancedVelocityAnalyzer$new()
      data <- create_decay_data()

      result <- analyzer$analyze_velocity_decay(data)
      acceleration <- result$decay_acceleration

      expect_true("slope" %in% names(acceleration))
      expect_true("p_value" %in% names(acceleration))
      expect_true("accelerating" %in% names(acceleration))
      expect_true("interpretation" %in% names(acceleration))
    })

    it("detects decay breakpoint", {
      analyzer <- AdvancedVelocityAnalyzer$new()
      data <- create_decay_data()

      result <- analyzer$analyze_velocity_decay(data)
      breakpoint <- result$breakpoint

      expect_true("breakpoint_rep" %in% names(breakpoint))
      expect_true("avg_decay" %in% names(breakpoint))
      expect_true("interpretation" %in% names(breakpoint))
    })

    it("throws error when required columns missing", {
      analyzer <- AdvancedVelocityAnalyzer$new()
      data <- data.frame(set_id = "S1", rep_number = 1:5)  # Missing mean_velocity

      expect_error(
        analyzer$analyze_velocity_decay(data),
        "Missing required columns"
      )
    })
  })

  # ============================================================================
  # H6: Failure Prediction Tests
  # ============================================================================

  describe("build_failure_predictor (H6)", {

    it("returns FailurePredictionResult object", {
      analyzer <- AdvancedVelocityAnalyzer$new()
      data <- create_prediction_data()

      result <- analyzer$build_failure_predictor(data)

      expect_s3_class(result, "FailurePredictionResult")
    })

    it("fits multiple prediction models", {
      analyzer <- AdvancedVelocityAnalyzer$new()
      data <- create_prediction_data()

      result <- analyzer$build_failure_predictor(data)
      models <- result$prediction_models

      expect_true("v1_only" %in% names(models))
      expect_true("v1_load" %in% names(models))
    })

    it("performs leave-one-out cross-validation", {
      analyzer <- AdvancedVelocityAnalyzer$new()
      data <- create_prediction_data()

      result <- analyzer$build_failure_predictor(data)
      cv <- result$cv_results

      expect_true("mae" %in% names(cv))
      expect_true("rmse" %in% names(cv))
      expect_true("r2" %in% names(cv))
      expect_true("within_1_rep_pct" %in% names(cv))
      expect_true("within_2_reps_pct" %in% names(cv))
      expect_true("n_sets" %in% names(cv))

      # MAE should be reasonable (< 5 reps)
      expect_true(cv$mae < 5)
    })

    it("creates practical lookup table", {
      analyzer <- AdvancedVelocityAnalyzer$new()
      data <- create_prediction_data()

      result <- analyzer$build_failure_predictor(data)
      lookup <- result$lookup_table

      expect_true("first_rep_velocity" %in% names(lookup))
      expect_true("predicted_reps" %in% names(lookup))
      expect_true("lower_95" %in% names(lookup))
      expect_true("upper_95" %in% names(lookup))

      # Predictions should be at least 1 rep
      expect_true(all(lookup$predicted_reps >= 1))
    })

    it("higher v1 predicts more reps to failure", {
      analyzer <- AdvancedVelocityAnalyzer$new()
      data <- create_prediction_data()

      result <- analyzer$build_failure_predictor(data)
      lookup <- result$lookup_table

      # Correlation between v1 and predicted reps should be positive
      correlation <- cor(lookup$first_rep_velocity, lookup$predicted_reps)
      expect_true(correlation > 0)
    })

    it("throws error when required columns missing", {
      analyzer <- AdvancedVelocityAnalyzer$new()
      data <- data.frame(
        set_id = "S1",
        rep_number = 1:5,
        mean_velocity = seq(0.4, 0.2, length.out = 5)
      )  # Missing load_percentage, reps_to_failure

      expect_error(
        analyzer$build_failure_predictor(data),
        "Missing required columns"
      )
    })
  })

  # ============================================================================
  # Edge Cases and Robustness Tests
  # ============================================================================

  describe("edge cases", {

    it("handles minimum valid data for MVT", {
      analyzer <- AdvancedVelocityAnalyzer$new()
      data <- data.frame(
        id = c("P1", "P2", "P3"),
        sex = c("male", "male", "female"),
        load_percentage = c("80%", "90%", "80%"),
        rir = c(0, 0, 0),
        mean_velocity = c(0.18, 0.20, 0.22)
      )

      result <- analyzer$analyze_mvt_variability(data)

      expect_s3_class(result, "MvtAnalysisResult")
      expect_equal(result$population_stats$n, 3)
    })

    it("handles insufficient data for ICC gracefully", {
      analyzer <- AdvancedVelocityAnalyzer$new()

      # Data with only 2 participants - ICC needs at least 3
      data <- data.frame(
        id = rep(c("P1", "P2"), each = 16),
        day = rep(c("Day 1", "Day 2"), each = 8, times = 2),
        rir = rep(0:7, 4),
        load_percentage = "90%",
        mean_velocity = 0.2 + 0.04 * rep(0:7, 4) + rnorm(32, 0, 0.02)
      )

      result <- analyzer$calculate_day_reliability(data)

      # Should return result with NA ICCs or "Insufficient data" interpretation
      expect_s3_class(result, "ReliabilityResult")
    })

    it("handles single-set decay data", {
      analyzer <- AdvancedVelocityAnalyzer$new()
      data <- data.frame(
        set_id = rep("S1", 5),
        rep_number = 1:5,
        mean_velocity = c(0.40, 0.35, 0.30, 0.25, 0.20)
      )

      result <- analyzer$analyze_velocity_decay(data)

      expect_s3_class(result, "VelocityDecayResult")
      # With only one set, breakpoint detection should note insufficient data
      expect_true(
        is.na(result$breakpoint$breakpoint_rep) ||
          grepl("Insufficient", result$breakpoint$note %||% "")
      )
    })
  })
})
