# tests/testthat/test-diagnostic_plotter.R
# Unit tests for DiagnosticPlotter and related classes
#
# =============================================================================
# TESTING PRINCIPLES (from CLAUDE.md)
# =============================================================================
# - Classicist approach: Real objects, no mocks, verify state through public API
# - Property-based where applicable: Verify invariants across inputs
# - Realistic fake data: Mimic actual data distributions
# =============================================================================

box::use(
  testthat[describe, it, expect_equal, expect_true, expect_false,
           expect_s3_class, expect_error, expect_type, expect_length],
  ../../R/visualizers/diagnostic_plotter[
    BlandAltmanResult,
    SpaghettiDataGenerator,
    DiagnosticPlotter,
    PLOT_COLORS
  ]
)

# =============================================================================
# TEST FIXTURES
# =============================================================================

create_test_reliability_data <- function(n = 20, seed = 42) {
  set.seed(seed)
  day1 <- rnorm(n, mean = 0.03, sd = 0.01)
  day2 <- day1 + rnorm(n, mean = 0.002, sd = 0.005)  # Small systematic bias
  list(day1 = day1, day2 = day2)
}

create_test_velocity_data <- function(n_subjects = 10, n_obs_per_subject = 8, seed = 42) {
  set.seed(seed)

  data <- do.call(rbind, lapply(1:n_subjects, function(i) {
    # Individual intercept and slope
    intercept <- rnorm(1, mean = 0.2, sd = 0.03)
    slope <- rnorm(1, mean = 0.03, sd = 0.008)

    rir <- 0:(n_obs_per_subject - 1)
    velocity <- intercept + slope * rir + rnorm(n_obs_per_subject, 0, 0.02)

    data.frame(
      id = paste0("P", sprintf("%02d", i)),
      rir = rir,
      mean_velocity = velocity
    )
  }))

  data
}

create_test_random_effects <- function(n = 15, seed = 42) {
  set.seed(seed)
  estimates <- rnorm(n, mean = 0, sd = 0.01)
  se <- runif(n, 0.002, 0.005)

  data.frame(
    id = paste0("P", sprintf("%02d", 1:n)),
    estimate = estimates,
    lower = estimates - 1.96 * se,
    upper = estimates + 1.96 * se
  )
}

create_test_bootstrap_samples <- function(n = 1000, observed = 0.029, seed = 42) {
  set.seed(seed)
  rnorm(n, mean = observed, sd = 0.004)
}

# =============================================================================
# BLAND-ALTMAN RESULT TESTS
# =============================================================================

describe("BlandAltmanResult", {

  describe("initialization", {

    it("calculates correct statistics for known values", {
      # Known values for easy verification
      day1 <- c(1, 2, 3, 4, 5)
      day2 <- c(1.1, 2.2, 2.9, 4.1, 4.8)

      result <- BlandAltmanResult$new(day1, day2)

      # Expected: differences = c(-0.1, -0.2, 0.1, -0.1, 0.2)
      # mean_diff = -0.02
      expect_equal(result$mean_diff, -0.02, tolerance = 1e-10)
      expect_equal(result$n, 5)
    })

    it("calculates correct limits of agreement", {
      data <- create_test_reliability_data(n = 100)
      result <- BlandAltmanResult$new(data$day1, data$day2)

      # LOA should be symmetric around mean_diff
      expected_width <- 2 * 1.96 * result$sd_diff
      actual_width <- result$upper_loa - result$lower_loa

      expect_equal(actual_width, expected_width, tolerance = 1e-10)
    })

    it("rejects non-numeric inputs", {
      expect_error(
        BlandAltmanResult$new(c("a", "b"), c(1, 2)),
        "day1 must be numeric"
      )
    })

    it("rejects mismatched lengths", {
      expect_error(
        BlandAltmanResult$new(c(1, 2, 3), c(1, 2)),
        "day1 and day2 must have same length"
      )
    })

    it("rejects too few observations", {
      expect_error(
        BlandAltmanResult$new(c(1, 2), c(1, 2)),
        "Need at least 3 observations"
      )
    })
  })

  describe("is_bias_acceptable", {

    it("returns TRUE when bias below threshold", {
      day1 <- c(1, 2, 3, 4, 5)
      day2 <- c(1.01, 2.01, 3.01, 4.01, 5.01)  # Small consistent bias

      result <- BlandAltmanResult$new(day1, day2)

      expect_true(result$is_bias_acceptable(threshold = 0.05))
    })

    it("returns FALSE when bias above threshold", {
      day1 <- c(1, 2, 3, 4, 5)
      day2 <- c(0.5, 1.5, 2.5, 3.5, 4.5)  # Large consistent bias

      result <- BlandAltmanResult$new(day1, day2)

      expect_false(result$is_bias_acceptable(threshold = 0.1))
    })
  })

  describe("get_loa_width", {

    it("returns correct width", {
      data <- create_test_reliability_data()
      result <- BlandAltmanResult$new(data$day1, data$day2)

      expected_width <- result$upper_loa - result$lower_loa

      expect_equal(result$get_loa_width(), expected_width)
    })
  })

  describe("to_list", {

    it("returns all fields", {
      data <- create_test_reliability_data()
      result <- BlandAltmanResult$new(data$day1, data$day2)

      as_list <- result$to_list()

      expect_type(as_list, "list")
      expect_true("mean_diff" %in% names(as_list))
      expect_true("sd_diff" %in% names(as_list))
      expect_true("upper_loa" %in% names(as_list))
      expect_true("lower_loa" %in% names(as_list))
      expect_true("n" %in% names(as_list))
      expect_true("loa_width" %in% names(as_list))
    })
  })
})

# =============================================================================
# SPAGHETTI DATA GENERATOR TESTS
# =============================================================================

describe("SpaghettiDataGenerator", {

  describe("generate", {

    it("returns correct structure", {
      data <- create_test_velocity_data()
      generator <- SpaghettiDataGenerator$new()

      result <- generator$generate(data)

      expect_type(result, "list")
      expect_true("individual_lines" %in% names(result))
      expect_true("population_line" %in% names(result))
      expect_true("n_individuals" %in% names(result))
      expect_true("rir_range" %in% names(result))
    })

    it("generates correct number of individuals", {
      data <- create_test_velocity_data(n_subjects = 15)
      generator <- SpaghettiDataGenerator$new()

      result <- generator$generate(data)

      expect_equal(result$n_individuals, 15)
    })

    it("generates individual lines for all subjects", {
      data <- create_test_velocity_data(n_subjects = 10, n_obs_per_subject = 8)
      generator <- SpaghettiDataGenerator$new()

      result <- generator$generate(data, n_points = 20)

      # Should have n_subjects * n_points rows
      n_unique_ids <- length(unique(result$individual_lines$id))
      expect_equal(n_unique_ids, 10)
    })

    it("respects n_points parameter", {
      data <- create_test_velocity_data(n_subjects = 5)
      generator <- SpaghettiDataGenerator$new()

      result <- generator$generate(data, n_points = 100)

      expect_equal(nrow(result$population_line), 100)
    })

    it("rejects data without required columns", {
      bad_data <- data.frame(x = 1:10, y = 1:10)
      generator <- SpaghettiDataGenerator$new()

      expect_error(generator$generate(bad_data), "data must have 'id' column")
    })
  })
})

# =============================================================================
# DIAGNOSTIC PLOTTER TESTS
# =============================================================================

describe("DiagnosticPlotter", {

  describe("initialization", {

    it("uses default colors when none provided", {
      plotter <- DiagnosticPlotter$new()
      # Can't directly access private, but can verify plots work
      expect_s3_class(plotter, "DiagnosticPlotter")
    })

    it("accepts custom colors", {
      custom_colors <- list(
        primary = "#000000",
        secondary = "#FFFFFF",
        tertiary = "#FF0000",
        gray = "#888888"
      )

      plotter <- DiagnosticPlotter$new(colors = custom_colors)
      expect_s3_class(plotter, "DiagnosticPlotter")
    })
  })

  describe("plot_spaghetti", {

    it("returns ggplot object", {
      data <- create_test_velocity_data()
      plotter <- DiagnosticPlotter$new()

      result <- plotter$plot_spaghetti(data)

      expect_s3_class(result, "ggplot")
    })

    it("works with highlight_population = FALSE", {
      data <- create_test_velocity_data()
      plotter <- DiagnosticPlotter$new()

      result <- plotter$plot_spaghetti(data, highlight_population = FALSE)

      expect_s3_class(result, "ggplot")
    })
  })

  describe("plot_diagnostics_panel", {

    it("returns patchwork object", {
      residuals <- rnorm(100)
      fitted <- rnorm(100, mean = 0.3)
      plotter <- DiagnosticPlotter$new()

      result <- plotter$plot_diagnostics_panel(residuals, fitted)

      expect_s3_class(result, "patchwork")
    })

    it("rejects mismatched lengths", {
      plotter <- DiagnosticPlotter$new()

      expect_error(
        plotter$plot_diagnostics_panel(rnorm(10), rnorm(20)),
        "residuals and fitted must have same length"
      )
    })
  })

  describe("plot_bland_altman", {

    it("returns ggplot object", {
      data <- create_test_reliability_data()
      plotter <- DiagnosticPlotter$new()

      result <- plotter$plot_bland_altman(data$day1, data$day2)

      expect_s3_class(result, "ggplot")
    })

    it("accepts custom title and label", {
      data <- create_test_reliability_data()
      plotter <- DiagnosticPlotter$new()

      result <- plotter$plot_bland_altman(
        data$day1, data$day2,
        title = "Custom Title",
        y_label = "Custom Label"
      )

      expect_s3_class(result, "ggplot")
    })
  })

  describe("plot_caterpillar", {

    it("returns ggplot object", {
      random_effects <- create_test_random_effects()
      plotter <- DiagnosticPlotter$new()

      result <- plotter$plot_caterpillar(random_effects)

      expect_s3_class(result, "ggplot")
    })

    it("rejects data without required columns", {
      bad_data <- data.frame(x = 1:10, y = 1:10)
      plotter <- DiagnosticPlotter$new()

      expect_error(
        plotter$plot_caterpillar(bad_data),
        "random_effects must have 'id' column"
      )
    })
  })

  describe("plot_bootstrap_distribution", {

    it("returns ggplot object", {
      samples <- create_test_bootstrap_samples()
      plotter <- DiagnosticPlotter$new()

      result <- plotter$plot_bootstrap_distribution(
        bootstrap_samples = samples,
        observed_value = 0.029,
        ci_lower = 0.022,
        ci_upper = 0.036
      )

      expect_s3_class(result, "ggplot")
    })

    it("rejects too few bootstrap samples", {
      plotter <- DiagnosticPlotter$new()

      expect_error(
        plotter$plot_bootstrap_distribution(
          bootstrap_samples = rnorm(50),
          observed_value = 0.029,
          ci_lower = 0.022,
          ci_upper = 0.036
        ),
        "Need at least 100 bootstrap samples"
      )
    })

    it("rejects invalid CI bounds", {
      plotter <- DiagnosticPlotter$new()

      expect_error(
        plotter$plot_bootstrap_distribution(
          bootstrap_samples = rnorm(1000),
          observed_value = 0.029,
          ci_lower = 0.036,  # Wrong order
          ci_upper = 0.022
        ),
        "ci_lower must be less than ci_upper"
      )
    })
  })
})

# =============================================================================
# VELOCITY ZONES TESTS
# =============================================================================

create_test_velocity_thresholds <- function() {
  data.frame(
    rir = 0:7,
    velocity = c(0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50),
    lower = c(0.12, 0.17, 0.22, 0.27, 0.32, 0.37, 0.42, 0.47),
    upper = c(0.18, 0.23, 0.28, 0.33, 0.38, 0.43, 0.48, 0.53)
  )
}

describe("DiagnosticPlotter", {

  describe("plot_velocity_zones", {

    it("returns ggplot object", {
      thresholds <- create_test_velocity_thresholds()
      plotter <- DiagnosticPlotter$new()

      result <- plotter$plot_velocity_zones(thresholds)

      expect_s3_class(result, "ggplot")
    })

    it("works without CI columns", {
      thresholds <- data.frame(
        rir = 0:5,
        velocity = c(0.15, 0.20, 0.25, 0.30, 0.35, 0.40)
      )
      plotter <- DiagnosticPlotter$new()

      result <- plotter$plot_velocity_zones(thresholds)

      expect_s3_class(result, "ggplot")
    })

    it("rejects data without rir column", {
      bad_data <- data.frame(x = 1:5, velocity = 1:5)
      plotter <- DiagnosticPlotter$new()

      expect_error(
        plotter$plot_velocity_zones(bad_data),
        "velocity_thresholds must have 'rir' column"
      )
    })
  })

  describe("plot_coverage_calibration", {

    it("returns ggplot object", {
      calibration <- data.frame(
        target_coverage = seq(0.5, 0.95, by = 0.05),
        empirical_coverage = seq(0.48, 0.93, by = 0.05)
      )
      plotter <- DiagnosticPlotter$new()

      result <- plotter$plot_coverage_calibration(calibration)

      expect_s3_class(result, "ggplot")
    })

    it("rejects data without required columns", {
      bad_data <- data.frame(x = 1:10, y = 1:10)
      plotter <- DiagnosticPlotter$new()

      expect_error(
        plotter$plot_coverage_calibration(bad_data),
        "calibration_data must have 'target_coverage' column"
      )
    })
  })

  describe("plot_exercise_comparison", {

    it("returns ggplot object", {
      comparison <- data.frame(
        exercise = c("Deadlift", "Squat"),
        slope = c(0.029, 0.037),
        slope_se = c(0.003, 0.004)
      )
      plotter <- DiagnosticPlotter$new()

      result <- plotter$plot_exercise_comparison(comparison)

      expect_s3_class(result, "ggplot")
    })

    it("rejects data without exercise column", {
      bad_data <- data.frame(slope = 0.03, slope_se = 0.003)
      plotter <- DiagnosticPlotter$new()

      expect_error(
        plotter$plot_exercise_comparison(bad_data),
        "comparison_data must have 'exercise' column"
      )
    })
  })
})

# =============================================================================
# SCIENTIFIC VALIDITY TESTS
# =============================================================================

describe("Scientific Validity", {

  describe("Bland-Altman LOA calculation", {

    it("uses correct 1.96 multiplier for 95% LOA", {
      # With normal data, ~95% of differences should fall within LOA
      set.seed(42)
      n <- 1000
      day1 <- rnorm(n, mean = 0.3, sd = 0.05)
      day2 <- day1 + rnorm(n, mean = 0, sd = 0.02)  # No systematic bias

      result <- BlandAltmanResult$new(day1, day2)
      differences <- day1 - day2

      # Count how many fall within LOA
      within_loa <- sum(
        differences >= result$lower_loa & differences <= result$upper_loa
      )
      coverage <- within_loa / n

      # Should be approximately 95% (allow some sampling variability)
      expect_true(coverage > 0.93 && coverage < 0.97)
    })
  })
})
