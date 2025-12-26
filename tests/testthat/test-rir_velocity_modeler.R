# tests/testthat/test-rir_velocity_modeler.R
# Unit tests for RirVelocityModeler

box::use(
  testthat[...],
  ../../R/calculators/rir_velocity_modeler[
    RirVelocityModeler,
    RirVelocityModelResult,
    RirVelocityModelSummary
  ]
)

# ==============================================================================
# Test Fixtures
# ==============================================================================

create_synthetic_rir_data <- function(n = 50, r_squared_target = 0.8) {
  set.seed(42)
  # Create velocity data (typically 0.2 to 0.8 m/s)
  mean_velocity <- runif(n, 0.2, 0.8)

  # Create RIR with polynomial relationship + noise
  # Higher velocity = more RIR
  noise_sd <- sqrt((1 - r_squared_target) / r_squared_target) * 3
  rir <- 5 + 20 * mean_velocity - 10 * mean_velocity^2 + rnorm(n, 0, noise_sd)
  rir <- pmax(0, round(rir))  # RIR must be >= 0

  data.frame(
    mean_velocity = mean_velocity,
    rir = rir
  )
}

create_grouped_data <- function() {
  set.seed(42)
  data.frame(
    id = rep(1:5, each = 30),
    set_type = rep(c("RTF70", "RTF80", "RTF90"), times = 50),
    mean_velocity = runif(150, 0.2, 0.8),
    rir = sample(0:15, 150, replace = TRUE)
  )
}

# ==============================================================================
# RirVelocityModelResult Tests
# ==============================================================================

describe("RirVelocityModelResult", {

  describe("initialization", {
    it("stores model metadata", {
      result <- RirVelocityModelResult$new(
        model_type = "linear",
        r_squared = 0.75,
        rse = 2.5,
        coefficients = c(intercept = 5, slope = -10),
        n = 50,
        group_id = "participant_1"
      )

      expect_equal(result$model_type, "linear")
      expect_equal(result$r_squared, 0.75)
      expect_equal(result$rse, 2.5)
      expect_equal(result$n, 50)
      expect_equal(result$group_id, "participant_1")
    })
  })

  describe("to_list", {
    it("returns all fields as list", {
      result <- RirVelocityModelResult$new(
        model_type = "polynomial",
        r_squared = 0.85,
        rse = 1.8,
        coefficients = c(a = 1, b = 2, c = 3),
        n = 100
      )
      as_list <- result$to_list()

      expect_equal(as_list$model_type, "polynomial")
      expect_equal(as_list$r_squared, 0.85)
      expect_equal(as_list$rse, 1.8)
      expect_equal(as_list$n, 100)
    })
  })
})

# ==============================================================================
# RirVelocityModelSummary Tests
# ==============================================================================

describe("RirVelocityModelSummary", {

  describe("initialization", {
    it("calculates summary statistics from results", {
      results <- list(
        RirVelocityModelResult$new("poly", 0.80, 2.0, NA, 50),
        RirVelocityModelResult$new("poly", 0.85, 1.8, NA, 50),
        RirVelocityModelResult$new("poly", 0.90, 1.5, NA, 50)
      )

      summary <- RirVelocityModelSummary$new(results, "polynomial")

      expect_equal(summary$n_models, 3)
      expect_equal(summary$median_r_squared, 0.85)
      expect_equal(summary$mean_r_squared, 0.85)
      expect_equal(summary$min_r_squared, 0.80)
      expect_equal(summary$max_r_squared, 0.90)
    })

    it("handles NA values", {
      results <- list(
        RirVelocityModelResult$new("linear", 0.70, 2.5, NA, 50),
        RirVelocityModelResult$new("linear", NA_real_, NA_real_, NA, 5),
        RirVelocityModelResult$new("linear", 0.80, 2.0, NA, 50)
      )

      summary <- RirVelocityModelSummary$new(results, "linear")

      expect_equal(summary$n_models, 3)
      expect_equal(summary$median_r_squared, 0.75)
    })
  })

  describe("to_list", {
    it("returns summary as list", {
      results <- list(
        RirVelocityModelResult$new("poly", 0.85, 1.8, NA, 50)
      )

      summary <- RirVelocityModelSummary$new(results, "polynomial")
      as_list <- summary$to_list()

      expect_true("model_type" %in% names(as_list))
      expect_true("n_models" %in% names(as_list))
      expect_true("median_r_squared" %in% names(as_list))
    })
  })
})

# ==============================================================================
# RirVelocityModeler Tests
# ==============================================================================

describe("RirVelocityModeler", {

  describe("fit_linear", {
    it("fits linear model and returns result", {
      modeler <- RirVelocityModeler$new()
      data <- create_synthetic_rir_data(50)

      result <- modeler$fit_linear(data)

      expect_s3_class(result, "RirVelocityModelResult")
      expect_equal(result$model_type, "linear")
      expect_true(result$r_squared >= 0 && result$r_squared <= 1)
      expect_true(result$rse > 0)
      expect_equal(result$n, 50)
    })

    it("handles insufficient data", {
      modeler <- RirVelocityModeler$new()
      small_data <- data.frame(mean_velocity = c(0.5), rir = c(5))

      result <- modeler$fit_linear(small_data)

      expect_true(is.na(result$r_squared))
    })
  })

  describe("fit_polynomial", {
    it("fits polynomial model and returns result", {
      modeler <- RirVelocityModeler$new()
      data <- create_synthetic_rir_data(50)

      result <- modeler$fit_polynomial(data)

      expect_s3_class(result, "RirVelocityModelResult")
      expect_equal(result$model_type, "polynomial")
      expect_true(result$r_squared >= 0 && result$r_squared <= 1)
    })

    it("polynomial fits better than linear for curved data", {
      modeler <- RirVelocityModeler$new()

      # Create data with clear curvature
      set.seed(123)
      velocity <- seq(0.2, 0.8, length.out = 50)
      rir <- 10 + 30 * velocity - 25 * velocity^2 + rnorm(50, 0, 1)
      data <- data.frame(mean_velocity = velocity, rir = pmax(0, rir))

      linear_result <- modeler$fit_linear(data)
      poly_result <- modeler$fit_polynomial(data)

      expect_true(poly_result$r_squared >= linear_result$r_squared)
    })
  })

  describe("fit_both", {
    it("returns both linear and polynomial results", {
      modeler <- RirVelocityModeler$new()
      data <- create_synthetic_rir_data(50)

      results <- modeler$fit_both(data)

      expect_true("linear" %in% names(results))
      expect_true("polynomial" %in% names(results))
      expect_s3_class(results$linear, "RirVelocityModelResult")
      expect_s3_class(results$polynomial, "RirVelocityModelResult")
    })
  })

  describe("fit_general_by_load", {
    it("fits models for each load level", {
      modeler <- RirVelocityModeler$new()
      data <- create_grouped_data()

      results <- modeler$fit_general_by_load(data, "set_type")

      expect_true("RTF70" %in% names(results))
      expect_true("RTF80" %in% names(results))
      expect_true("RTF90" %in% names(results))
    })
  })

  describe("fit_individual", {
    it("fits models for each participant and load", {
      modeler <- RirVelocityModeler$new()
      data <- create_grouped_data()

      results <- modeler$fit_individual(data, "id", "set_type")

      expect_true("1" %in% names(results))
      expect_true(length(results) == 5)  # 5 participants
    })
  })

  describe("summarize_individual", {
    it("aggregates individual results", {
      modeler <- RirVelocityModeler$new()
      data <- create_grouped_data()

      individual_results <- modeler$fit_individual(data, "id", "set_type")
      summary <- modeler$summarize_individual(individual_results, "polynomial")

      expect_s3_class(summary, "RirVelocityModelSummary")
      expect_true(summary$n_models > 0)
    })
  })

  describe("calculate_prediction_accuracy", {
    it("computes prediction error", {
      modeler <- RirVelocityModeler$new()

      set.seed(42)
      train_data <- create_synthetic_rir_data(100, r_squared_target = 0.9)
      test_data <- create_synthetic_rir_data(50, r_squared_target = 0.9)

      accuracy <- modeler$calculate_prediction_accuracy(
        train_data, test_data, "polynomial"
      )

      expect_true("actual_rir" %in% names(accuracy))
      expect_true("predicted_rir" %in% names(accuracy))
      expect_true("error" %in% names(accuracy))
      expect_true("absolute_error" %in% names(accuracy))
      expect_equal(nrow(accuracy), 50)
    })
  })
})
