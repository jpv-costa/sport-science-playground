# tests/testthat/test-model_plotter.R
# Tests for ModelPlotter R6 class
# Following CLAUDE.md testing principles: TDD, classicist approach, realistic data

box::use(
  testthat[describe, it, expect_true, expect_s3_class, expect_type,
           expect_named, expect_error, expect_gte, expect_lte],
  lme4[lmer]
)

# Set box path relative to project root
options(box.path = c("../../R", getOption("box.path")))

# Load the module under test
box::use(
  visualizers/model_plotter[ModelPlotter, MODEL_COLORS]
)

# Create realistic test data
create_test_data <- function(n_participants = 10, obs_per_participant = 5) {
  set.seed(42)
  participants <- paste0("P", seq_len(n_participants))

  do.call(rbind, lapply(participants, function(pid) {
    # Each participant has slightly different intercept and slope
    intercept <- 0.20 + rnorm(1, 0, 0.03)
    slope <- 0.030 + rnorm(1, 0, 0.005)

    rir <- sample(0:7, obs_per_participant, replace = TRUE)
    velocity <- intercept + slope * rir + rnorm(obs_per_participant, 0, 0.02)

    data.frame(
      id = pid,
      rir = rir,
      mean_velocity = velocity,
      stringsAsFactors = FALSE
    )
  }))
}

# Create a test LMM model
create_test_model <- function(data) {
  lmer(mean_velocity ~ rir + (1 + rir | id), data = data)
}

describe("ModelPlotter", {

  describe("initialization", {
    it("creates instance with default settings", {
      plotter <- ModelPlotter$new()
      expect_s3_class(plotter, "ModelPlotter")
    })

    it("accepts custom base size", {
      plotter <- ModelPlotter$new(base_size = 16)
      expect_s3_class(plotter, "ModelPlotter")
    })

    it("accepts custom colors", {
      custom_colors <- list(primary = "#000000", secondary = "#FFFFFF")
      plotter <- ModelPlotter$new(colors = custom_colors)
      expect_s3_class(plotter, "ModelPlotter")
    })
  })

  describe("extract_caterpillar_data", {
    data <- create_test_data()
    model <- create_test_model(data)
    plotter <- ModelPlotter$new()

    it("extracts random effects data from lmer model", {
      result <- plotter$extract_caterpillar_data(model)
      expect_s3_class(result, "data.frame")
    })

    it("returns data with required columns", {
      result <- plotter$extract_caterpillar_data(model)
      expect_true(all(c("id", "effect_type", "estimate", "lower", "upper") %in% names(result)))
    })

    it("returns both intercept and slope by default", {
      result <- plotter$extract_caterpillar_data(model)
      expect_true("intercept" %in% result$effect_type)
      expect_true("slope" %in% result$effect_type)
    })

    it("can extract only intercepts", {
      result <- plotter$extract_caterpillar_data(model, effect_type = "intercept")
      expect_true(all(result$effect_type == "intercept"))
    })

    it("can extract only slopes", {
      result <- plotter$extract_caterpillar_data(model, effect_type = "slope")
      expect_true(all(result$effect_type == "slope"))
    })

    it("errors for non-lmer models", {
      fake_model <- lm(mean_velocity ~ rir, data = data)
      expect_error(plotter$extract_caterpillar_data(fake_model))
    })
  })

  describe("plot_caterpillar_dual", {
    data <- create_test_data()
    model <- create_test_model(data)
    plotter <- ModelPlotter$new()

    it("returns a ggplot object", {
      p <- plotter$plot_caterpillar_dual(model)
      expect_s3_class(p, "ggplot")
    })

    it("accepts custom title", {
      p <- plotter$plot_caterpillar_dual(model, title = "Custom Title")
      expect_s3_class(p, "ggplot")
    })
  })

  describe("plot_model", {
    data <- create_test_data()
    model <- create_test_model(data)
    plotter <- ModelPlotter$new()

    it("returns a ggplot object", {
      p <- plotter$plot_model(data, model)
      expect_s3_class(p, "ggplot")
    })

    it("accepts custom title", {
      p <- plotter$plot_model(data, model, title = "Custom Model Plot")
      expect_s3_class(p, "ggplot")
    })

    it("errors when data missing required columns", {
      bad_data <- data.frame(x = 1:10, y = 1:10)
      expect_error(plotter$plot_model(bad_data, model))
    })

    it("errors for non-lmer models", {
      fake_model <- lm(mean_velocity ~ rir, data = data)
      expect_error(plotter$plot_model(data, fake_model))
    })
  })

  describe("plot_prediction_comparison", {
    data <- create_test_data()
    model <- create_test_model(data)
    plotter <- ModelPlotter$new()

    it("returns a ggplot object", {
      p <- plotter$plot_prediction_comparison(data, model)
      expect_s3_class(p, "ggplot")
    })

    it("accepts custom title", {
      p <- plotter$plot_prediction_comparison(data, model, title = "Custom Comparison")
      expect_s3_class(p, "ggplot")
    })
  })

  describe("calculate_prediction_errors", {
    data <- create_test_data()
    model <- create_test_model(data)
    plotter <- ModelPlotter$new()

    it("returns a list with required elements", {
      result <- plotter$calculate_prediction_errors(data, model)
      expect_type(result, "list")
      expect_named(result, c("pop_mae_mm", "ind_mae_mm", "pop_rmse_mm", "ind_rmse_mm", "improvement_pct"))
    })

    it("individual MAE is less than population MAE", {
      result <- plotter$calculate_prediction_errors(data, model)
      expect_lte(result$ind_mae_mm, result$pop_mae_mm)
    })

    it("improvement percentage is positive", {
      result <- plotter$calculate_prediction_errors(data, model)
      expect_gte(result$improvement_pct, 0)
    })

    it("returns values in mm/s", {
      result <- plotter$calculate_prediction_errors(data, model)
      # MAE should be in reasonable mm/s range (not m/s)
      expect_gte(result$pop_mae_mm, 1)  # At least 1 mm/s
      expect_lte(result$pop_mae_mm, 100)  # Less than 100 mm/s
    })
  })
})

describe("MODEL_COLORS", {
  it("is a named list", {
    expect_type(MODEL_COLORS, "list")
    expect_true(length(names(MODEL_COLORS)) > 0)
  })

  it("contains essential colors", {
    expect_true("primary" %in% names(MODEL_COLORS))
    expect_true("secondary" %in% names(MODEL_COLORS))
    expect_true("population" %in% names(MODEL_COLORS))
    expect_true("individual" %in% names(MODEL_COLORS))
  })
})
