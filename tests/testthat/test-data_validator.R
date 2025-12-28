# tests/testthat/test-data_validator.R
# Unit tests for DataValidator and ValidationResult

box::use(
  testthat[...],
  ../../R/validators/run_validation[
    DataValidator,
    ValidationResult
  ]
)

# ==============================================================================
# Test Fixtures
# ==============================================================================

#' Create valid test data matching DeadliftRirDataLoader output
create_valid_data <- function(n_rows = 50, seed = 42) {
  set.seed(seed)

  data.frame(
    id = rep(paste0("P", 1:5), each = n_rows / 5),
    sex = rep(c("male", "male", "female", "male", "female"), each = n_rows / 5),
    day = sample(c("Day 1", "Day 2"), n_rows, replace = TRUE),
    load_percentage = sample(c("80%", "90%"), n_rows, replace = TRUE),
    rir = sample(0:7, n_rows, replace = TRUE),
    mean_velocity = runif(n_rows, 0.15, 0.55),
    set_id = paste0("S", seq_len(n_rows)),
    rep_number = rep(1:10, length.out = n_rows),
    stringsAsFactors = FALSE
  )
}

# ==============================================================================
# ValidationResult Tests
# ==============================================================================

describe("ValidationResult", {

  describe("initialization", {
    it("stores all fields correctly", {
      result <- ValidationResult$new(
        valid = TRUE,
        errors = character(0),
        warnings = c("Minor issue"),
        checks = list(test = list(valid = TRUE)),
        n_rows = 100
      )

      expect_true(result$valid)
      expect_length(result$errors, 0)
      expect_length(result$warnings, 1)
      expect_equal(result$n_rows, 100)
    })
  })

  describe("summary", {
    it("returns PASSED for valid data", {
      result <- ValidationResult$new(
        valid = TRUE,
        errors = character(0),
        warnings = character(0),
        checks = list(),
        n_rows = 50
      )

      summary <- result$summary()
      expect_true(grepl("PASSED", summary))
      expect_true(grepl("50", summary))
    })

    it("returns FAILED for invalid data", {
      result <- ValidationResult$new(
        valid = FALSE,
        errors = c("Error 1", "Error 2"),
        warnings = character(0),
        checks = list(),
        n_rows = 50
      )

      summary <- result$summary()
      expect_true(grepl("FAILED", summary))
      expect_true(grepl("Errors: 2", summary))
    })
  })
})

# ==============================================================================
# DataValidator Tests
# ==============================================================================

describe("DataValidator", {

  describe("initialization", {
    it("creates validator without error", {
      validator <- DataValidator$new()
      expect_s3_class(validator, "DataValidator")
    })
  })

  describe("validate", {

    it("returns ValidationResult for valid data", {
      validator <- DataValidator$new()
      data <- create_valid_data()

      result <- validator$validate(data)

      expect_s3_class(result, "ValidationResult")
      expect_true(result$valid)
      expect_length(result$errors, 0)
    })

    it("detects missing required columns", {
      validator <- DataValidator$new()
      data <- data.frame(id = "P1", rir = 0)  # Missing most columns

      result <- validator$validate(data)

      expect_false(result$valid)
      expect_true(any(grepl("Missing required columns", result$errors)))
    })

    it("detects RIR out of range", {
      validator <- DataValidator$new()
      data <- create_valid_data()
      data$rir[1] <- 10  # Invalid RIR

      result <- validator$validate(data)

      expect_false(result$valid)
      expect_true(any(grepl("RIR out of range", result$errors)))
    })

    it("detects negative RIR values", {
      validator <- DataValidator$new()
      data <- create_valid_data()
      data$rir[1] <- -1  # Invalid RIR

      result <- validator$validate(data)

      expect_false(result$valid)
      expect_true(any(grepl("RIR out of range", result$errors)))
    })

    it("detects non-positive velocity", {
      validator <- DataValidator$new()
      data <- create_valid_data()
      data$mean_velocity[1] <- -0.1  # Invalid velocity

      result <- validator$validate(data)

      expect_false(result$valid)
      expect_true(any(grepl("Non-positive velocity", result$errors)))
    })

    it("detects zero velocity", {
      validator <- DataValidator$new()
      data <- create_valid_data()
      data$mean_velocity[1] <- 0  # Invalid velocity

      result <- validator$validate(data)

      expect_false(result$valid)
    })

    it("detects invalid load percentages", {
      validator <- DataValidator$new()
      data <- create_valid_data()
      data$load_percentage[1] <- "75%"  # Invalid load

      result <- validator$validate(data)

      expect_false(result$valid)
      expect_true(any(grepl("Invalid load percentages", result$errors)))
    })

    it("detects invalid day values", {
      validator <- DataValidator$new()
      data <- create_valid_data()
      data$day[1] <- "Day 3"  # Invalid day

      result <- validator$validate(data)

      expect_false(result$valid)
      expect_true(any(grepl("Invalid day values", result$errors)))
    })

    it("detects empty participant IDs", {
      validator <- DataValidator$new()
      data <- create_valid_data()
      data$id[1] <- ""  # Empty ID

      result <- validator$validate(data)

      # This should be a warning, not error
      expect_true(length(result$warnings) > 0 || !result$valid)
    })

    it("warns about velocity outliers", {
      validator <- DataValidator$new()
      data <- create_valid_data()
      data$mean_velocity[1] <- 1.5  # Outlier velocity

      result <- validator$validate(data)

      # Should pass but with warning
      expect_true(result$valid)
      expect_true(any(grepl("outliers", result$warnings)))
    })

    it("strict mode treats warnings as errors", {
      validator <- DataValidator$new()
      data <- create_valid_data()
      data$mean_velocity[1] <- 1.5  # Outlier (warning)

      result <- validator$validate(data, strict = TRUE)

      expect_false(result$valid)
    })
  })

  describe("is_valid", {

    it("returns TRUE for valid data", {
      validator <- DataValidator$new()
      data <- create_valid_data()

      expect_true(validator$is_valid(data))
    })

    it("returns FALSE for invalid data", {
      validator <- DataValidator$new()
      data <- create_valid_data()
      data$rir[1] <- 99

      expect_false(validator$is_valid(data))
    })
  })

  describe("edge cases", {

    it("handles empty data frame", {
      validator <- DataValidator$new()
      data <- data.frame()

      result <- validator$validate(data)

      expect_false(result$valid)
      expect_true(any(grepl("Missing required columns", result$errors)))
    })

    it("handles single row of valid data", {
      validator <- DataValidator$new()
      data <- data.frame(
        id = "P1",
        sex = "male",
        day = "Day 1",
        load_percentage = "80%",
        rir = 3,
        mean_velocity = 0.35,
        set_id = "S1",
        rep_number = 1
      )

      result <- validator$validate(data)

      expect_true(result$valid)
    })

    it("handles NA values in RIR", {
      validator <- DataValidator$new()
      data <- create_valid_data()
      data$rir[1] <- NA

      result <- validator$validate(data)

      expect_false(result$valid)
    })

    it("handles NA values in velocity", {
      validator <- DataValidator$new()
      data <- create_valid_data()
      data$mean_velocity[1] <- NA

      result <- validator$validate(data)

      expect_false(result$valid)
    })
  })

  describe("validation checks object", {

    it("includes all check results", {
      validator <- DataValidator$new()
      data <- create_valid_data()

      result <- validator$validate(data)

      expect_true("required_columns" %in% names(result$checks))
      expect_true("rir_range" %in% names(result$checks))
      expect_true("velocity_positive" %in% names(result$checks))
      expect_true("load_values" %in% names(result$checks))
      expect_true("day_values" %in% names(result$checks))
    })

    it("each check has valid flag and message", {
      validator <- DataValidator$new()
      data <- create_valid_data()

      result <- validator$validate(data)

      for (check_name in names(result$checks)) {
        check <- result$checks[[check_name]]
        expect_true("valid" %in% names(check))
        expect_true("message" %in% names(check))
      }
    })
  })
})
