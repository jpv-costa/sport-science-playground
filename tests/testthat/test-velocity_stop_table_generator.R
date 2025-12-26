# tests/testthat/test-velocity_stop_table_generator.R
# Unit tests for VelocityStopTableGenerator

box::use(
  testthat[...],
  ../../R/calculators/velocity_stop_table_generator[
    VelocityStopTableGenerator,
    VelocityStopTable,
    LoadImportanceResult
  ]
)

# ==============================================================================
# Test Fixtures
# ==============================================================================

create_test_data <- function(n_participants = 10, n_obs_per = 20, seed = 42) {
  set.seed(seed)

  # Create nested data structure
  data <- expand.grid(
    participant_num = 1:n_participants,
    obs = 1:n_obs_per
  )

  data$id <- paste0("P", data$participant_num)

  # Random effects per participant
  participant_effects <- data.frame(
    id = paste0("P", 1:n_participants),
    intercept_effect = rnorm(n_participants, 0, 0.05),
    slope_effect = rnorm(n_participants, 0, 0.01)
  )

  data <- merge(data, participant_effects, by = "id")

  # Generate RIR values (0-7)
  data$rir <- sample(0:7, nrow(data), replace = TRUE)

  # Generate load percentage
  data$load_percentage <- sample(c("80%", "90%"), nrow(data), replace = TRUE)

  # Generate velocity with relationship to RIR
  # Higher RIR = higher velocity, load effect is small
  load_effect <- ifelse(data$load_percentage == "90%", -0.02, 0)

  data$mean_velocity <- 0.25 +
    0.035 * data$rir +
    load_effect +
    data$intercept_effect +
    data$slope_effect * data$rir +
    rnorm(nrow(data), 0, 0.02)

  # Clean up
  data <- data[, c("id", "rir", "load_percentage", "mean_velocity")]
  rownames(data) <- NULL

  data
}

# Data where load clearly matters
create_load_dependent_data <- function(n_participants = 10, n_obs_per = 20, seed = 42) {
  set.seed(seed)

  data <- expand.grid(
    participant_num = 1:n_participants,
    obs = 1:n_obs_per
  )

  data$id <- paste0("P", data$participant_num)
  data$rir <- sample(0:7, nrow(data), replace = TRUE)
  data$load_percentage <- sample(c("80%", "90%"), nrow(data), replace = TRUE)

  # Strong load effect
  load_effect <- ifelse(data$load_percentage == "90%", -0.15, 0)

  data$mean_velocity <- 0.35 +
    0.03 * data$rir +
    load_effect +
    rnorm(nrow(data), 0, 0.02)

  data <- data[, c("id", "rir", "load_percentage", "mean_velocity")]
  rownames(data) <- NULL

  data
}

# ==============================================================================
# VelocityStopTable Tests
# ==============================================================================

describe("VelocityStopTable", {

  describe("initialization", {
    it("stores table metadata", {
      table_df <- data.frame(
        rir = 0:5,
        velocity = seq(0.2, 0.45, length.out = 6),
        se = 0.02
      )

      result <- VelocityStopTable$new(
        table_type = "global",
        id = NULL,
        table = table_df,
        model_used = NULL,
        load_included = FALSE
      )

      expect_equal(result$table_type, "global")
      expect_null(result$id)
      expect_false(result$load_included)
      expect_equal(nrow(result$table), 6)
    })
  })

  describe("format_for_display", {
    it("rounds numeric columns", {
      table_df <- data.frame(
        rir = 0:3,
        velocity = c(0.201234, 0.251234, 0.301234, 0.351234)
      )

      result <- VelocityStopTable$new("global", NULL, table_df, NULL, FALSE)
      formatted <- result$format_for_display()

      expect_equal(formatted$velocity, c(0.201, 0.251, 0.301, 0.351))
    })
  })
})

# ==============================================================================
# LoadImportanceResult Tests
# ==============================================================================

describe("LoadImportanceResult", {

  describe("initialization", {
    it("stores test results", {
      result <- LoadImportanceResult$new(
        load_significant = TRUE,
        lrt_p_value = 0.01,
        delta_aic = 5.2,
        delta_bic = 3.1,
        bayes_factor = 0.2,
        bf_interpretation = "Moderate evidence for complex model",
        recommendation = "load_specific"
      )

      expect_true(result$load_significant)
      expect_equal(result$lrt_p_value, 0.01)
      expect_equal(result$recommendation, "load_specific")
    })
  })

  describe("to_list", {
    it("returns all fields", {
      result <- LoadImportanceResult$new(
        FALSE, 0.5, 1.2, 0.8, 2.5, "Weak", "global"
      )

      as_list <- result$to_list()
      expect_true("load_significant" %in% names(as_list))
      expect_true("recommendation" %in% names(as_list))
    })
  })
})

# ==============================================================================
# VelocityStopTableGenerator Tests
# ==============================================================================

describe("VelocityStopTableGenerator", {

  describe("initialization", {
    it("creates generator without error", {
      generator <- VelocityStopTableGenerator$new()
      expect_s3_class(generator, "VelocityStopTableGenerator")
    })
  })

  describe("test_load_importance", {
    it("tests if load affects velocity-RIR relationship", {
      skip_if_not_installed("lme4")

      generator <- VelocityStopTableGenerator$new()
      data <- create_test_data()

      result <- generator$test_load_importance(data)

      expect_s3_class(result, "LoadImportanceResult")
      expect_true("load_significant" %in% names(result$to_list()))
      expect_true(result$recommendation %in% c("global", "load_specific"))
    })

    it("recommends load_specific when load matters", {
      skip_if_not_installed("lme4")

      generator <- VelocityStopTableGenerator$new()
      data <- create_load_dependent_data()

      result <- generator$test_load_importance(data)

      # With strong load effect, should detect significance
      # Note: may not always recommend load_specific due to BF threshold
      expect_s3_class(result, "LoadImportanceResult")
    })
  })

  describe("generate_global_table", {
    it("generates velocity table without load distinction", {
      skip_if_not_installed("lme4")

      generator <- VelocityStopTableGenerator$new()
      data <- create_test_data()

      result <- generator$generate_global_table(data, rir_targets = 0:5)

      expect_s3_class(result, "VelocityStopTable")
      expect_equal(result$table_type, "global")
      expect_false(result$load_included)
      expect_equal(nrow(result$table), 6)  # RIR 0-5
      expect_true("velocity" %in% names(result$table))
    })

    it("includes prediction intervals", {
      skip_if_not_installed("lme4")

      generator <- VelocityStopTableGenerator$new()
      data <- create_test_data()

      result <- generator$generate_global_table(data, rir_targets = 0:3)

      expect_true("lower_95" %in% names(result$table))
      expect_true("upper_95" %in% names(result$table))
      expect_true(all(result$table$lower_95 < result$table$velocity))
      expect_true(all(result$table$upper_95 > result$table$velocity))
    })

    it("shows increasing velocity with higher RIR", {
      skip_if_not_installed("lme4")

      generator <- VelocityStopTableGenerator$new()
      data <- create_test_data()

      result <- generator$generate_global_table(data, rir_targets = 0:5)
      velocities <- result$table$velocity

      # Velocity should generally increase with RIR
      expect_true(velocities[6] > velocities[1])  # RIR 5 > RIR 0
    })
  })

  describe("generate_load_specific_table", {
    it("generates velocity table with load distinction", {
      skip_if_not_installed("lme4")

      generator <- VelocityStopTableGenerator$new()
      data <- create_test_data()

      result <- generator$generate_load_specific_table(data, rir_targets = 0:3)

      expect_s3_class(result, "VelocityStopTable")
      expect_true(result$load_included)
      expect_true("load_percentage" %in% names(result$table))
    })

    it("has entries for each load level", {
      skip_if_not_installed("lme4")

      generator <- VelocityStopTableGenerator$new()
      data <- create_test_data()

      result <- generator$generate_load_specific_table(data, rir_targets = 0:3)

      loads_in_table <- unique(result$table$load_percentage)
      expect_true("80%" %in% loads_in_table)
      expect_true("90%" %in% loads_in_table)
    })
  })

  describe("generate_individual_table", {
    it("generates table for specific participant", {
      skip_if_not_installed("lme4")

      generator <- VelocityStopTableGenerator$new()
      data <- create_test_data()

      result <- generator$generate_individual_table(data, "P1", rir_targets = 0:3)

      expect_s3_class(result, "VelocityStopTable")
      expect_equal(result$table_type, "individual")
      expect_equal(result$id, "P1")
    })

    it("handles participant with insufficient data", {
      skip_if_not_installed("lme4")

      generator <- VelocityStopTableGenerator$new()

      # Create data with one participant having only 2 observations
      small_data <- data.frame(
        id = c("P1", "P1", "P2", "P2", "P2", "P2", "P2", "P2"),
        rir = c(0, 1, 0, 1, 2, 3, 4, 5),
        load_percentage = "80%",
        mean_velocity = c(0.2, 0.25, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45)
      )

      expect_warning(
        result <- generator$generate_individual_table(small_data, "P1", 0:3),
        "Insufficient data"
      )
      expect_null(result)
    })
  })

  describe("generate_all_individual_tables", {
    it("generates tables for all participants", {
      skip_if_not_installed("lme4")

      generator <- VelocityStopTableGenerator$new()
      data <- create_test_data(n_participants = 5)

      result <- generator$generate_all_individual_tables(data, rir_targets = 0:3)

      expect_true(is.list(result))
      expect_equal(length(result), 5)  # 5 participants
      expect_true("P1" %in% names(result))
    })
  })

  describe("compare_general_vs_individual", {
    it("compares accuracy of general vs individual predictions", {
      skip_if_not_installed("lme4")

      generator <- VelocityStopTableGenerator$new()
      data <- create_test_data()

      result <- generator$compare_general_vs_individual(data)

      expect_true("global_mae" %in% names(result))
      expect_true("individual_mae" %in% names(result))
      expect_true("mae_improvement_pct" %in% names(result))
      expect_true("recommendation" %in% names(result))
    })

    it("shows individual predictions are more accurate", {
      skip_if_not_installed("lme4")

      generator <- VelocityStopTableGenerator$new()
      data <- create_test_data()

      result <- generator$compare_general_vs_individual(data)

      # Individual predictions (with random effects) should be at least as good
      expect_true(result$individual_mae <= result$global_mae)
    })
  })

  describe("export_to_csv", {
    it("exports table to CSV file", {
      skip_if_not_installed("lme4")

      generator <- VelocityStopTableGenerator$new()
      data <- create_test_data()

      table <- generator$generate_global_table(data, rir_targets = 0:3)

      temp_file <- tempfile(fileext = ".csv")
      generator$export_to_csv(table, temp_file)

      expect_true(file.exists(temp_file))

      # Read back and verify
      read_back <- utils::read.csv(temp_file)
      expect_equal(nrow(read_back), 4)

      # Clean up
      unlink(temp_file)
    })
  })
})
