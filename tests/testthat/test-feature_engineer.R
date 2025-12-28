# tests/testthat/test-feature_engineer.R
# Unit tests for FeatureEngineer (comprehensive feature engineering)

box::use(
  testthat[...],
  ../../R/calculators/feature_engineer[
    FeatureEngineer, FeatureEngineeringResult
  ]
)

# Helper function to create test data matching the loader output
create_test_data <- function(n_participants = 5, n_sets_per_participant = 4) {
  set.seed(42)

  data_list <- lapply(seq_len(n_participants), function(p) {
    participant_id <- sprintf("P%02d", p)
    base_velocity <- runif(1, 0.4, 0.6)

    sets <- lapply(seq_len(n_sets_per_participant), function(s) {
      load <- sample(c("80%", "90%"), 1)
      day <- sample(c("Day 1", "Day 2"), 1)
      reps <- sample(5:10, 1)

      data.frame(
        id = participant_id,
        sex = sample(c("male", "female"), 1),
        day = day,
        load_percentage = load,
        weight_kg = ifelse(load == "80%", 120, 140),
        rir = rev(seq(0, reps - 1)),  # RIR decreases from high to 0
        mean_velocity = base_velocity - seq(0, reps - 1) * 0.02 + rnorm(reps, 0, 0.02),
        set_id = paste(participant_id, gsub(" ", "", day),
                       gsub("%", "pct", load), paste0("S", s), sep = "_"),
        rep_number = seq_len(reps),
        reps_to_failure = reps,
        stringsAsFactors = FALSE
      )
    })

    do.call(rbind, sets)
  })

  do.call(rbind, data_list)
}

describe("FeatureEngineeringResult", {

  it("stores all feature levels", {
    obs_df <- data.frame(id = "P01", velocity_z_global = 0.5)
    set_df <- data.frame(set_id = "S1", set_mean_velocity = 0.5)
    part_df <- data.frame(id = "P01", p_velocity_cv = 0.1)
    ps_df <- data.frame(id = "P01", set_id = "S1", ps_velocity_deviation = 0.01)
    pl_df <- data.frame(id = "P01", load_percentage = "80%", pl_mean_velocity = 0.5)
    pr_df <- data.frame(id = "P01", rir = 0, pr_mean_velocity = 0.4)
    names_list <- list(observation = "velocity_z_global", set = "set_mean_velocity",
                       participant = "p_velocity_cv")
    global_stats <- list(global_mean_velocity = 0.45)

    result <- FeatureEngineeringResult$new(
      obs_df, set_df, part_df, ps_df, pl_df, pr_df, names_list, global_stats
    )

    expect_equal(nrow(result$observation_features), 1)
    expect_equal(nrow(result$set_features), 1)
    expect_equal(nrow(result$participant_features), 1)
    expect_equal(nrow(result$participant_set_features), 1)
    expect_equal(nrow(result$participant_load_features), 1)
    expect_equal(nrow(result$participant_rir_features), 1)
  })

  it("provides summary", {
    obs_df <- data.frame(id = c("P01", "P02"), velocity_z_global = c(0.5, -0.3))
    set_df <- data.frame(set_id = c("S1", "S2"), set_mean_velocity = c(0.5, 0.45))
    part_df <- data.frame(id = c("P01", "P02"), p_velocity_cv = c(0.1, 0.15))
    ps_df <- data.frame(id = "P01", set_id = "S1", ps_velocity_deviation = 0.01)
    pl_df <- data.frame(id = "P01", load_percentage = "80%", pl_mean_velocity = 0.5)
    pr_df <- data.frame(id = "P01", rir = 0, pr_mean_velocity = 0.4)
    names_list <- list(
      observation = "velocity_z_global",
      set = "set_mean_velocity",
      participant = "p_velocity_cv"
    )
    global_stats <- list(global_mean_velocity = 0.45)

    result <- FeatureEngineeringResult$new(
      obs_df, set_df, part_df, ps_df, pl_df, pr_df, names_list, global_stats
    )
    summary <- result$summarize()

    expect_equal(summary$n_observations, 2)
    expect_equal(summary$n_sets, 2)
    expect_equal(summary$n_participants, 2)
  })

  it("extracts observation matrix", {
    fe <- FeatureEngineer$new()
    data <- create_test_data()

    result <- fe$engineer_features(data)
    matrix <- result$get_observation_matrix()

    expect_true("velocity_z_global" %in% names(matrix))
    expect_true(is.data.frame(matrix))
  })

  it("extracts participant matrix", {
    fe <- FeatureEngineer$new()
    data <- create_test_data()

    result <- fe$engineer_features(data)
    matrix <- result$get_participant_matrix()

    expect_true("p_mean_velocity" %in% names(matrix) || "p_velocity_cv" %in% names(matrix))
    expect_true(is.data.frame(matrix))
  })
})

describe("FeatureEngineer", {

  describe("engineer_features", {

    it("validates required columns", {
      fe <- FeatureEngineer$new()

      # Missing columns
      bad_data <- data.frame(id = "P01", velocity = 0.5)

      expect_error(
        fe$engineer_features(bad_data),
        "Missing required columns"
      )
    })

    it("returns FeatureEngineeringResult object", {
      fe <- FeatureEngineer$new()
      data <- create_test_data()

      result <- fe$engineer_features(data)

      expect_s3_class(result, "FeatureEngineeringResult")
    })

    it("creates observation features", {
      fe <- FeatureEngineer$new()
      data <- create_test_data()

      result <- fe$engineer_features(data)
      obs <- result$observation_features

      expect_true("velocity_z_global" %in% names(obs))
      expect_true("velocity_z_within_participant" %in% names(obs))
      expect_true("set_position_normalized" %in% names(obs))
      expect_true("velocity_z_within_set" %in% names(obs))
    })

    it("creates set features", {
      fe <- FeatureEngineer$new()
      data <- create_test_data()

      result <- fe$engineer_features(data)
      sets <- result$set_features

      expect_true("set_mean_velocity" %in% names(sets))
      expect_true("set_velocity_sd" %in% names(sets))
      expect_true("set_velocity_cv" %in% names(sets))
      expect_true("set_decay_slope" %in% names(sets))
    })

    it("creates participant features", {
      fe <- FeatureEngineer$new()
      data <- create_test_data()

      result <- fe$engineer_features(data)
      part <- result$participant_features

      expect_true("p_mean_velocity" %in% names(part))
      expect_true("p_velocity_sd" %in% names(part))
      expect_true("p_velocity_cv" %in% names(part))
      expect_true("p_rir_slope" %in% names(part))
      expect_true("p_load_sensitivity" %in% names(part))
    })

    it("creates cross-level features", {
      fe <- FeatureEngineer$new()
      data <- create_test_data()

      result <- fe$engineer_features(data)

      # Participant × Set
      expect_true("ps_deviation_from_p_mean" %in% names(result$participant_set_features))

      # Participant × Load
      expect_true("pl_mean_velocity" %in% names(result$participant_load_features))

      # Participant × RIR
      expect_true("pr_mean_velocity" %in% names(result$participant_rir_features))
    })
  })

  describe("velocity_z_global", {

    it("has approximately mean 0 and sd 1", {
      fe <- FeatureEngineer$new()
      data <- create_test_data(n_participants = 10)

      result <- fe$engineer_features(data)
      z_scores <- result$observation_features$velocity_z_global

      expect_equal(mean(z_scores, na.rm = TRUE), 0, tolerance = 0.1)
      expect_equal(sd(z_scores, na.rm = TRUE), 1, tolerance = 0.1)
    })

    it("correctly identifies high and low velocities", {
      fe <- FeatureEngineer$new()
      data <- create_test_data()

      result <- fe$engineer_features(data)
      obs <- result$observation_features

      # High velocities should have positive z-scores
      high_velocity_obs <- obs[obs$mean_velocity > median(obs$mean_velocity, na.rm = TRUE), ]
      low_velocity_obs <- obs[obs$mean_velocity < median(obs$mean_velocity, na.rm = TRUE), ]

      expect_gt(mean(high_velocity_obs$velocity_z_global, na.rm = TRUE), 0)
      expect_lt(mean(low_velocity_obs$velocity_z_global, na.rm = TRUE), 0)
    })
  })

  describe("set_position_normalized", {

    it("ranges from 0 to 1 within sets", {
      fe <- FeatureEngineer$new()
      data <- create_test_data()

      result <- fe$engineer_features(data)
      pos <- result$observation_features$set_position_normalized

      expect_true(all(pos >= 0 & pos <= 1, na.rm = TRUE))
    })

    it("equals 1 for last rep in set", {
      fe <- FeatureEngineer$new()
      data <- create_test_data()

      result <- fe$engineer_features(data)
      obs <- result$observation_features

      # Last rep should have position = 1 (is_last_rep is logical TRUE/FALSE)
      last_reps <- obs[obs$is_last_rep == TRUE, ]

      expect_true(all(last_reps$set_position_normalized == 1, na.rm = TRUE))
    })
  })

  describe("p_velocity_cv", {

    it("is calculated correctly", {
      fe <- FeatureEngineer$new()
      data <- create_test_data()

      result <- fe$engineer_features(data)
      part <- result$participant_features

      # CV = sd / mean, should be positive
      expect_true(all(part$p_velocity_cv > 0, na.rm = TRUE))
      expect_true(all(part$p_velocity_cv < 1, na.rm = TRUE))  # Reasonable range
    })
  })

  describe("p_load_sensitivity", {

    it("is calculated for participants with multiple loads", {
      fe <- FeatureEngineer$new()
      data <- create_test_data()

      result <- fe$engineer_features(data)
      part <- result$participant_features

      # Most participants should have load sensitivity calculated
      non_na_count <- sum(!is.na(part$p_load_sensitivity))
      expect_gt(non_na_count, 0)
    })

    it("is negative when velocity decreases with load", {
      fe <- FeatureEngineer$new()

      # Create data with clear load effect
      data <- data.frame(
        id = rep("P01", 10),
        sex = "male",
        day = "Day 1",
        load_percentage = rep(c("80%", "90%"), each = 5),
        weight_kg = rep(c(120, 140), each = 5),
        rir = rep(c(4, 3, 2, 1, 0), 2),
        mean_velocity = c(
          0.5, 0.48, 0.45, 0.42, 0.38,  # 80% load
          0.4, 0.38, 0.35, 0.32, 0.28   # 90% load (slower)
        ),
        set_id = rep(c("P01_Day1_80pct_S1", "P01_Day1_90pct_S1"), each = 5),
        rep_number = rep(1:5, 2),
        reps_to_failure = 5
      )

      result <- fe$engineer_features(data)

      # Sensitivity should be negative (higher load = lower velocity)
      expect_lt(result$participant_features$p_load_sensitivity, 0)
    })
  })

  describe("p_rir_slope", {

    it("is positive for typical velocity-RIR relationship", {
      fe <- FeatureEngineer$new()
      data <- create_test_data()

      result <- fe$engineer_features(data)
      part <- result$participant_features

      # Velocity increases with RIR, so slope should be positive
      slopes <- part$p_rir_slope[!is.na(part$p_rir_slope)]
      expect_gt(mean(slopes), 0)
    })
  })

  describe("velocity_z_within_participant", {

    it("captures within-participant variability", {
      fe <- FeatureEngineer$new()
      data <- create_test_data()

      result <- fe$engineer_features(data)
      obs <- result$observation_features

      # Within each participant, z-scores should have mean ~0
      for (pid in unique(obs$id)) {
        p_data <- obs[obs$id == pid, ]
        p_mean <- mean(p_data$velocity_z_within_participant, na.rm = TRUE)
        expect_equal(p_mean, 0, tolerance = 0.5)
      }
    })
  })

  describe("get_observation_feature_names", {

    it("returns expected feature names", {
      fe <- FeatureEngineer$new()

      names <- fe$get_observation_feature_names()

      expect_true("velocity_z_global" %in% names)
      expect_true("velocity_z_within_participant" %in% names)
      expect_true("velocity_z_within_set" %in% names)
      expect_true("set_position_normalized" %in% names)
      expect_gt(length(names), 5)
    })
  })

  describe("get_participant_feature_names", {

    it("returns expected feature names", {
      fe <- FeatureEngineer$new()

      names <- fe$get_participant_feature_names()

      expect_true("p_mean_velocity" %in% names)
      expect_true("p_velocity_cv" %in% names)
      expect_gt(length(names), 5)
    })
  })

  describe("edge cases", {

    it("handles single participant", {
      fe <- FeatureEngineer$new()
      data <- create_test_data(n_participants = 1)

      result <- fe$engineer_features(data)

      expect_equal(nrow(result$participant_features), 1)
      expect_gt(nrow(result$observation_features), 0)
    })

    it("handles missing velocity values gracefully", {
      fe <- FeatureEngineer$new()
      data <- create_test_data()

      # Introduce some NAs
      data$mean_velocity[c(1, 5, 10)] <- NA

      result <- fe$engineer_features(data)

      # Should complete without error
      expect_s3_class(result, "FeatureEngineeringResult")
    })

    it("handles constant velocity within participant", {
      fe <- FeatureEngineer$new()

      data <- data.frame(
        id = rep("P01", 5),
        sex = "male",
        day = "Day 1",
        load_percentage = "80%",
        weight_kg = 120,
        rir = c(4, 3, 2, 1, 0),
        mean_velocity = rep(0.5, 5),  # Constant velocity
        set_id = "P01_Day1_80pct_S1",
        rep_number = 1:5,
        reps_to_failure = 5
      )

      result <- fe$engineer_features(data)

      # velocity_z_within_participant should be 0 (or handled gracefully)
      expect_true(all(
        result$observation_features$velocity_z_within_participant == 0 |
        is.na(result$observation_features$velocity_z_within_participant)
      ))
    })
  })

  describe("reproducibility", {

    it("produces same results for same input", {
      fe <- FeatureEngineer$new()
      data <- create_test_data()

      result1 <- fe$engineer_features(data)
      result2 <- fe$engineer_features(data)

      expect_equal(
        result1$observation_features$velocity_z_global,
        result2$observation_features$velocity_z_global
      )

      expect_equal(
        result1$participant_features$p_velocity_cv,
        result2$participant_features$p_velocity_cv
      )
    })
  })

  describe("global_stats", {

    it("contains population-level statistics", {
      fe <- FeatureEngineer$new()
      data <- create_test_data()

      result <- fe$engineer_features(data)
      stats <- result$global_stats

      expect_true("mean_velocity" %in% names(stats))
      expect_true("sd_velocity" %in% names(stats))
    })
  })
})
