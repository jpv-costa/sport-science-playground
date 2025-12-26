# tests/testthat/test-rir_velocity_data_loader.R
# Unit tests for RirVelocityDataLoader

box::use(
  testthat[...],
  ../../R/loaders/rir_velocity_data_loader[RirVelocityDataLoader]
)

describe("RirVelocityDataLoader", {

  describe("initialization", {
    it("accepts data path", {
      loader <- RirVelocityDataLoader$new("path/to/data.xlsx")

      expect_s3_class(loader, "RirVelocityDataLoader")
    })
  })

  describe("load", {
    it("loads and preprocesses RIR-velocity data", {
      skip_if_not_installed("readxl")
      skip_if_not(
        file.exists("data/external/zourdos_rir/rir_velocity_dataset.xlsx"),
        "Test data file not available"
      )

      loader <- RirVelocityDataLoader$new(
        "data/external/zourdos_rir/rir_velocity_dataset.xlsx"
      )

      data <- loader$load()

      expect_s3_class(data, "data.frame")
      expect_true("id" %in% names(data))
      expect_true("sex" %in% names(data))
      expect_true("mean_velocity" %in% names(data))
      expect_true("rir" %in% names(data))
      expect_true("set_type" %in% names(data))
      expect_true("day" %in% names(data))
    })

    it("derives velocity metrics", {
      skip_if_not_installed("readxl")
      skip_if_not(
        file.exists("data/external/zourdos_rir/rir_velocity_dataset.xlsx"),
        "Test data file not available"
      )

      loader <- RirVelocityDataLoader$new(
        "data/external/zourdos_rir/rir_velocity_dataset.xlsx"
      )

      data <- loader$load()

      expect_true("max_velocity" %in% names(data))
      expect_true("velocity_percentage" %in% names(data))
      expect_true("max_reps" %in% names(data))
    })
  })

  describe("filter_by_day", {
    it("filters to specified day", {
      skip_if_not_installed("readxl")
      skip_if_not(
        file.exists("data/external/zourdos_rir/rir_velocity_dataset.xlsx"),
        "Test data file not available"
      )

      loader <- RirVelocityDataLoader$new(
        "data/external/zourdos_rir/rir_velocity_dataset.xlsx"
      )

      data <- loader$load()
      day1_data <- loader$filter_by_day(data, "Day 1")

      expect_true(all(day1_data$day == "Day 1"))
      expect_gt(nrow(day1_data), 0)
    })
  })

  describe("filter_by_load", {
    it("filters to specified load", {
      skip_if_not_installed("readxl")
      skip_if_not(
        file.exists("data/external/zourdos_rir/rir_velocity_dataset.xlsx"),
        "Test data file not available"
      )

      loader <- RirVelocityDataLoader$new(
        "data/external/zourdos_rir/rir_velocity_dataset.xlsx"
      )

      data <- loader$load()
      rtf90_data <- loader$filter_by_load(data, "RTF90")

      expect_true(all(rtf90_data$set_type == "RTF90"))
      expect_gt(nrow(rtf90_data), 0)
    })
  })

  describe("summarize", {
    it("computes summary statistics", {
      skip_if_not_installed("readxl")
      skip_if_not(
        file.exists("data/external/zourdos_rir/rir_velocity_dataset.xlsx"),
        "Test data file not available"
      )

      loader <- RirVelocityDataLoader$new(
        "data/external/zourdos_rir/rir_velocity_dataset.xlsx"
      )

      data <- loader$load()
      summary <- loader$summarize(data)

      expect_true("n_observations" %in% names(summary))
      expect_true("n_participants" %in% names(summary))
      expect_true("load_types" %in% names(summary))
      expect_true("velocity_range" %in% names(summary))
      expect_true("rir_range" %in% names(summary))

      expect_gt(summary$n_observations, 0)
      expect_gt(summary$n_participants, 0)
      expect_length(summary$velocity_range, 2)
    })

    it("identifies correct load types", {
      skip_if_not_installed("readxl")
      skip_if_not(
        file.exists("data/external/zourdos_rir/rir_velocity_dataset.xlsx"),
        "Test data file not available"
      )

      loader <- RirVelocityDataLoader$new(
        "data/external/zourdos_rir/rir_velocity_dataset.xlsx"
      )

      data <- loader$load()
      summary <- loader$summarize(data)

      # Should have RTF70, RTF80, RTF90
      expect_true(length(summary$load_types) >= 3)
    })
  })
})
