# tests/testthat/test-velocity_data_loader.R
# Unit tests for VelocityDataLoader

box::use(
  testthat[...],
  ../../R/loaders/velocity_data_loader[VelocityDataLoader]
)

describe("VelocityDataLoader", {

  describe("initialization", {
    it("accepts data path", {
      loader <- VelocityDataLoader$new("path/to/data.xlsx")

      expect_s3_class(loader, "VelocityDataLoader")
    })
  })

  describe("load", {
    it("loads and cleans data from Excel file", {
      skip_if_not_installed("readxl")
      skip_if_not(
        file.exists("data/external/peerj_2025_velocity/peerj_velocity_rir_data.xlsx"),
        "Test data file not available"
      )

      loader <- VelocityDataLoader$new(
        "data/external/peerj_2025_velocity/peerj_velocity_rir_data.xlsx"
      )

      data <- loader$load()

      expect_s3_class(data, "data.frame")
      expect_true("id" %in% names(data))
      expect_true("exercise" %in% names(data))
      expect_true("mean_velocity" %in% names(data))
      expect_true("perceived_rir" %in% names(data))
    })

    it("cleans column names to snake_case", {
      skip_if_not_installed("readxl")
      skip_if_not(
        file.exists("data/external/peerj_2025_velocity/peerj_velocity_rir_data.xlsx"),
        "Test data file not available"
      )

      loader <- VelocityDataLoader$new(
        "data/external/peerj_2025_velocity/peerj_velocity_rir_data.xlsx"
      )

      data <- loader$load()

      # Column names should be lowercase with underscores
      expect_true(all(names(data) == tolower(names(data))))
      expect_false(any(grepl(" ", names(data))))
    })
  })

  describe("summarize", {
    it("computes summary statistics", {
      skip_if_not_installed("readxl")
      skip_if_not(
        file.exists("data/external/peerj_2025_velocity/peerj_velocity_rir_data.xlsx"),
        "Test data file not available"
      )

      loader <- VelocityDataLoader$new(
        "data/external/peerj_2025_velocity/peerj_velocity_rir_data.xlsx"
      )

      data <- loader$load()
      summary <- loader$summarize(data)

      expect_true("n_observations" %in% names(summary))
      expect_true("n_participants" %in% names(summary))
      expect_true("exercises" %in% names(summary))
      expect_true("velocity_range" %in% names(summary))
      expect_true("rir_range" %in% names(summary))

      expect_gt(summary$n_observations, 0)
      expect_gt(summary$n_participants, 0)
      expect_length(summary$velocity_range, 2)
    })

    it("identifies correct exercises from PeerJ data", {
      skip_if_not_installed("readxl")
      skip_if_not(
        file.exists("data/external/peerj_2025_velocity/peerj_velocity_rir_data.xlsx"),
        "Test data file not available"
      )

      loader <- VelocityDataLoader$new(
        "data/external/peerj_2025_velocity/peerj_velocity_rir_data.xlsx"
      )

      data <- loader$load()
      summary <- loader$summarize(data)

      # Paper includes squat and bench press
      expect_true(length(summary$exercises) >= 2)
    })
  })
})
