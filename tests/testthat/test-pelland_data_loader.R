# tests/testthat/test-pelland_data_loader.R
# Unit tests for PellandDataLoader

box::use(
  testthat[...],
  ../../R/loaders/pelland_data_loader[PellandDataLoader]
)

describe("PellandDataLoader", {

  describe("initialization", {
    it("accepts data path", {
      loader <- PellandDataLoader$new("path/to/data.xlsx")

      expect_s3_class(loader, "PellandDataLoader")
    })
  })

  describe("load", {
    it("loads and preprocesses Pelland meta-analysis data", {
      skip_if_not_installed("readxl")
      skip_if_not(
        file.exists("data/external/pelland_2022_rtf/41598_2022_16909_MOESM1_ESM.xlsx"),
        "Test data file not available"
      )

      loader <- PellandDataLoader$new(
        "data/external/pelland_2022_rtf/41598_2022_16909_MOESM1_ESM.xlsx"
      )

      data <- loader$load()

      expect_s3_class(data, "data.frame")
      expect_true("study" %in% names(data))
      expect_true("group" %in% names(data))
      expect_true("outcome" %in% names(data))
      expect_true("avg.rir" %in% names(data))
    })

    it("creates observation index column", {
      skip_if_not_installed("readxl")
      skip_if_not(
        file.exists("data/external/pelland_2022_rtf/41598_2022_16909_MOESM1_ESM.xlsx"),
        "Test data file not available"
      )

      loader <- PellandDataLoader$new(
        "data/external/pelland_2022_rtf/41598_2022_16909_MOESM1_ESM.xlsx"
      )

      data <- loader$load()

      expect_true("obs" %in% names(data))
      expect_equal(data$obs, seq_len(nrow(data)))
    })
  })

  describe("filter_by_outcome", {
    it("filters to strength outcomes", {
      skip_if_not_installed("readxl")
      skip_if_not(
        file.exists("data/external/pelland_2022_rtf/41598_2022_16909_MOESM1_ESM.xlsx"),
        "Test data file not available"
      )

      loader <- PellandDataLoader$new(
        "data/external/pelland_2022_rtf/41598_2022_16909_MOESM1_ESM.xlsx"
      )

      data <- loader$load()
      strength_data <- loader$filter_by_outcome(data, "Strength")

      expect_true(all(strength_data$outcome == "Strength"))
      expect_gt(nrow(strength_data), 0)
    })

    it("filters to hypertrophy outcomes", {
      skip_if_not_installed("readxl")
      skip_if_not(
        file.exists("data/external/pelland_2022_rtf/41598_2022_16909_MOESM1_ESM.xlsx"),
        "Test data file not available"
      )

      loader <- PellandDataLoader$new(
        "data/external/pelland_2022_rtf/41598_2022_16909_MOESM1_ESM.xlsx"
      )

      data <- loader$load()
      hypertrophy_data <- loader$filter_by_outcome(data, "Hypertrophy")

      expect_true(all(hypertrophy_data$outcome == "Hypertrophy"))
      expect_gt(nrow(hypertrophy_data), 0)
    })
  })

  describe("summarize", {
    it("computes summary statistics", {
      skip_if_not_installed("readxl")
      skip_if_not(
        file.exists("data/external/pelland_2022_rtf/41598_2022_16909_MOESM1_ESM.xlsx"),
        "Test data file not available"
      )

      loader <- PellandDataLoader$new(
        "data/external/pelland_2022_rtf/41598_2022_16909_MOESM1_ESM.xlsx"
      )

      data <- loader$load()
      summary <- loader$summarize(data)

      expect_true("n_total_effects" %in% names(summary))
      expect_true("n_studies" %in% names(summary))
      expect_true("n_strength_effects" %in% names(summary))
      expect_true("n_hypertrophy_effects" %in% names(summary))
      expect_true("rir_range" %in% names(summary))

      expect_gt(summary$n_total_effects, 0)
      expect_gt(summary$n_studies, 0)
    })

    it("correctly counts by outcome", {
      skip_if_not_installed("readxl")
      skip_if_not(
        file.exists("data/external/pelland_2022_rtf/41598_2022_16909_MOESM1_ESM.xlsx"),
        "Test data file not available"
      )

      loader <- PellandDataLoader$new(
        "data/external/pelland_2022_rtf/41598_2022_16909_MOESM1_ESM.xlsx"
      )

      data <- loader$load()
      summary <- loader$summarize(data)

      strength_count <- sum(data$outcome == "Strength")
      hypertrophy_count <- sum(data$outcome == "Hypertrophy")

      expect_equal(summary$n_strength_effects, strength_count)
      expect_equal(summary$n_hypertrophy_effects, hypertrophy_count)
    })
  })
})
