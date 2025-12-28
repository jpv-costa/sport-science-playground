# tests/testthat/test-deadlift_rir_data_loader.R
# Unit tests for DeadliftRirDataLoader

box::use(
  testthat[...],
  ../../R/loaders/deadlift_rir_data_loader[DeadliftRirDataLoader]
)

describe("DeadliftRirDataLoader", {

  describe("initialization", {
    it("accepts data path", {
      loader <- DeadliftRirDataLoader$new("path/to/data.xlsx")

      expect_s3_class(loader, "DeadliftRirDataLoader")
    })
  })

  describe("load", {
    it("loads and preprocesses deadlift data from multi-sheet Excel", {
      skip_if_not_installed("readxl")
      test_file <- "../../deadlift-study/TESE - Análise estatística (Preliminar).xlsx"
      skip_if_not(
        file.exists(test_file),
        "Test data file not available"
      )

      loader <- DeadliftRirDataLoader$new(test_file)

      data <- loader$load()

      expect_s3_class(data, "data.frame")
      expect_true("id" %in% names(data))
      expect_true("sex" %in% names(data))
      expect_true("mean_velocity" %in% names(data))
      expect_true("rir" %in% names(data))
      expect_true("load_percentage" %in% names(data))
      expect_true("day" %in% names(data))
      expect_true("weight_kg" %in% names(data))
    })

    it("returns correct number of participants", {
      skip_if_not_installed("readxl")
      test_file <- "../../deadlift-study/TESE - Análise estatística (Preliminar).xlsx"
      skip_if_not(
        file.exists(test_file),
        "Test data file not available"
      )

      loader <- DeadliftRirDataLoader$new(test_file)

      data <- loader$load()

      expect_equal(length(unique(data$id)), 19)
    })

    it("correctly identifies sex from sheet names", {
      skip_if_not_installed("readxl")
      test_file <- "../../deadlift-study/TESE - Análise estatística (Preliminar).xlsx"
      skip_if_not(
        file.exists(test_file),
        "Test data file not available"
      )

      loader <- DeadliftRirDataLoader$new(test_file)

      data <- loader$load()

      male_participants <- unique(data$id[data$sex == "male"])
      female_participants <- unique(data$id[data$sex == "female"])

      expect_equal(length(male_participants), 15)
      expect_equal(length(female_participants), 4)
    })

    it("extracts both load percentages", {
      skip_if_not_installed("readxl")
      test_file <- "../../deadlift-study/TESE - Análise estatística (Preliminar).xlsx"
      skip_if_not(
        file.exists(test_file),
        "Test data file not available"
      )

      loader <- DeadliftRirDataLoader$new(test_file)

      data <- loader$load()

      expect_true("80%" %in% data$load_percentage)
      expect_true("90%" %in% data$load_percentage)
    })

    it("extracts both days", {
      skip_if_not_installed("readxl")
      test_file <- "../../deadlift-study/TESE - Análise estatística (Preliminar).xlsx"
      skip_if_not(
        file.exists(test_file),
        "Test data file not available"
      )

      loader <- DeadliftRirDataLoader$new(test_file)

      data <- loader$load()

      expect_true("Day 1" %in% data$day)
      expect_true("Day 2" %in% data$day)
    })
  })

  describe("filter_by_day", {
    it("filters to specified day", {
      skip_if_not_installed("readxl")
      test_file <- "../../deadlift-study/TESE - Análise estatística (Preliminar).xlsx"
      skip_if_not(
        file.exists(test_file),
        "Test data file not available"
      )

      loader <- DeadliftRirDataLoader$new(test_file)

      data <- loader$load()
      day1_data <- loader$filter_by_day(data, "Day 1")

      expect_true(all(day1_data$day == "Day 1"))
      expect_gt(nrow(day1_data), 0)
    })
  })

  describe("filter_by_load", {
    it("filters to specified load", {
      skip_if_not_installed("readxl")
      test_file <- "../../deadlift-study/TESE - Análise estatística (Preliminar).xlsx"
      skip_if_not(
        file.exists(test_file),
        "Test data file not available"
      )

      loader <- DeadliftRirDataLoader$new(test_file)

      data <- loader$load()
      load_80_data <- loader$filter_by_load(data, "80%")

      expect_true(all(load_80_data$load_percentage == "80%"))
      expect_gt(nrow(load_80_data), 0)
    })
  })

  describe("summarize", {
    it("computes summary statistics", {
      skip_if_not_installed("readxl")
      test_file <- "../../deadlift-study/TESE - Análise estatística (Preliminar).xlsx"
      skip_if_not(
        file.exists(test_file),
        "Test data file not available"
      )

      loader <- DeadliftRirDataLoader$new(test_file)

      data <- loader$load()
      summary <- loader$summarize(data)

      expect_true("n_observations" %in% names(summary))
      expect_true("n_participants" %in% names(summary))
      expect_true("load_types" %in% names(summary))
      expect_true("velocity_range" %in% names(summary))
      expect_true("rir_range" %in% names(summary))
      expect_true("weight_range" %in% names(summary))

      expect_gt(summary$n_observations, 0)
      expect_equal(summary$n_participants, 19)
      expect_length(summary$velocity_range, 2)
    })

    it("identifies correct load types", {
      skip_if_not_installed("readxl")
      test_file <- "../../deadlift-study/TESE - Análise estatística (Preliminar).xlsx"
      skip_if_not(
        file.exists(test_file),
        "Test data file not available"
      )

      loader <- DeadliftRirDataLoader$new(test_file)

      data <- loader$load()
      summary <- loader$summarize(data)

      expect_true(length(summary$load_types) == 2)
      expect_true("80%" %in% summary$load_types)
      expect_true("90%" %in% summary$load_types)
    })
  })

  describe("anonymization", {
    it("anonymizes participant IDs when anonymize = TRUE", {
      skip_if_not_installed("readxl")
      test_file <- "../../deadlift-study/TESE - Análise estatística (Preliminar).xlsx"
      skip_if_not(
        file.exists(test_file),
        "Test data file not available"
      )

      loader <- DeadliftRirDataLoader$new(test_file)

      data <- loader$load(anonymize = TRUE)

      # All IDs should be in P01-P19 format
      expect_true(all(grepl("^P[0-9]{2}$", data$id)))
    })

    it("anonymizes set_id column when anonymize = TRUE", {
      skip_if_not_installed("readxl")
      test_file <- "../../deadlift-study/TESE - Análise estatística (Preliminar).xlsx"
      skip_if_not(
        file.exists(test_file),
        "Test data file not available"
      )

      loader <- DeadliftRirDataLoader$new(test_file)

      data <- loader$load(anonymize = TRUE)

      # All set_ids should start with P01-P19 format
      expect_true(all(grepl("^P[0-9]{2}_", data$set_id)))
    })

    it("keeps original names when anonymize = FALSE (default)", {
      skip_if_not_installed("readxl")
      test_file <- "../../deadlift-study/TESE - Análise estatística (Preliminar).xlsx"
      skip_if_not(
        file.exists(test_file),
        "Test data file not available"
      )

      loader <- DeadliftRirDataLoader$new(test_file)

      data <- loader$load(anonymize = FALSE)

      # IDs should NOT be in P01-P19 format
      expect_false(all(grepl("^P[0-9]{2}$", data$id)))
    })

    it("maintains same number of participants after anonymization", {
      skip_if_not_installed("readxl")
      test_file <- "../../deadlift-study/TESE - Análise estatística (Preliminar).xlsx"
      skip_if_not(
        file.exists(test_file),
        "Test data file not available"
      )

      loader <- DeadliftRirDataLoader$new(test_file)

      data_original <- loader$load(anonymize = FALSE)
      data_anon <- loader$load(anonymize = TRUE)

      expect_equal(length(unique(data_anon$id)), length(unique(data_original$id)))
      expect_equal(nrow(data_anon), nrow(data_original))
    })
  })

  describe("data quality", {
    it("velocity values are in valid range", {
      skip_if_not_installed("readxl")
      test_file <- "../../deadlift-study/TESE - Análise estatística (Preliminar).xlsx"
      skip_if_not(
        file.exists(test_file),
        "Test data file not available"
      )

      loader <- DeadliftRirDataLoader$new(test_file)

      data <- loader$load()

      expect_true(all(data$mean_velocity > 0))
      expect_true(all(data$mean_velocity < 1.5))  # Reasonable deadlift velocity range
    })

    it("RIR values are in valid range", {
      skip_if_not_installed("readxl")
      test_file <- "../../deadlift-study/TESE - Análise estatística (Preliminar).xlsx"
      skip_if_not(
        file.exists(test_file),
        "Test data file not available"
      )

      loader <- DeadliftRirDataLoader$new(test_file)

      data <- loader$load()

      expect_true(all(data$rir >= 0))
      expect_true(all(data$rir <= 10))  # Reasonable RIR range
    })

    it("weight values are in valid range for deadlift", {
      skip_if_not_installed("readxl")
      test_file <- "../../deadlift-study/TESE - Análise estatística (Preliminar).xlsx"
      skip_if_not(
        file.exists(test_file),
        "Test data file not available"
      )

      loader <- DeadliftRirDataLoader$new(test_file)

      data <- loader$load()

      expect_true(all(data$weight_kg >= 50))   # Minimum reasonable deadlift
      expect_true(all(data$weight_kg <= 300))  # Maximum reasonable deadlift
    })
  })
})
