# tests/testthat/test-treatment_group.R
# Unit tests for TreatmentGroup domain entity

box::use(
  testthat[...],
  ../../R/domain/treatment_group[TreatmentGroup]
)

describe("TreatmentGroup", {

  describe("constructor", {
    it("creates a valid treatment group with all parameters", {
      group <- TreatmentGroup$new(
        group_id = "study1_group1",
        repetitions_in_reserve = 3,
        mean_pre = 100,
        mean_post = 110,
        standard_deviation_pre = 15,
        sample_size = 30,
        pre_post_correlation = 0.5
      )

      expect_equal(group$group_id, "study1_group1")
      expect_equal(group$repetitions_in_reserve, 3)
      expect_equal(group$mean_pre, 100)
      expect_equal(group$mean_post, 110)
      expect_equal(group$standard_deviation_pre, 15)
      expect_equal(group$sample_size, 30)
      expect_equal(group$pre_post_correlation, 0.5)
    })

    it("allows null correlation", {
      group <- TreatmentGroup$new(
        group_id = "test",
        repetitions_in_reserve = 0,
        mean_pre = 100,
        mean_post = 110,
        standard_deviation_pre = 15,
        sample_size = 20,
        pre_post_correlation = NULL
      )

      expect_null(group$pre_post_correlation)
      expect_false(group$has_correlation)
    })
  })

  describe("calculate_mean_change", {
    it("returns correct mean change", {
      group <- TreatmentGroup$new(
        group_id = "test",
        repetitions_in_reserve = 2,
        mean_pre = 100,
        mean_post = 115,
        standard_deviation_pre = 10,
        sample_size = 25
      )

      expect_equal(group$calculate_mean_change(), 15)
    })

    it("handles negative change (decrease)", {
      group <- TreatmentGroup$new(
        group_id = "test",
        repetitions_in_reserve = 5,
        mean_pre = 100,
        mean_post = 95,
        standard_deviation_pre = 10,
        sample_size = 25
      )

      expect_equal(group$calculate_mean_change(), -5)
    })
  })

  describe("validation", {
    it("passes validation for valid data", {
      group <- TreatmentGroup$new(
        group_id = "test",
        repetitions_in_reserve = 3,
        mean_pre = 100,
        mean_post = 110,
        standard_deviation_pre = 15,
        sample_size = 30
      )

      result <- group$validate()

      expect_true(result$is_valid)
      expect_length(result$errors, 0)
    })

    it("fails validation for negative RIR", {
      group <- TreatmentGroup$new(
        group_id = "test",
        repetitions_in_reserve = -1,
        mean_pre = 100,
        mean_post = 110,
        standard_deviation_pre = 15,
        sample_size = 30
      )

      result <- group$validate()

      expect_false(result$is_valid)
      expect_true(any(grepl("Repetitions in reserve", result$errors)))
    })

    it("fails validation for non-positive standard deviation", {
      group <- TreatmentGroup$new(
        group_id = "test",
        repetitions_in_reserve = 3,
        mean_pre = 100,
        mean_post = 110,
        standard_deviation_pre = 0,
        sample_size = 30
      )

      result <- group$validate()

      expect_false(result$is_valid)
      expect_true(any(grepl("Standard deviation", result$errors)))
    })

    it("fails validation for sample size less than 2", {
      group <- TreatmentGroup$new(
        group_id = "test",
        repetitions_in_reserve = 3,
        mean_pre = 100,
        mean_post = 110,
        standard_deviation_pre = 15,
        sample_size = 1
      )

      result <- group$validate()

      expect_false(result$is_valid)
      expect_true(any(grepl("Sample size", result$errors)))
    })
  })

  describe("active bindings", {
    it("correctly identifies training to failure", {
      failure_group <- TreatmentGroup$new(
        group_id = "test",
        repetitions_in_reserve = 0,
        mean_pre = 100,
        mean_post = 110,
        standard_deviation_pre = 15,
        sample_size = 30
      )

      rir_group <- TreatmentGroup$new(
        group_id = "test",
        repetitions_in_reserve = 3,
        mean_pre = 100,
        mean_post = 110,
        standard_deviation_pre = 15,
        sample_size = 30
      )

      expect_true(failure_group$is_training_to_failure)
      expect_false(rir_group$is_training_to_failure)
    })
  })

  describe("to_list", {
    it("returns complete list representation", {
      group <- TreatmentGroup$new(
        group_id = "test",
        repetitions_in_reserve = 2,
        mean_pre = 100,
        mean_post = 115,
        standard_deviation_pre = 10,
        sample_size = 25,
        pre_post_correlation = 0.6
      )

      result <- group$to_list()

      expect_type(result, "list")
      expect_equal(result$group_id, "test")
      expect_equal(result$mean_change, 15)
    })
  })
})
