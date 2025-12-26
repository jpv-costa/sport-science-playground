# tests/testthat/test-effect_size_calculator.R
# Unit tests for effect size calculators

box::use(
  testthat[...],
  ../../R/domain/treatment_group[TreatmentGroup],
  ../../R/domain/effect_size[EffectSize],
  ../../R/calculators/effect_size_calculator[SMCRCalculator, ROMCCalculator]
)

describe("SMCRCalculator", {

  describe("calculate", {
    it("computes correct effect size for known values", {
      calculator <- SMCRCalculator$new(
        default_correlation = 0.5,
        apply_bias_correction = TRUE
      )

      group <- TreatmentGroup$new(
        group_id = "test",
        repetitions_in_reserve = 3,
        mean_pre = 100,
        mean_post = 110,
        standard_deviation_pre = 15,
        sample_size = 30,
        pre_post_correlation = 0.5
      )

      result <- calculator$calculate(group)

      # Expected: (110-100)/15 = 0.667
      # With Hedge's J correction for n=30: J ≈ 0.974
      # Corrected ES ≈ 0.667 * 0.974 ≈ 0.650
      expect_s3_class(result, "EffectSize")
      expect_equal(result$effect_estimate, 0.650, tolerance = 0.01)
      expect_equal(result$effect_type, "smcr")
      expect_true(result$is_bias_corrected)
    })

    it("produces positive effect for improvement", {
      calculator <- SMCRCalculator$new()
      group <- TreatmentGroup$new(
        group_id = "test",
        repetitions_in_reserve = 0,
        mean_pre = 100,
        mean_post = 120,
        standard_deviation_pre = 20,
        sample_size = 25
      )

      result <- calculator$calculate(group)

      expect_gt(result$effect_estimate, 0)
    })

    it("produces negative effect for decline", {
      calculator <- SMCRCalculator$new()
      group <- TreatmentGroup$new(
        group_id = "test",
        repetitions_in_reserve = 5,
        mean_pre = 100,
        mean_post = 90,
        standard_deviation_pre = 15,
        sample_size = 30
      )

      result <- calculator$calculate(group)

      expect_lt(result$effect_estimate, 0)
    })

    it("uses default correlation when not provided", {
      calculator <- SMCRCalculator$new(default_correlation = 0.7)
      group <- TreatmentGroup$new(
        group_id = "test",
        repetitions_in_reserve = 2,
        mean_pre = 100,
        mean_post = 110,
        standard_deviation_pre = 15,
        sample_size = 30,
        pre_post_correlation = NULL
      )

      result <- calculator$calculate(group)

      expect_equal(result$pre_post_correlation, 0.7)
    })

    it("applies smaller correction for larger samples", {
      calculator <- SMCRCalculator$new()

      small_sample <- TreatmentGroup$new(
        group_id = "small",
        repetitions_in_reserve = 3,
        mean_pre = 100,
        mean_post = 110,
        standard_deviation_pre = 15,
        sample_size = 10
      )

      large_sample <- TreatmentGroup$new(
        group_id = "large",
        repetitions_in_reserve = 3,
        mean_pre = 100,
        mean_post = 110,
        standard_deviation_pre = 15,
        sample_size = 100
      )

      small_result <- calculator$calculate(small_sample)
      large_result <- calculator$calculate(large_sample)

      # Larger sample should have effect closer to raw (less correction)
      expect_lt(abs(small_result$effect_estimate), abs(large_result$effect_estimate))
    })

    it("rejects invalid treatment group", {
      calculator <- SMCRCalculator$new()
      group <- TreatmentGroup$new(
        group_id = "invalid",
        repetitions_in_reserve = 3,
        mean_pre = 100,
        mean_post = 110,
        standard_deviation_pre = 0,  # Invalid: zero SD
        sample_size = 30
      )

      expect_error(
        calculator$calculate(group),
        "Invalid treatment group"
      )
    })

    it("rejects non-TreatmentGroup input", {
      calculator <- SMCRCalculator$new()

      expect_error(
        calculator$calculate(list(mean_pre = 100, mean_post = 110)),
        "Input must be a TreatmentGroup object"
      )
    })
  })

  describe("calculate_batch", {
    it("processes multiple groups", {
      calculator <- SMCRCalculator$new()

      groups <- list(
        TreatmentGroup$new("g1", 0, 100, 115, 15, 30),
        TreatmentGroup$new("g2", 3, 100, 110, 15, 25),
        TreatmentGroup$new("g3", 5, 100, 105, 15, 20)
      )

      results <- calculator$calculate_batch(groups)

      expect_length(results, 3)
      expect_true(all(sapply(results, inherits, "EffectSize")))
    })
  })
})

describe("ROMCCalculator", {

  describe("calculate", {
    it("computes correct log ratio", {
      calculator <- ROMCCalculator$new()
      group <- TreatmentGroup$new(
        group_id = "test",
        repetitions_in_reserve = 2,
        mean_pre = 100,
        mean_post = 110,
        standard_deviation_pre = 15,
        sample_size = 30
      )

      result <- calculator$calculate(group)

      # Expected: log(110/100) = log(1.1) ≈ 0.0953
      expect_equal(result$effect_estimate, log(1.1), tolerance = 0.001)
      expect_equal(result$effect_type, "romc")
    })

    it("rejects non-positive pre-mean", {
      calculator <- ROMCCalculator$new()
      group <- TreatmentGroup$new(
        group_id = "test",
        repetitions_in_reserve = 2,
        mean_pre = 0,  # Invalid for log ratio
        mean_post = 110,
        standard_deviation_pre = 15,
        sample_size = 30
      )

      expect_error(
        calculator$calculate(group),
        "Pre-intervention mean must be positive"
      )
    })
  })
})
