# tests/testthat/test-model_comparator.R
# Tests for ModelComparator R6 class
# Following CLAUDE.md testing principles

box::use(
 testthat[describe, it, expect_true, expect_s3_class, expect_type,
           expect_named, expect_error, expect_gte, expect_lte, expect_equal],
  lme4[lmer]
)

# Set box path relative to project root
options(box.path = c("../../R", getOption("box.path")))

# Load the module under test
box::use(
  calculators/model_comparator[ModelComparator, ComparisonResult, compare_models]
)

# Create realistic test data
create_test_data <- function(n_participants = 15, obs_per_participant = 8) {
  set.seed(42)
  participants <- paste0("P", seq_len(n_participants))

  do.call(rbind, lapply(participants, function(pid) {
    intercept <- 0.20 + rnorm(1, 0, 0.03)
    slope <- 0.030 + rnorm(1, 0, 0.005)

    rir <- sample(0:7, obs_per_participant, replace = TRUE)
    velocity <- intercept + slope * rir + rnorm(obs_per_participant, 0, 0.02)

    data.frame(
      id = pid,
      rir = rir,
      mean_velocity = velocity,
      load = sample(c(60, 70, 80), obs_per_participant, replace = TRUE),
      stringsAsFactors = FALSE
    )
  }))
}

describe("ComparisonResult", {
  it("creates result with required fields", {
    result <- ComparisonResult$new(
      model1_name = "Simple",
      model2_name = "Complex",
      aic_diff = 10,
      bic_diff = 12
    )
    expect_s3_class(result, "ComparisonResult")
    expect_equal(result$model1_name, "Simple")
    expect_equal(result$aic_diff, 10)
  })

  it("determines preferred model from Bayes factor", {
    # BF > 1 favors model1
    result <- ComparisonResult$new(
      model1_name = "Null",
      model2_name = "Alternative",
      aic_diff = 5,
      bic_diff = 8,
      bayes_factor = 15
    )
    expect_equal(result$preferred_model, "Null")
    expect_equal(result$evidence_strength, "strong")
  })

  it("interprets weak evidence correctly", {
    result <- ComparisonResult$new(
      model1_name = "A",
      model2_name = "B",
      aic_diff = 1,
      bic_diff = 1,
      bayes_factor = 2
    )
    expect_equal(result$evidence_strength, "weak/anecdotal")
  })

  it("generates interpretation text", {
    result <- ComparisonResult$new(
      model1_name = "Base",
      model2_name = "Full",
      aic_diff = 10,
      bic_diff = 12,
      bayes_factor = 20,
      lrt_chisq = 15.5,
      lrt_df = 2,
      lrt_pvalue = 0.0004
    )
    interp <- result$interpret()
    expect_type(interp, "character")
    expect_true(grepl("Bayes factor", interp))
    expect_true(grepl("ΔAIC", interp))
  })

  it("converts to list", {
    result <- ComparisonResult$new(
      model1_name = "A",
      model2_name = "B",
      aic_diff = 5,
      bic_diff = 7
    )
    lst <- result$to_list()
    expect_type(lst, "list")
    expect_named(lst, c("model1", "model2", "aic_diff", "bic_diff",
                        "bayes_factor", "lrt_chisq", "lrt_df", "lrt_pvalue",
                        "preferred_model", "evidence_strength"))
  })
})

describe("ModelComparator", {
  data <- create_test_data()

  # Fit test models
  model_simple <- lmer(mean_velocity ~ rir + (1 | id), data = data)
  model_complex <- lmer(mean_velocity ~ rir + (1 + rir | id), data = data)

  describe("initialization", {
    it("creates instance", {
      comparator <- ModelComparator$new()
      expect_s3_class(comparator, "ModelComparator")
    })
  })

  describe("compare", {
    comparator <- ModelComparator$new()

    it("compares two models and returns ComparisonResult", {
      result <- comparator$compare(model_simple, model_complex)
      expect_s3_class(result, "ComparisonResult")
    })

    it("calculates AIC and BIC differences", {
      result <- comparator$compare(model_simple, model_complex)
      expect_type(result$aic_diff, "double")
      expect_type(result$bic_diff, "double")
    })

    it("calculates Bayes factor", {
      result <- comparator$compare(model_simple, model_complex)
      expect_type(result$bayes_factor, "double")
      expect_gte(result$bayes_factor, 0)
    })

    it("performs LRT for nested models", {
      result <- comparator$compare(model_simple, model_complex, nested = TRUE)
      # LRT should be available for nested lmer models
      expect_type(result$lrt_chisq, "double")
      expect_type(result$lrt_pvalue, "double")
    })

    it("accepts custom model names", {
      result <- comparator$compare(
        model_simple, model_complex,
        model1_name = "Random Intercept",
        model2_name = "Random Slope"
      )
      expect_equal(result$model1_name, "Random Intercept")
      expect_equal(result$model2_name, "Random Slope")
    })
  })

  describe("compare_multiple", {
    comparator <- ModelComparator$new()
    models <- list(
      simple = model_simple,
      complex = model_complex
    )

    it("compares multiple models", {
      result <- comparator$compare_multiple(models)
      expect_s3_class(result, "data.frame")
      expect_equal(nrow(result), 2)
    })

    it("includes required columns", {
      result <- comparator$compare_multiple(models)
      expect_true(all(c("model", "aic", "bic", "delta_aic", "bayes_factor") %in% names(result)))
    })

    it("uses first model as reference by default", {
      result <- comparator$compare_multiple(models)
      expect_equal(result$delta_aic[1], 0)  # Reference model
    })

    it("accepts custom reference model", {
      result <- comparator$compare_multiple(models, reference_model = "complex")
      expect_equal(result$evidence[result$model == "complex"], "reference")
    })
  })

  describe("calculate_prediction_difference", {
    comparator <- ModelComparator$new()

    it("calculates prediction difference for a predictor", {
      result <- comparator$calculate_prediction_difference(
        model_complex, data, "rir"
      )
      expect_type(result, "list")
      expect_named(result, c("predictor", "low_value", "high_value",
                             "pred_at_low", "pred_at_high",
                             "prediction_difference", "interpretation"))
    })

    it("returns positive difference for positive slope", {
      result <- comparator$calculate_prediction_difference(
        model_complex, data, "rir"
      )
      # RIR has positive slope with velocity
      expect_gte(result$prediction_difference, 0)
    })

    it("provides interpretation text", {
      result <- comparator$calculate_prediction_difference(
        model_complex, data, "rir"
      )
      expect_type(result$interpretation, "character")
      expect_true(grepl("SD increase", result$interpretation))
    })

    it("errors for missing predictor", {
      expect_error(
        comparator$calculate_prediction_difference(model_complex, data, "nonexistent")
      )
    })
  })

  describe("interpret_bayes_factor", {
    comparator <- ModelComparator$new()

    it("interprets decisive evidence for", {
      expect_equal(comparator$interpret_bayes_factor(150), "Decisive evidence for")
    })

    it("interprets strong evidence for", {
      expect_equal(comparator$interpret_bayes_factor(20), "Strong evidence for")
    })

    it("interprets moderate evidence for", {
      expect_equal(comparator$interpret_bayes_factor(5), "Moderate evidence for")
    })

    it("interprets weak evidence", {
      expect_equal(comparator$interpret_bayes_factor(2), "Weak evidence for")
    })

    it("interprets evidence against", {
      expect_equal(comparator$interpret_bayes_factor(0.05), "Strong evidence against")
    })
  })
})

describe("compare_models convenience function", {
  data <- create_test_data()
  model1 <- lmer(mean_velocity ~ rir + (1 | id), data = data)
  model2 <- lmer(mean_velocity ~ rir + (1 + rir | id), data = data)

  it("returns ComparisonResult", {
    result <- compare_models(model1, model2)
    expect_s3_class(result, "ComparisonResult")
  })
})
