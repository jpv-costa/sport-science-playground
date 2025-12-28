# tests/testthat/test-statistical_methods_guide.R
# Tests for StatisticalMethodsGuide R6 class
# Following CLAUDE.md testing principles: TDD, classicist approach

box::use(
  testthat[describe, it, expect_true, expect_s3_class, expect_type,
           expect_length, expect_gt]
)

# Set box path relative to project root
options(box.path = c("../../R", getOption("box.path")))

# Load the module under test
box::use(
  guides/statistical_methods_guide[
    StatisticalMethodsGuide,
    print_statistical_guide,
    get_guide_topics
  ]
)

describe("StatisticalMethodsGuide", {

  describe("initialization", {
    it("creates instance", {
      guide <- StatisticalMethodsGuide$new()
      expect_s3_class(guide, "StatisticalMethodsGuide")
    })
  })

  describe("explain_why_mixed_models", {
    guide <- StatisticalMethodsGuide$new()

    it("returns character vector", {
      result <- guide$explain_why_mixed_models()
      expect_type(result, "character")
    })

    it("includes key concepts", {
      result <- guide$explain_why_mixed_models()
      content <- paste(result, collapse = " ")
      expect_true(grepl("Non-Independent", content, ignore.case = TRUE))
      expect_true(grepl("random", content, ignore.case = TRUE))
      expect_true(grepl("intercept", content, ignore.case = TRUE))
    })

    it("has substantial content", {
      result <- guide$explain_why_mixed_models()
      expect_gt(length(result), 10)
    })
  })

  describe("explain_random_effects_structure", {
    guide <- StatisticalMethodsGuide$new()

    it("returns character vector", {
      result <- guide$explain_random_effects_structure()
      expect_type(result, "character")
    })

    it("includes three-step framework", {
      result <- guide$explain_random_effects_structure()
      content <- paste(result, collapse = " ")
      expect_true(grepl("Step 1", content))
      expect_true(grepl("Step 2", content))
      expect_true(grepl("Step 3", content))
    })
  })

  describe("explain_bayes_factors", {
    guide <- StatisticalMethodsGuide$new()

    it("returns character vector", {
      result <- guide$explain_bayes_factors()
      expect_type(result, "character")
    })

    it("includes Jeffreys scale", {
      result <- guide$explain_bayes_factors()
      content <- paste(result, collapse = " ")
      expect_true(grepl("Jeffreys", content))
      expect_true(grepl("decisive", content, ignore.case = TRUE))
      expect_true(grepl("strong", content, ignore.case = TRUE))
    })

    it("explains BIC approximation", {
      result <- guide$explain_bayes_factors()
      content <- paste(result, collapse = " ")
      expect_true(grepl("BIC", content))
    })
  })

  describe("explain_variance_components", {
    guide <- StatisticalMethodsGuide$new()

    it("returns character vector", {
      result <- guide$explain_variance_components()
      expect_type(result, "character")
    })

    it("includes ICC interpretation", {
      result <- guide$explain_variance_components()
      content <- paste(result, collapse = " ")
      expect_true(grepl("ICC", content))
      expect_true(grepl("0.05", content))
      expect_true(grepl("0.25", content))
    })

    it("explains design effect", {
      result <- guide$explain_variance_components()
      content <- paste(result, collapse = " ")
      expect_true(grepl("design effect", content, ignore.case = TRUE))
    })
  })

  describe("explain_visual_diagnostics", {
    guide <- StatisticalMethodsGuide$new()

    it("returns character vector", {
      result <- guide$explain_visual_diagnostics()
      expect_type(result, "character")
    })

    it("mentions QQ plot", {
      result <- guide$explain_visual_diagnostics()
      content <- paste(result, collapse = " ")
      expect_true(grepl("QQ", content))
    })

    it("emphasizes visual over tests", {
      result <- guide$explain_visual_diagnostics()
      content <- paste(result, collapse = " ")
      expect_true(grepl("visual", content, ignore.case = TRUE))
    })
  })

  describe("explain_robustness_checks", {
    guide <- StatisticalMethodsGuide$new()

    it("returns character vector", {
      result <- guide$explain_robustness_checks()
      expect_type(result, "character")
    })

    it("includes pit stop philosophy", {
      result <- guide$explain_robustness_checks()
      content <- paste(result, collapse = " ")
      expect_true(grepl("pit stop", content, ignore.case = TRUE))
    })

    it("covers three key checks", {
      result <- guide$explain_robustness_checks()
      content <- paste(result, collapse = " ")
      expect_true(grepl("cluster-robust", content, ignore.case = TRUE))
      expect_true(grepl("bootstrap", content, ignore.case = TRUE))
      expect_true(grepl("sensitivity", content, ignore.case = TRUE))
    })
  })

  describe("explain_effect_sizes", {
    guide <- StatisticalMethodsGuide$new()

    it("returns character vector", {
      result <- guide$explain_effect_sizes()
      expect_type(result, "character")
    })

    it("includes Cohen's d guidelines", {
      result <- guide$explain_effect_sizes()
      content <- paste(result, collapse = " ")
      expect_true(grepl("Cohen", content))
      expect_true(grepl("0.2", content))
      expect_true(grepl("0.5", content))
      expect_true(grepl("0.8", content))
    })

    it("explains R-squared types", {
      result <- guide$explain_effect_sizes()
      content <- paste(result, collapse = " ")
      expect_true(grepl("marginal", content, ignore.case = TRUE))
      expect_true(grepl("conditional", content, ignore.case = TRUE))
    })
  })

  describe("get_all_topics", {
    guide <- StatisticalMethodsGuide$new()

    it("returns named list", {
      result <- guide$get_all_topics()
      expect_type(result, "list")
      expect_true(length(names(result)) > 0)
    })

    it("contains all expected topics", {
      result <- guide$get_all_topics()
      expected_topics <- c(
        "why_mixed_models",
        "random_effects_structure",
        "bayes_factors",
        "variance_components",
        "visual_diagnostics",
        "robustness_checks",
        "effect_sizes"
      )
      expect_true(all(expected_topics %in% names(result)))
    })
  })

  describe("get_quick_reference", {
    guide <- StatisticalMethodsGuide$new()

    it("returns character vector", {
      result <- guide$get_quick_reference()
      expect_type(result, "character")
    })

    it("is concise but comprehensive", {
      result <- guide$get_quick_reference()
      # Should be shorter than full explanations but still useful
      expect_gt(length(result), 20)
      content <- paste(result, collapse = " ")
      expect_true(grepl("ICC", content))
      expect_true(grepl("BAYES", content, ignore.case = TRUE))
      expect_true(grepl("Cohen", content))
    })
  })
})

describe("get_guide_topics", {
  it("returns character vector of topics", {
    topics <- get_guide_topics()
    expect_type(topics, "character")
    expect_length(topics, 7)
  })

  it("includes all expected topics", {
    topics <- get_guide_topics()
    expect_true("bayes_factors" %in% topics)
    expect_true("robustness_checks" %in% topics)
    expect_true("effect_sizes" %in% topics)
  })
})
