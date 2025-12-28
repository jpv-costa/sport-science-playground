# tests/testthat/test-glm_family_guide.R
# Unit tests for GlmFamilyGuide and GlmFamilyInfo

box::use(
  testthat[describe, it, expect_equal, expect_true, expect_false,
           expect_s3_class, expect_error, expect_type, expect_length],
  ../../R/guides/glm_family_guide[
    GlmFamilyInfo,
    GlmFamilyGuide
  ]
)

# =============================================================================
# GLM FAMILY INFO TESTS
# =============================================================================

describe("GlmFamilyInfo", {

  describe("initialization", {

    it("stores all fields correctly", {
      info <- GlmFamilyInfo$new(
        name = "Test",
        link = "identity",
        outcome_type = "Continuous",
        r_function = "glm()",
        interpretation = "Direct effect",
        assumptions = c("Normality", "Independence"),
        example = "Test example"
      )

      expect_equal(info$name, "Test")
      expect_equal(info$link, "identity")
      expect_equal(info$outcome_type, "Continuous")
      expect_equal(info$r_function, "glm()")
      expect_equal(info$interpretation, "Direct effect")
      expect_length(info$assumptions, 2)
      expect_equal(info$example, "Test example")
    })
  })

  describe("to_list", {

    it("returns complete list representation", {
      info <- GlmFamilyInfo$new(
        name = "Binomial",
        link = "logit",
        outcome_type = "Binary",
        r_function = "glm(family = binomial())",
        interpretation = "Odds ratios",
        assumptions = c("Binary outcome"),
        example = "Success/failure"
      )

      result <- info$to_list()

      expect_type(result, "list")
      expect_true("name" %in% names(result))
      expect_true("link" %in% names(result))
      expect_true("assumptions" %in% names(result))
    })
  })
})

# =============================================================================
# GLM FAMILY GUIDE TESTS
# =============================================================================

describe("GlmFamilyGuide", {

  describe("initialization", {

    it("creates guide without error", {
      guide <- GlmFamilyGuide$new()
      expect_s3_class(guide, "GlmFamilyGuide")
    })
  })

  describe("list_families", {

    it("returns all available families", {
      guide <- GlmFamilyGuide$new()
      families <- guide$list_families()

      expect_type(families, "character")
      expect_true("gaussian" %in% families)
      expect_true("binomial" %in% families)
      expect_true("poisson" %in% families)
      expect_true("gamma" %in% families)
    })
  })

  describe("describe_family", {

    it("returns GlmFamilyInfo for valid family", {
      guide <- GlmFamilyGuide$new()

      result <- guide$describe_family("binomial")

      expect_s3_class(result, "GlmFamilyInfo")
      expect_equal(result$name, "Binomial")
      expect_equal(result$link, "logit")
    })

    it("is case insensitive", {
      guide <- GlmFamilyGuide$new()

      result1 <- guide$describe_family("GAMMA")
      result2 <- guide$describe_family("gamma")

      expect_equal(result1$name, result2$name)
    })

    it("throws error for unknown family", {
      guide <- GlmFamilyGuide$new()

      expect_error(
        guide$describe_family("unknown_family"),
        "Unknown family"
      )
    })
  })

  describe("suggest_family", {

    it("suggests binomial for binary outcomes", {
      guide <- GlmFamilyGuide$new()

      result <- guide$suggest_family("binary")

      expect_s3_class(result, "GlmFamilyInfo")
      expect_equal(result$name, "Binomial")
    })

    it("suggests poisson for count data with mean = variance", {
      guide <- GlmFamilyGuide$new()

      result <- guide$suggest_family("count", mean_equals_variance = TRUE)

      expect_equal(result$name, "Poisson")
    })

    it("suggests negative binomial for overdispersed counts", {
      guide <- GlmFamilyGuide$new()

      result <- guide$suggest_family("count", mean_equals_variance = FALSE)

      expect_equal(result$name, "Negative Binomial")
    })

    it("suggests gaussian for normal continuous data", {
      guide <- GlmFamilyGuide$new()

      result <- guide$suggest_family("continuous_normal")

      expect_equal(result$name, "Gaussian")
    })

    it("suggests gamma for skewed continuous data", {
      guide <- GlmFamilyGuide$new()

      result <- guide$suggest_family("continuous_skewed")

      expect_equal(result$name, "Gamma")
    })

    it("is case insensitive", {
      guide <- GlmFamilyGuide$new()

      result1 <- guide$suggest_family("BINARY")
      result2 <- guide$suggest_family("binary")

      expect_equal(result1$name, result2$name)
    })

    it("throws error for unknown outcome type", {
      guide <- GlmFamilyGuide$new()

      expect_error(
        guide$suggest_family("invalid_type"),
        "Unknown outcome type"
      )
    })
  })

  describe("analyze_outcome", {

    it("detects binary outcomes from 0/1 vector", {
      guide <- GlmFamilyGuide$new()
      y <- c(0, 1, 0, 1, 1, 0, 1, 0, 0, 1)

      result <- guide$analyze_outcome(y)

      expect_equal(result$detected_type, "binary")
      expect_equal(result$suggestion$name, "Binomial")
    })

    it("detects count data", {
      guide <- GlmFamilyGuide$new()
      set.seed(42)
      y <- rpois(100, lambda = 5)

      result <- guide$analyze_outcome(y)

      expect_equal(result$detected_type, "count")
      expect_true(result$suggestion$name %in% c("Poisson", "Negative Binomial"))
    })

    it("detects overdispersed counts", {
      guide <- GlmFamilyGuide$new()
      set.seed(42)
      # Create overdispersed data (variance >> mean)
      y <- rnbinom(100, size = 2, mu = 5)

      result <- guide$analyze_outcome(y)

      expect_equal(result$detected_type, "count")
      expect_true(!is.null(result$dispersion_ratio))
    })

    it("detects zero-inflated data", {
      guide <- GlmFamilyGuide$new()
      # Create data with many zeros
      y <- c(rep(0, 60), rpois(40, lambda = 3))

      result <- guide$analyze_outcome(y)

      expect_equal(result$detected_type, "count")
      expect_true(result$zero_proportion > 0.5)
      expect_true(grepl("zero", result$note, ignore.case = TRUE))
    })

    it("detects continuous normal data", {
      guide <- GlmFamilyGuide$new()
      set.seed(42)
      y <- rnorm(100, mean = 0.3, sd = 0.05)

      result <- guide$analyze_outcome(y)

      expect_equal(result$detected_type, "continuous_normal")
    })

    it("detects skewed continuous data", {
      guide <- GlmFamilyGuide$new()
      set.seed(42)
      # Create positively skewed data (gamma distributed)
      y <- rgamma(100, shape = 2, rate = 10)

      result <- guide$analyze_outcome(y)

      expect_equal(result$detected_type, "continuous_skewed")
      expect_true(result$skewness > 0)
    })

    it("handles factor outcomes", {
      guide <- GlmFamilyGuide$new()
      y <- factor(c("yes", "no", "yes", "no", "yes"))

      result <- guide$analyze_outcome(y)

      expect_equal(result$detected_type, "binary")
    })

    it("handles ordered factors as ordinal", {
      guide <- GlmFamilyGuide$new()
      y <- ordered(c("low", "medium", "high", "medium", "low"),
                   levels = c("low", "medium", "high"))

      result <- guide$analyze_outcome(y)

      expect_equal(result$detected_type, "ordinal")
    })
  })

  describe("get_interpretation_guide", {

    it("returns guide for gaussian family", {
      guide <- GlmFamilyGuide$new()

      result <- guide$get_interpretation_guide("gaussian")

      expect_type(result, "character")
      expect_true(grepl("Identity", result))
    })

    it("returns guide for binomial family", {
      guide <- GlmFamilyGuide$new()

      result <- guide$get_interpretation_guide("binomial")

      expect_true(grepl("odds ratio", result, ignore.case = TRUE))
    })

    it("returns guide for poisson family", {
      guide <- GlmFamilyGuide$new()

      result <- guide$get_interpretation_guide("poisson")

      expect_true(grepl("Exponentiate", result))
      expect_true(grepl("rate ratio", result, ignore.case = TRUE))
    })

    it("handles unknown family gracefully", {
      guide <- GlmFamilyGuide$new()

      result <- guide$get_interpretation_guide("unknown")

      expect_true(grepl("No interpretation guide available", result))
    })
  })

  describe("VBT-specific recommendations", {

    it("recommends appropriate family for velocity data", {
      guide <- GlmFamilyGuide$new()
      # Simulated velocity data (positive, possibly skewed)
      set.seed(42)
      velocity <- abs(rnorm(100, mean = 0.3, sd = 0.05))

      result <- guide$analyze_outcome(velocity)

      # Should be continuous (normal or skewed depending on actual distribution)
      expect_true(result$detected_type %in%
                    c("continuous_normal", "continuous_skewed"))
    })

    it("recommends poisson for rep counts", {
      guide <- GlmFamilyGuide$new()
      # Simulated rep counts
      set.seed(42)
      reps <- rpois(100, lambda = 8)

      result <- guide$analyze_outcome(reps)

      expect_equal(result$detected_type, "count")
    })
  })
})
