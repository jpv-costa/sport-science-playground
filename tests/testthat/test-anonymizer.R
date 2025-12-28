# tests/testthat/test-anonymizer.R
# Unit tests for Anonymizer

box::use(
  testthat[...],
  ../../R/utils/anonymizer[Anonymizer]
)

describe("Anonymizer", {

  describe("initialization", {
    it("creates instance with mapping", {
      anon <- Anonymizer$new()

      expect_s3_class(anon, "Anonymizer")
    })

    it("has 19 participants in mapping", {
      anon <- Anonymizer$new()

      mapping <- anon$get_mapping()

      expect_equal(nrow(mapping), 19)
    })
  })

  describe("anonymize_id", {
    it("maps first participant to P01", {
      anon <- Anonymizer$new()

      result <- anon$anonymize_id("Alexia")

      expect_equal(result, "P01")
    })

    it("maps last participant to P19", {
      anon <- Anonymizer$new()

      result <- anon$anonymize_id("Toni")

      expect_equal(result, "P19")
    })

    it("maps Sebastião to P17", {
      anon <- Anonymizer$new()

      result <- anon$anonymize_id("Sebastião")

      expect_equal(result, "P17")
    })

    it("handles vector of IDs", {
      anon <- Anonymizer$new()

      result <- anon$anonymize_id(c("Alexia", "Toni", "Sebastião"))

      expect_equal(result, c("P01", "P19", "P17"))
    })

    it("handles empty input", {
      anon <- Anonymizer$new()

      result <- anon$anonymize_id(character(0))

      expect_length(result, 0)
    })

    it("throws error for unknown ID", {
      anon <- Anonymizer$new()

      expect_error(
        anon$anonymize_id("Unknown"),
        "Unknown participant ID"
      )
    })
  })

  describe("anonymize_set_id", {
    it("replaces name in set_id format", {
      anon <- Anonymizer$new()

      result <- anon$anonymize_set_id("Sebastião_Day1_90pct_S1")

      expect_equal(result, "P17_Day1_90pct_S1")
    })

    it("handles complex set_id with spaces in name", {
      anon <- Anonymizer$new()

      # Note: Manuel M. contains a space and period
      result <- anon$anonymize_set_id("Manuel M._Day2_80pct_S2")

      expect_equal(result, "P05_Day2_80pct_S2")
    })

    it("handles vector of set_ids", {
      anon <- Anonymizer$new()

      result <- anon$anonymize_set_id(c(
        "Alexia_Day1_90pct_S1",
        "Toni_Day2_80pct_S3"
      ))

      expect_equal(result, c("P01_Day1_90pct_S1", "P19_Day2_80pct_S3"))
    })

    it("handles empty input", {
      anon <- Anonymizer$new()

      result <- anon$anonymize_set_id(character(0))

      expect_length(result, 0)
    })

    it("throws error for invalid format", {
      anon <- Anonymizer$new()

      expect_error(
        anon$anonymize_set_id("InvalidFormat"),
        "Invalid set_id format"
      )
    })
  })

  describe("anonymize_data", {
    it("anonymizes id column in data frame", {
      anon <- Anonymizer$new()
      data <- data.frame(
        id = c("Alexia", "Toni"),
        velocity = c(0.5, 0.6),
        stringsAsFactors = FALSE
      )

      result <- anon$anonymize_data(data, id_col = "id", set_id_col = "nonexistent")

      expect_equal(result$id, c("P01", "P19"))
      expect_equal(result$velocity, c(0.5, 0.6))
    })

    it("anonymizes both id and set_id columns", {
      anon <- Anonymizer$new()
      data <- data.frame(
        id = c("Sebastião", "Mário"),
        set_id = c("Sebastião_Day1_90pct_S1", "Mário_Day2_80pct_S2"),
        velocity = c(0.5, 0.6),
        stringsAsFactors = FALSE
      )

      result <- anon$anonymize_data(data)

      expect_equal(result$id, c("P17", "P06"))
      expect_equal(result$set_id, c("P17_Day1_90pct_S1", "P06_Day2_80pct_S2"))
    })

    it("throws error for non-data-frame input", {
      anon <- Anonymizer$new()

      expect_error(
        anon$anonymize_data("not a data frame"),
        "Input must be a data frame"
      )
    })
  })

  describe("get_mapping", {
    it("returns data frame with original and pseudonym columns", {
      anon <- Anonymizer$new()

      mapping <- anon$get_mapping()

      expect_true("original" %in% names(mapping))
      expect_true("pseudonym" %in% names(mapping))
    })

    it("has pseudonyms in P01-P19 format", {
      anon <- Anonymizer$new()

      mapping <- anon$get_mapping()

      expect_true(all(grepl("^P[0-9]{2}$", mapping$pseudonym)))
    })

    it("pseudonyms are sorted P01 to P19", {
      anon <- Anonymizer$new()

      mapping <- anon$get_mapping()

      expected <- sprintf("P%02d", 1:19)
      expect_equal(mapping$pseudonym, expected)
    })
  })

  describe("deanonymize_id", {
    it("reverses P01 to Alexia", {
      anon <- Anonymizer$new()

      result <- anon$deanonymize_id("P01")

      expect_equal(result, "Alexia")
    })

    it("reverses P17 to Sebastião", {
      anon <- Anonymizer$new()

      result <- anon$deanonymize_id("P17")

      expect_equal(result, "Sebastião")
    })

    it("handles vector of pseudonyms", {
      anon <- Anonymizer$new()

      result <- anon$deanonymize_id(c("P01", "P17", "P19"))

      expect_equal(result, c("Alexia", "Sebastião", "Toni"))
    })

    it("throws error for unknown pseudonym", {
      anon <- Anonymizer$new()

      expect_error(
        anon$deanonymize_id("P99"),
        "Unknown pseudonym"
      )
    })
  })

  describe("is_anonymized", {
    it("returns TRUE for valid pseudonym format", {
      anon <- Anonymizer$new()

      expect_true(anon$is_anonymized("P01"))
      expect_true(anon$is_anonymized("P19"))
    })

    it("returns FALSE for original names", {
      anon <- Anonymizer$new()

      expect_false(anon$is_anonymized("Sebastião"))
      expect_false(anon$is_anonymized("Alexia"))
    })

    it("returns FALSE for invalid formats", {
      anon <- Anonymizer$new()

      expect_false(anon$is_anonymized("P1"))  # Too short
      expect_false(anon$is_anonymized("P001"))  # Too long
      expect_false(anon$is_anonymized("A01"))  # Wrong prefix
    })
  })

  describe("round-trip consistency", {
    it("anonymize then deanonymize returns original", {
      anon <- Anonymizer$new()
      original <- "Sebastião"

      pseudonym <- anon$anonymize_id(original)
      recovered <- anon$deanonymize_id(pseudonym)

      expect_equal(recovered, original)
    })

    it("all 19 participants round-trip correctly", {
      anon <- Anonymizer$new()
      mapping <- anon$get_mapping()

      for (i in seq_len(nrow(mapping))) {
        original <- mapping$original[i]
        pseudonym <- anon$anonymize_id(original)
        recovered <- anon$deanonymize_id(pseudonym)

        expect_equal(recovered, original,
                     info = sprintf("Failed for %s", original))
      }
    })
  })
})
