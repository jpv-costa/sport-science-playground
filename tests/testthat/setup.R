# tests/testthat/setup.R
# Test fixtures and setup

# Set box path for module imports
options(box.path = file.path(getwd(), "R"))

# Suppress logging during tests
options(lgr.default_threshold = "error")

# Create a sample TreatmentGroup for testing
create_sample_treatment_group <- function(
    group_id = "test_group_1",
    rir = 3,
    mean_pre = 100,
    mean_post = 110,
    sd_pre = 15,
    n = 30,
    correlation = 0.5
) {
  box::use(../../R/domain/treatment_group[TreatmentGroup])

  TreatmentGroup$new(
    group_id = group_id,
    repetitions_in_reserve = rir,
    mean_pre = mean_pre,
    mean_post = mean_post,
    standard_deviation_pre = sd_pre,
    sample_size = n,
    pre_post_correlation = correlation
  )
}

# Create a sample Study for testing
create_sample_study <- function(
    study_id = "test_study_2023",
    author = "Smith",
    year = 2023,
    outcome = "strength"
) {
  box::use(../../R/domain/study[Study])

  Study$new(
    study_id = study_id,
    first_author = author,
    publication_year = year,
    outcome_category = outcome
  )
}
