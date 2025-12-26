# tests/testthat/test-correlation_analyzer.R
# Unit tests for CorrelationAnalyzer

box::use(
 testthat[...],
 ../../R/calculators/correlation_analyzer[
   CorrelationAnalyzer,
   CorrelationResult,
   CorrelationSummary
 ]
)

describe("CorrelationResult", {

 describe("initialization", {
   it("stores correlation and sample size", {
     result <- CorrelationResult$new(r = 0.75, n = 100)

     expect_equal(result$r, 0.75)
     expect_equal(result$n, 100)
     expect_equal(result$r_squared, 0.5625)
   })

   it("stores optional group_id", {
     result <- CorrelationResult$new(r = 0.5, n = 50, group_id = "bench")

     expect_equal(result$group_id, "bench")
   })
 })

 describe("to_list", {
   it("returns all fields as list", {
     result <- CorrelationResult$new(r = 0.6, n = 80, group_id = "squat")
     as_list <- result$to_list()

     expect_equal(as_list$r, 0.6)
     expect_equal(as_list$r_squared, 0.36)
     expect_equal(as_list$n, 80)
     expect_equal(as_list$group_id, "squat")
   })
 })
})

describe("CorrelationSummary", {

 describe("initialization", {
   it("calculates summary statistics from correlations", {
     correlations <- c(0.3, 0.5, 0.7, 0.6, 0.4)
     summary <- CorrelationSummary$new(correlations)

     expect_equal(summary$n_groups, 5)
     expect_equal(summary$mean_r, 0.5)
     expect_equal(summary$min_r, 0.3)
     expect_equal(summary$max_r, 0.7)
   })

   it("handles single correlation", {
     summary <- CorrelationSummary$new(c(0.8))

     expect_equal(summary$n_groups, 1)
     expect_equal(summary$mean_r, 0.8)
     expect_true(is.na(summary$sd_r))  # SD undefined for n=1
   })
 })

 describe("to_list", {
   it("returns summary as list", {
     correlations <- c(0.4, 0.6)
     summary <- CorrelationSummary$new(correlations)
     as_list <- summary$to_list()

     expect_equal(as_list$n_groups, 2)
     expect_equal(as_list$mean_r, 0.5)
   })
 })
})

describe("CorrelationAnalyzer", {

 describe("calculate_overall", {
   it("computes Pearson correlation correctly", {
     analyzer <- CorrelationAnalyzer$new(method = "pearson")

     x <- c(1, 2, 3, 4, 5)
     y <- c(2, 4, 6, 8, 10)

     result <- analyzer$calculate_overall(x, y)

     expect_equal(result$r, 1.0, tolerance = 0.0001)
     expect_equal(result$n, 5)
   })

   it("handles negative correlations", {
     analyzer <- CorrelationAnalyzer$new()

     x <- c(1, 2, 3, 4, 5)
     y <- c(10, 8, 6, 4, 2)

     result <- analyzer$calculate_overall(x, y)

     expect_equal(result$r, -1.0, tolerance = 0.0001)
   })

   it("handles NA values", {
     analyzer <- CorrelationAnalyzer$new()

     x <- c(1, 2, NA, 4, 5)
     y <- c(2, 4, 6, 8, NA)

     result <- analyzer$calculate_overall(x, y)

     expect_equal(result$n, 3)  # Only 3 complete cases
   })
 })

 describe("calculate_by_group", {
   it("computes correlations for each group", {
     analyzer <- CorrelationAnalyzer$new()

     data <- data.frame(
       id = rep(c("A", "B"), each = 5),
       velocity = c(0.3, 0.4, 0.5, 0.6, 0.7,   # Group A
                    0.2, 0.3, 0.4, 0.5, 0.6),  # Group B
       rir = c(1, 2, 3, 4, 5,                   # Group A: perfect correlation
               5, 4, 3, 2, 1)                   # Group B: negative correlation
     )

     result <- analyzer$calculate_by_group(data, "velocity", "rir", "id")

     expect_equal(result$n_groups, 2)
     # Group A should have r = 1, Group B should have r = -1
     correlations <- result$correlations
     expect_true(any(abs(correlations - 1) < 0.0001))
     expect_true(any(abs(correlations + 1) < 0.0001))
   })
 })

 describe("calculate_by_subgroup", {
   it("returns named list of results by subgroup", {
     analyzer <- CorrelationAnalyzer$new()

     data <- data.frame(
       exercise = rep(c("Squat", "Bench"), each = 10),
       velocity = runif(20, 0.2, 0.8),
       rir = runif(20, 0, 5)
     )

     results <- analyzer$calculate_by_subgroup(data, "velocity", "rir", "exercise")

     expect_true("Squat" %in% names(results))
     expect_true("Bench" %in% names(results))
     expect_s3_class(results$Squat, "CorrelationResult")
     expect_s3_class(results$Bench, "CorrelationResult")
   })
 })

 describe("analyze_velocity_rir", {
   it("performs full analysis returning all components", {
     analyzer <- CorrelationAnalyzer$new()

     # Create realistic test data
     set.seed(42)
     data <- data.frame(
       id = rep(1:5, each = 20),
       exercise = rep(c("Squat", "Bench"), times = 50),
       mean_velocity = runif(100, 0.2, 0.7),
       perceived_rir = runif(100, 0, 5)
     )

     results <- analyzer$analyze_velocity_rir(data)

     expect_true("overall" %in% names(results))
     expect_true("by_participant" %in% names(results))
     expect_true("by_exercise" %in% names(results))

     expect_s3_class(results$overall, "CorrelationResult")
     expect_s3_class(results$by_participant, "CorrelationSummary")
     expect_type(results$by_exercise, "list")
   })
 })
})
