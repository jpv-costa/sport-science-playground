# tests/testthat.R
# Test runner for Deadlift Study

library(testthat)

# Set box path for module imports
options(box.path = file.path(getwd(), "R"))

# Run all tests
test_check("deadlift.study")
