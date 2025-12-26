# .Rprofile
# Project-specific R configuration for Deadlift Study Meta-Analysis

# Set CRAN mirror to Posit Public Package Manager for reproducibility
# Locked to specific date for consistent package versions
options(repos = c(CRAN = "https://packagemanager.posit.co/cran/2024-12-01"))

# Load renv for package management (if available)
if (file.exists("renv/activate.R")) {
  source("renv/activate.R")
}

# Box module options for R6 class organization
options(box.path = "R")

# Logging configuration
options(lgr.default_threshold = "info")

# Strict mode: catch partial matching issues early
options(
  warnPartialMatchArgs = TRUE,
  warnPartialMatchAttr = TRUE,
  warnPartialMatchDollar = TRUE
)

# Scientific notation threshold
options(scipen = 10)

# Default seed for reproducibility (can be overridden in analysis)
options(deadlift.default_seed = 42L)

# Development helpers (only in interactive sessions)
if (interactive()) {
  # Load devtools for development workflow
  if (requireNamespace("devtools", quietly = TRUE)) {
    suppressMessages(require(devtools))
  }

  # Load usethis for project setup helpers
  if (requireNamespace("usethis", quietly = TRUE)) {
    suppressMessages(require(usethis))
  }

  # Set prompt with project name for clarity
  options(prompt = "[deadlift-study] > ")

  # Colorful output
  options(crayon.enabled = TRUE)

  # Show more output in console
  options(max.print = 200)
}

# Print confirmation message
cat("Loaded .Rprofile for deadlift-study project\n")
cat("  - CRAN mirror: packagemanager.posit.co (2024-12-01 snapshot)\n")
cat("  - Box module path: R/\n")
cat("  - Default seed: 42\n")
