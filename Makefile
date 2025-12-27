# Makefile
# Deadlift Study Meta-Analysis Project Automation
# Usage: make <target>

.PHONY: all setup install pipeline test lint check clean help

# Default configuration
# rig is REQUIRED for R version management (reads .R-version)
# Install rig: https://github.com/r-lib/rig
R_VERSION := $(shell cat .R-version 2>/dev/null)
WORKERS := 4

# Quarto detection: prefer system quarto, fallback to RStudio-bundled
QUARTO_SYSTEM := $(shell command -v quarto 2>/dev/null)
QUARTO_RSTUDIO := /usr/lib/rstudio/resources/app/bin/quarto/bin/quarto
QUARTO := $(if $(QUARTO_SYSTEM),$(QUARTO_SYSTEM),$(if $(wildcard $(QUARTO_RSTUDIO)),$(QUARTO_RSTUDIO),quarto))

# Verify rig is installed (fail fast if not)
RIG_CHECK := $(shell command -v rig 2>/dev/null || echo "NOT_FOUND")
ifeq ($(RIG_CHECK),NOT_FOUND)
  $(error rig is required but not installed. Install from: https://github.com/r-lib/rig)
endif

# Use rig run with explicit version from .R-version file
# Note: Commands should use $(RSCRIPT) "expression" (no -e flag needed)
RSCRIPT := rig run --r-version $(R_VERSION) --no-echo -e
RSCRIPT_FILE := rig run --r-version $(R_VERSION) -f

# Colors for output
BLUE := \033[0;34m
GREEN := \033[0;32m
YELLOW := \033[0;33m
RED := \033[0;31m
NC := \033[0m  # No Color

#------------------------------------------------------------------------------
# Help
#------------------------------------------------------------------------------
r-version:  ## Show R version being used
	@echo "$(BLUE)R Configuration:$(NC)"
	@echo "  .R-version: $(R_VERSION)"
	@echo "  Available versions (rig list):"
	@rig list
	@echo ""
	@echo "$(GREEN)Using:$(NC)"
	@$(RSCRIPT) "cat(R.version.string, '\n')"

help:  ## Show this help message
	@echo "$(BLUE)Deadlift Study Meta-Analysis$(NC)"
	@echo "=============================="
	@echo ""
	@echo "Available targets:"
	@grep -E '^[a-zA-Z_-]+:.*?## .*$$' $(MAKEFILE_LIST) | sort | \
		awk 'BEGIN {FS = ":.*?## "}; {printf "  $(GREEN)%-20s$(NC) %s\n", $$1, $$2}'
	@echo ""
	@echo "$(YELLOW)Quick start:$(NC)"
	@echo "  make setup     # First time setup"
	@echo "  make pipeline  # Run analysis"

#------------------------------------------------------------------------------
# Setup & Installation
#------------------------------------------------------------------------------
setup: install-r install-hooks  ## Complete project setup

install-r:  ## Install R packages via renv
	@echo "$(BLUE)Installing R packages...$(NC)"
	$(RSCRIPT) "if (!requireNamespace('renv', quietly = TRUE)) install.packages('renv')"
	$(RSCRIPT) "renv::restore()"
	@echo "$(GREEN)R packages installed!$(NC)"

install-hooks:  ## Install pre-commit hooks
	@echo "$(BLUE)Installing pre-commit hooks...$(NC)"
	@if command -v pre-commit >/dev/null 2>&1; then \
		pre-commit install; \
		echo "$(GREEN)Pre-commit hooks installed!$(NC)"; \
	else \
		echo "$(YELLOW)Warning: pre-commit not found. Install with: pip install pre-commit$(NC)"; \
	fi

install-quarto:  ## Check Quarto installation (auto-detects RStudio-bundled)
	@echo "$(BLUE)Checking Quarto installation...$(NC)"
	@echo "Detected QUARTO path: $(QUARTO)"
	@if [ -x "$(QUARTO)" ]; then \
		echo "$(GREEN)Quarto version: $$($(QUARTO) --version)$(NC)"; \
	else \
		echo "$(RED)Quarto not found. Please install from: https://quarto.org/docs/get-started/$(NC)"; \
	fi

#------------------------------------------------------------------------------
# Pipeline Execution
#------------------------------------------------------------------------------
all: check pipeline  ## Run checks and full pipeline

pipeline:  ## Run targets pipeline
	@echo "$(BLUE)Running analysis pipeline...$(NC)"
	$(RSCRIPT) "targets::tar_make()"
	@echo "$(GREEN)Pipeline complete!$(NC)"

pipeline-parallel:  ## Run pipeline with parallel workers
	@echo "$(BLUE)Running pipeline with $(WORKERS) workers...$(NC)"
	$(RSCRIPT) "targets::tar_make_future(workers = $(WORKERS))"

pipeline-viz:  ## Visualize pipeline DAG
	@echo "$(BLUE)Generating pipeline visualization...$(NC)"
	$(RSCRIPT) "targets::tar_visnetwork()"

pipeline-status:  ## Check pipeline status
	$(RSCRIPT) "print(targets::tar_progress())"

pipeline-clean:  ## Clean pipeline cache (invalidate all targets)
	@echo "$(YELLOW)Cleaning pipeline cache...$(NC)"
	$(RSCRIPT) "targets::tar_invalidate(everything())"
	rm -rf _targets/
	@echo "$(GREEN)Cache cleaned!$(NC)"

#------------------------------------------------------------------------------
# Testing
#------------------------------------------------------------------------------
test:  ## Run all tests
	@echo "$(BLUE)Running tests...$(NC)"
	$(RSCRIPT) "testthat::test_dir('tests/testthat')"

test-coverage:  ## Run tests with coverage report
	@echo "$(BLUE)Running tests with coverage...$(NC)"
	$(RSCRIPT) "covr::report()"

test-property:  ## Run property-based tests only
	@echo "$(BLUE)Running property-based tests...$(NC)"
	$(RSCRIPT) "testthat::test_file('tests/testthat/properties/')"

#------------------------------------------------------------------------------
# Code Quality
#------------------------------------------------------------------------------
lint:  ## Lint R code with lintr
	@echo "$(BLUE)Linting R code...$(NC)"
	$(RSCRIPT) "lintr::lint_dir('R/')"
	$(RSCRIPT) "lintr::lint_dir('tests/')"

lint-fast:  ## Fast lint with Jarl (Rust-based, if installed)
	@echo "$(BLUE)Fast linting with Jarl...$(NC)"
	@if command -v jarl >/dev/null 2>&1; then \
		jarl R/ tests/; \
	else \
		echo "$(YELLOW)Jarl not installed. Falling back to lintr...$(NC)"; \
		$(MAKE) lint; \
	fi

format:  ## Format R code with Air (Rust-based, 300x faster than styler)
	@echo "$(BLUE)Formatting R code with Air...$(NC)"
	@if command -v air >/dev/null 2>&1; then \
		air format R/; \
		air format tests/; \
		air format analyses/; \
	else \
		echo "$(YELLOW)Air not installed. Falling back to styler...$(NC)"; \
		$(MAKE) style; \
	fi

format-check:  ## Check formatting without modifying files
	@echo "$(BLUE)Checking format with Air...$(NC)"
	@if command -v air >/dev/null 2>&1; then \
		air format --check R/ tests/; \
	else \
		echo "$(YELLOW)Air not installed. Install from: https://github.com/posit-dev/air$(NC)"; \
	fi

style:  ## Style R code with styler (fallback if Air not available)
	@echo "$(BLUE)Styling R code with styler...$(NC)"
	$(RSCRIPT) "styler::style_dir('R/')"
	$(RSCRIPT) "styler::style_dir('tests/')"

type-check:  ## Run type checks with typed package
	@echo "$(BLUE)Running type checks...$(NC)"
	$(RSCRIPT) "source('scripts/type_check.R')"

check: format-check lint test  ## Run formatting, linting, and tests
	@echo "$(GREEN)All checks passed!$(NC)"

check-coverage:  ## Verify 90%+ code coverage
	@echo "$(BLUE)Checking code coverage...$(NC)"
	$(RSCRIPT) "\
		cov <- covr::package_coverage(); \
		pct <- covr::percent_coverage(cov); \
		cat(sprintf('Coverage: %.1f%%\n', pct)); \
		if (pct < 90) stop('Coverage below 90%!'); \
		cat('Coverage meets 90% threshold!\n')"

pre-commit:  ## Run pre-commit on all files
	pre-commit run --all-files

install-air:  ## Install Air formatter (Rust-based)
	@echo "$(BLUE)Installing Air formatter...$(NC)"
	@if command -v cargo >/dev/null 2>&1; then \
		cargo install air-r; \
	else \
		echo "$(YELLOW)Cargo not found. Install Rust first: https://rustup.rs$(NC)"; \
		echo "Or download pre-built binary from: https://github.com/posit-dev/air/releases"; \
	fi

#------------------------------------------------------------------------------
# Data Management
#------------------------------------------------------------------------------
validate-data:  ## Run data validation checks
	@echo "$(BLUE)Validating data...$(NC)"
	$(RSCRIPT) "source('R/validators/run_validation.R')"

dvc-push:  ## Push data to DVC remote
	dvc push

dvc-pull:  ## Pull data from DVC remote
	dvc pull

#------------------------------------------------------------------------------
# Replication Scripts
#------------------------------------------------------------------------------
replicate: replicate-pelland replicate-velocity  ## Run all replication scripts
	@echo "$(GREEN)All replications complete!$(NC)"

replicate-pelland:  ## Run Pelland meta-regression replication (original)
	@echo "$(BLUE)Running Pelland replication...$(NC)"
	$(RSCRIPT_FILE) scripts/replicate_pelland.R

replicate-pelland-refactored:  ## Run Pelland replication (OOP/SOLID)
	@echo "$(BLUE)Running refactored Pelland replication...$(NC)"
	$(RSCRIPT) "box::purge_cache()"
	$(RSCRIPT_FILE) scripts/replicate_pelland_refactored.R

replicate-velocity:  ## Run PeerJ velocity-RIR replication (original)
	@echo "$(BLUE)Running PeerJ velocity replication...$(NC)"
	$(RSCRIPT_FILE) scripts/replicate_peerj_velocity.R

replicate-velocity-refactored:  ## Run velocity replication (OOP/SOLID)
	@echo "$(BLUE)Running refactored velocity replication...$(NC)"
	$(RSCRIPT) "box::purge_cache()"
	$(RSCRIPT_FILE) scripts/replicate_peerj_velocity_refactored.R

replicate-rir-velocity:  ## Run RIR-velocity modeling replication (Jukic et al.)
	@echo "$(BLUE)Running RIR-velocity modeling replication...$(NC)"
	$(RSCRIPT) "box::purge_cache()"
	$(RSCRIPT) "source('scripts/replicate_rir_velocity_refactored.R')"

replicate-deadlift:  ## Run deadlift RIR-velocity analysis (Study 4)
	@echo "$(BLUE)Running deadlift RIR-velocity analysis...$(NC)"
	$(RSCRIPT) "box::purge_cache()"
	$(RSCRIPT) "source('scripts/replicate_deadlift_rir_velocity.R')"

replicate-all:  ## Run all replication scripts (including deadlift)
	@echo "$(BLUE)Running all replications...$(NC)"
	$(MAKE) replicate-pelland-refactored
	$(MAKE) replicate-velocity-refactored
	$(MAKE) replicate-rir-velocity
	$(MAKE) replicate-deadlift
	@echo "$(GREEN)All replications complete!$(NC)"

#------------------------------------------------------------------------------
# LMM Analysis (Study 5)
#------------------------------------------------------------------------------
analyze-lmm:  ## Run LMM analysis with velocity stop tables (Study 5)
	@echo "$(BLUE)Running LMM analysis...$(NC)"
	$(RSCRIPT) "box::purge_cache()"
	$(RSCRIPT) "source('scripts/analyze_deadlift_lmm.R')"
	@echo "$(GREEN)LMM analysis complete!$(NC)"

#------------------------------------------------------------------------------
# Advanced Velocity Analysis (Study 6)
#------------------------------------------------------------------------------
analyze-advanced:  ## Run advanced velocity analyses H2-H6 (Study 6)
	@echo "$(BLUE)Running advanced velocity analyses...$(NC)"
	$(RSCRIPT) "box::purge_cache()"
	$(RSCRIPT) "source('scripts/analyze_advanced_velocity.R')"
	@echo "$(GREEN)Advanced analyses complete!$(NC)"

analyze-all: replicate-all analyze-lmm analyze-advanced  ## Run all analyses
	@echo "$(GREEN)All analyses complete!$(NC)"

#------------------------------------------------------------------------------
# Reporting
#------------------------------------------------------------------------------
# Use R 4.4.2 for Quarto rendering (matches .R-version)
QUARTO_R := /opt/R/$(R_VERSION)/bin/R

reports:  ## Render all Quarto reports
	@echo "$(BLUE)Rendering reports...$(NC)"
	QUARTO_R=$(QUARTO_R) $(QUARTO) render analyses/

report-index:  ## Render index/overview page
	QUARTO_R=$(QUARTO_R) $(QUARTO) render analyses/index.qmd

report-pelland:  ## Render Pelland meta-regression replication
	QUARTO_R=$(QUARTO_R) $(QUARTO) render analyses/01_pelland_meta_regression.qmd

report-velocity:  ## Render velocity-RIR relationship report
	QUARTO_R=$(QUARTO_R) $(QUARTO) render analyses/02_velocity_rir_relationship.qmd

report-rir-velocity:  ## Render RIR-velocity modeling report (Jukic et al.)
	QUARTO_R=$(QUARTO_R) $(QUARTO) render analyses/03_rir_velocity_modeling.qmd

report-deadlift:  ## Render merged deadlift study (Thesis - HTML + PDF)
	QUARTO_R=$(QUARTO_R) $(QUARTO) render analyses/deadlift_study.qmd

reports-all: reports  ## Render all reports (alias for reports)
	@echo "$(GREEN)All reports rendered!$(NC)"

report-primary:  ## Render primary analysis report (legacy)
	$(QUARTO) render analyses/01_primary_analysis.qmd

report-sensitivity:  ## Render sensitivity analysis report (legacy)
	$(QUARTO) render analyses/02_sensitivity_analysis.qmd

#------------------------------------------------------------------------------
# Docker
#------------------------------------------------------------------------------
docker-build:  ## Build Docker image
	@echo "$(BLUE)Building Docker image...$(NC)"
	docker build -t deadlift-study -f docker/Dockerfile .

docker-run:  ## Run analysis in Docker container
	@echo "$(BLUE)Running in Docker...$(NC)"
	docker run --rm -v $(PWD)/data:/project/data:ro -v $(PWD)/reports:/project/reports deadlift-study

docker-dev:  ## Start development container with RStudio
	docker-compose -f docker/docker-compose.yml up dev

#------------------------------------------------------------------------------
# Cleanup
#------------------------------------------------------------------------------
clean: clean-cache clean-reports  ## Clean all generated files

clean-cache:  ## Clean R cache files
	@echo "$(YELLOW)Cleaning cache...$(NC)"
	rm -rf _targets/
	rm -rf .Rproj.user/
	rm -rf .quarto/
	find . -name "*.Rproj" -delete
	find . -name ".RData" -delete
	find . -name ".Rhistory" -delete

clean-reports:  ## Clean generated reports
	@echo "$(YELLOW)Cleaning reports...$(NC)"
	rm -rf reports/*.html
	rm -rf reports/*.pdf
	rm -rf analyses/*_cache/
	rm -rf analyses/*_files/

clean-renv:  ## Clean renv cache (use with caution)
	@echo "$(RED)Cleaning renv cache...$(NC)"
	rm -rf renv/library/
	rm -rf renv/staging/

#------------------------------------------------------------------------------
# Documentation
#------------------------------------------------------------------------------
docs:  ## Generate package documentation
	$(RSCRIPT) "devtools::document()"

readme:  ## Render README.Rmd
	$(RSCRIPT) "rmarkdown::render('README.Rmd')"
