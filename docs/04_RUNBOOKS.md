# Standard Operating Procedures & Runbooks

## Table of Contents

1. [Quick Reference](#quick-reference)
2. [Initial Setup](#initial-setup)
3. [Daily Development Workflow](#daily-development-workflow)
4. [Running Analyses](#running-analyses)
5. [Adding a New Study](#adding-a-new-study)
6. [Adding a New Analysis](#adding-a-new-analysis)
7. [Debugging Common Issues](#debugging-common-issues)
8. [Pre-Commit Checklist](#pre-commit-checklist)

---

## Quick Reference

```bash
# Most common commands
make setup           # First-time setup
make pipeline        # Run analysis
make test            # Run tests
make check           # Lint + test
make format          # Format code with Air
make reports         # Render Quarto reports
```

---

## Initial Setup

### Prerequisites

1. **R 4.4.2+** installed via rig:
   ```bash
   # Install rig (if not installed)
   # macOS: brew install rig
   # Linux: see https://github.com/r-lib/rig#installation

   rig install 4.4.2
   rig default 4.4.2
   ```

2. **Quarto** installed:
   ```bash
   # Download from https://quarto.org/docs/get-started/
   quarto --version  # Verify installation
   ```

3. **Pre-commit** (optional but recommended):
   ```bash
   pip install pre-commit
   ```

4. **Air formatter** (optional, much faster than styler):
   ```bash
   # If you have Rust/Cargo:
   cargo install air-r

   # Or download from: https://github.com/posit-dev/air/releases
   ```

### First-Time Setup

```bash
# 1. Clone repository
git clone <repository-url>
cd deadlift-study

# 2. Run setup (installs R packages and hooks)
make setup

# 3. Verify installation
make check
```

### Environment Configuration

Edit `.Renviron` for local settings:
```bash
# Set analysis mode (default, debug, quick, full)
R_CONFIG_ACTIVE=default

# Set seed for reproducibility
DEADLIFT_SEED=42

# Set log level (TRACE, DEBUG, INFO, WARN, ERROR)
DEADLIFT_LOG_LEVEL=INFO
```

---

## Daily Development Workflow

### Starting Work

```bash
# 1. Pull latest changes
git pull origin main

# 2. Sync dependencies (if renv.lock changed)
make install-r

# 3. Check pipeline status
make pipeline-status
```

### Development Cycle

```bash
# Edit code in R/ directory
# ...

# 4. Format code
make format

# 5. Run tests
make test

# 6. Run linting
make lint

# 7. If all passes, run pipeline
make pipeline
```

### Before Committing

```bash
# Run full check suite
make check

# Verify coverage (must be >= 90%)
make check-coverage

# Commit with descriptive message
git add .
git commit -m "feat: Add publication bias analysis"
```

---

## Running Analyses

### Full Pipeline

```bash
# Run complete analysis pipeline
make pipeline

# View pipeline visualization
make pipeline-viz
```

### Specific Analysis

```bash
# Primary analysis only
Rscript -e "targets::tar_make(names = c('model_strength', 'model_hypertrophy'))"

# Check what will run
Rscript -e "targets::tar_outdated()"
```

### Clean and Rerun

```bash
# Invalidate all targets and rerun
make pipeline-clean
make pipeline

# Invalidate specific target
Rscript -e "targets::tar_invalidate(matches('effect_sizes'))"
```

### Generate Reports

```bash
# All reports
make reports

# Specific report
make report-primary
```

---

## Adding a New Study

### Step 1: Update Raw Data

1. Add study data to `data/raw/` directory
2. Ensure columns match expected schema (see `config/default.yml` → `validation.required_columns`)

### Step 2: Create Study Configuration

Create `config/studies/<study_name>.yml`:

```yaml
<study_name>:
  inherits: default

  study:
    name: "Descriptive Study Name"
    hypothesis: "What we're testing"
    outcome_type: "strength"  # or "hypertrophy"

  data:
    filter:
      outcome_type: "strength"
```

### Step 3: Add to Pipeline

Edit `_targets.R` to include new study targets.

### Step 4: Run and Validate

```bash
# Validate data first
make validate-data

# Run pipeline
make pipeline

# Check results
make report-primary
```

---

## Adding a New Analysis

### Step 1: Create Analysis Document

```bash
# Create new Quarto document
touch analyses/05_new_analysis.qmd
```

### Step 2: Template Structure

```qmd
---
title: "New Analysis Title"
params:
  outcome: "strength"
  seed: 42
---

```{r setup}
#| include: false
box::use(
  R/utils/config[load_config],
  targets[tar_read]
)

config <- load_config()
set.seed(params$seed)
```

## Methods

Description of analysis methods...

## Results

```{r results}
# Load pre-computed results from pipeline
data <- tar_read(meta_analysis_data)
```

## Interpretation

Interpretation of findings...
```

### Step 3: Add to Quarto Project

Edit `_quarto.yml` to include new analysis in navbar.

### Step 4: Render and Review

```bash
quarto render analyses/05_new_analysis.qmd
```

---

## Debugging Common Issues

### Pipeline Fails to Start

**Symptom:** `targets::tar_make()` errors immediately

**Solutions:**
1. Check R version: `R --version` (must be 4.4.2+)
2. Restore packages: `make install-r`
3. Check for syntax errors: `make lint`

### Missing Package Dependencies

**Symptom:** `Error in library(X): there is no package called 'X'`

**Solutions:**
```bash
# Restore from lockfile
Rscript -e "renv::restore()"

# If package is new, add it
Rscript -e "renv::install('package_name')"
Rscript -e "renv::snapshot()"
```

### Effect Size Calculation Errors

**Symptom:** "Standard deviation must be positive" or similar

**Solutions:**
1. Check data validation: `make validate-data`
2. Inspect problematic rows in raw data
3. Fix data or add exclusion rule in config

### Memory Issues with Large Pipeline

**Symptom:** R session crashes or runs out of memory

**Solutions:**
```bash
# Run with garbage collection
Rscript -e "targets::tar_make(garbage_collection = TRUE)"

# Run specific targets sequentially
Rscript -e "targets::tar_make(names = target_name)"
```

### Box Module Import Errors

**Symptom:** `Error in box::use(): module 'R/...' not found`

**Solutions:**
1. Ensure working directory is project root
2. Check `.Rprofile` has `options(box.path = "R")`
3. Verify file exists at expected path

---

## Pre-Commit Checklist

Before every commit, verify:

- [ ] **Code formatted:** `make format` (or `make style` if Air not installed)
- [ ] **Linting passes:** `make lint`
- [ ] **All tests pass:** `make test`
- [ ] **Coverage >= 90%:** `make check-coverage`
- [ ] **Pipeline runs:** `make pipeline`
- [ ] **No secrets in code:** Check `.Renviron` not staged

### Quick Check Command

```bash
# Run all checks
make check

# If all passes, you're ready to commit
```

### Commit Message Format

```
<type>: <description>

[optional body]

[optional footer]
```

Types:
- `feat`: New feature
- `fix`: Bug fix
- `docs`: Documentation
- `refactor`: Code restructuring
- `test`: Adding tests
- `chore`: Maintenance tasks

Examples:
```bash
git commit -m "feat: Add Bayesian meta-analysis model"
git commit -m "fix: Correct Hedge's g calculation for small samples"
git commit -m "docs: Update runbook with troubleshooting section"
```

---

## Decision Trees

### When to Invalidate Pipeline Cache

```
Is there a new data file?
├── Yes → make pipeline-clean && make pipeline
└── No → Did you change R/ code?
    ├── Yes → Pipeline will detect changes automatically
    └── No → Did you change config?
        ├── Yes → Invalidate config target only
        └── No → No action needed
```

### When to Use Which Config Mode

```
What are you doing?
├── Debugging → R_CONFIG_ACTIVE=debug
├── Quick iteration → R_CONFIG_ACTIVE=quick
├── Final analysis → R_CONFIG_ACTIVE=full
└── Normal development → R_CONFIG_ACTIVE=default
```

---

*Last updated: 2024-12-25*
