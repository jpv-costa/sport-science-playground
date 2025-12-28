# Troubleshooting Guide

This document covers common issues and their solutions when running the deadlift study analysis pipeline.

## Table of Contents

1. [Model Convergence Issues](#model-convergence-issues)
2. [Robust LMM Failures](#robust-lmm-failures)
3. [Data Validation Errors](#data-validation-errors)
4. [Pipeline Execution Problems](#pipeline-execution-problems)
5. [Report Rendering Issues](#report-rendering-issues)

---

## Model Convergence Issues

### Full Interaction Model Fails to Converge

**Symptom**: Error message about singular fit or failure to converge when fitting the full LMM model with all interactions.

**Cause**: The full model `mean_velocity ~ rir * load_percentage * day + (1 + rir | id)` has too many parameters relative to the sample size.

**Solution**: The analysis script (`scripts/analyze_deadlift_lmm.R`) automatically falls back to simpler models:

1. First attempt: Full interaction model
2. Fallback 1: Two-way interactions only
3. Fallback 2: Main effects with random slopes
4. Fallback 3: Random intercepts only

**Manual intervention**: If all models fail, check data quality:

```r
# Check for sufficient observations per participant
table(data$id)

# Check for sufficient observations per condition
with(data, table(id, load_percentage, day))
```

### Singular Fit Warnings

**Symptom**: Warning "boundary (singular) fit: see help('isSingular')"

**Cause**: Random effects variance is estimated at or near zero.

**Impact**: Model is still usable, but random effects may not be meaningful.

**Solution**: Consider simplifying the random effects structure or pooling conditions.

---

## Robust LMM Failures

### robustlmm Package Errors

**Symptom**: Error when calling `robustlmm::rlmer()` in the robustness checks section.

**Cause**:
- Package not installed
- Convergence failure with robust estimator
- Insufficient data for robust estimation

**Solution**:

1. Install the package:
   ```r
   install.packages("robustlmm")
   ```

2. If robust LMM fails, the script proceeds with alternative robustness checks:
   - Cluster-robust standard errors (via `clubSandwich` package)
   - Bootstrap confidence intervals

3. Check the results object for the `$robustness$robust_lmm_note` field which explains any failures.

### Bootstrap Confidence Intervals Timeout

**Symptom**: Bootstrap CI calculation takes too long or hangs.

**Cause**: Default 1000 iterations with complex model.

**Solution**: The script uses 500 bootstrap iterations by default. If still slow:

```r
# Reduce iterations (in analyze_deadlift_lmm.R)
bootstrap_ci <- confint(model, method = "boot", nsim = 200)
```

---

## Data Validation Errors

### Missing Required Columns

**Symptom**: Error "Missing required columns: X, Y, Z"

**Cause**: Data loader output doesn't match expected schema.

**Solution**:

1. Check the Excel file format matches expected structure
2. Verify column names in the raw data
3. Run the validation manually:

```r
box::use(R/validators/run_validation[run_validation])
run_validation("path/to/data.xlsx")
```

### RIR Out of Range

**Symptom**: Validation error "RIR out of range [0,7]"

**Cause**: Data contains RIR values outside the expected range.

**Solution**:

1. Check raw data for data entry errors
2. RIR should be 0 (failure) to 7 (far from failure)
3. Values > 7 might indicate miscoded "reps remaining" vs "RIR"

### Invalid Load Percentages

**Symptom**: Validation error "Invalid load percentages: X"

**Cause**: Load values not matching study protocol (80% or 90%).

**Solution**:

1. Check data coding matches expected format ("80%" or "90%")
2. If using different loads, modify validation rules:

```r
# In R/validators/run_validation.R
valid_loads <- c("70%", "80%", "90%", "100%")  # Extend as needed
```

---

## Pipeline Execution Problems

### targets Not Building

**Symptom**: `make pipeline` completes but no outputs generated.

**Cause**: Targets already up-to-date or silent failures.

**Solution**:

1. Check target status:
   ```r
   targets::tar_visnetwork()  # Visualize pipeline
   targets::tar_progress()    # Check progress
   ```

2. Force rebuild:
   ```r
   targets::tar_make(callr_function = NULL)  # Run in session for debugging
   targets::tar_invalidate(everything())     # Force rebuild all
   ```

3. Check for errors in the targets log.

### Memory Issues

**Symptom**: R crashes or runs out of memory during analysis.

**Cause**: Large data or many bootstrap iterations.

**Solution**:

1. Process in chunks
2. Reduce bootstrap iterations
3. Use parallel workers more conservatively:
   ```bash
   make pipeline-parallel WORKERS=2
   ```

### renv Package Issues

**Symptom**: Packages not found or version conflicts.

**Solution**:

```r
renv::restore()           # Reinstall from lockfile
renv::status()            # Check synchronization
renv::snapshot()          # Update lockfile after changes
```

---

## Report Rendering Issues

### Quarto Not Found

**Symptom**: "quarto: command not found"

**Solution**: The Makefile auto-detects Quarto locations:

1. System Quarto: Install from https://quarto.org
2. RStudio-bundled: Available if RStudio is installed

Check detection:
```bash
make install-quarto
```

### Report Fails to Render

**Symptom**: Quarto report compilation fails.

**Cause**: Missing .rds results or R package issues.

**Solution**:

1. Ensure analysis ran first:
   ```bash
   make analyze-lmm analyze-advanced
   ```

2. Check .rds files exist:
   ```bash
   ls -la data/processed/*.rds
   ```

3. Render with verbose output:
   ```bash
   quarto render analyses/deadlift_study.qmd --log-level debug
   ```

### Images Not Displaying

**Symptom**: Missing images in rendered HTML.

**Cause**: Git LFS tracking HTML embedded images.

**Solution**: The `.gitattributes` file is configured to NOT use LFS for `docs/*_files/`. Verify:

```bash
git check-attr filter docs/deadlift_study_files/figure-html/plot.png
# Should return: unspecified
```

---

## Getting Help

If issues persist:

1. Check the GitHub Issues: https://github.com/anthropics/claude-code/issues
2. Ensure R version matches `.R-version` (4.4.2)
3. Verify all dependencies: `make check`
4. Run tests to identify broken components: `make test`

## Common Environment Setup Checklist

```bash
# 1. R version management
rig list                    # Available R versions
rig default 4.4.2          # Set correct version

# 2. Package installation
make install-r             # renv restore

# 3. Quarto check
make install-quarto        # Verify installation

# 4. Run full check
make check                 # Format + lint + test

# 5. Run pipeline
make analyze-all           # All analyses
make reports               # All reports
```
