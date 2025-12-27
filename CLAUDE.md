# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

Meta-analysis examining the relationship between Repetitions in Reserve (RIR) and resistance training outcomes. Includes replications of published studies (Pelland, Paulsen, Jukic) and original deadlift research.

## Commands

### Prerequisites
- R 4.4.2 (managed via `rig`)
- Quarto (auto-detects RStudio-bundled version)

### Essential Commands

```bash
# Run R code (uses rig for version management)
rig run --r-version 4.4.2 --no-echo -e "source('scripts/analyze_deadlift_lmm.R')"

# Run analysis scripts
make analyze-lmm          # Study 5: LMM analysis
make analyze-advanced     # Study 6: H2-H6 analyses

# Render reports
make report-lmm           # Single report
make reports              # All Quarto reports

# Run tests
make test                 # All tests
rig run --r-version 4.4.2 --no-echo -e "testthat::test_file('tests/testthat/test-specific.R')"

# Code quality
make check                # Format + lint + test
make format               # Format with Air (Rust-based)
make lint                 # Lint with lintr
```

## Architecture

### Module System (box)

Uses [box](https://klmr.me/box/) for Python-like explicit imports. Modules are in `R/`:

```r
box::use(
  R/loaders/deadlift_rir_data_loader[DeadliftRirDataLoader],
  R/calculators/lmm_analyzer[LmmAnalyzer]
)
```

The `.Rprofile` sets `options(box.path = "R")` so modules resolve from project root.

### R6 Class Structure

```
R/
├── loaders/              # Data loading (DeadliftRirDataLoader)
├── calculators/          # Analysis logic (LmmAnalyzer, AdvancedVelocityAnalyzer)
├── domain/               # Business entities (Study, TreatmentGroup, EffectSize)
├── models/               # RIR-velocity modeling (RirVelocityModeler)
└── utils/                # Configuration, logging, seeds
```

### Analysis Pipeline

1. **Data Loading**: `R/loaders/` - R6 classes that parse raw Excel/CSV
2. **Analysis**: `R/calculators/` - R6 classes with analysis methods
3. **Scripts**: `scripts/` - Orchestration scripts that use loaders and calculators
4. **Reports**: `analyses/*.qmd` - Quarto reports that load saved `.rds` results

Scripts save results to `data/processed/*.rds`, reports load from those files.

### Studies

| Study | Script | Report | Focus |
|-------|--------|--------|-------|
| 1 | `replicate_pelland.R` | `01_pelland_meta_regression.qmd` | Meta-regression |
| 2 | `replicate_peerj_velocity.R` | `02_velocity_rir_relationship.qmd` | Velocity-RIR |
| 3 | `replicate_rir_velocity_refactored.R` | `03_rir_velocity_modeling.qmd` | Squat modeling |
| 4 | `replicate_deadlift_rir_velocity.R` | `04_deadlift_rir_velocity.qmd` | Deadlift original |
| 5 | `analyze_deadlift_lmm.R` | `05_deadlift_lmm_analysis.qmd` | LMM & velocity tables |
| 6 | `analyze_advanced_velocity.R` | `06_advanced_velocity_analyses.qmd` | MVT, reliability, decay |

## Key Patterns

### Creating an R6 Analyzer Class

```r
box::use(
  R6[R6Class],
  stats[lm, anova, coef]
)

MyAnalyzer <- R6Class(
  classname = "MyAnalyzer",
  public = list(
    analyze = function(data) {
      # Return result object
    }
  ),
  private = list(
    .helper = function(x) { }
  )
)
```

### Data Loader Pattern

Loaders parse specific data formats and return standardized data frames:

```r
loader <- DeadliftRirDataLoader$new("path/to/data.xlsx")
data <- loader$load()
summary <- loader$summarize(data)
```

### Report Pattern

Reports load pre-computed results, not raw data:

```r
results <- readRDS("../data/processed/analysis_results.rds")
data <- results$data
# Use results for tables and plots
```

## Data Files

- Raw data: `deadlift-study/TESE - Análise estatística (Preliminar).xlsx`
- Processed results: `data/processed/*.rds`
- Reports output: `docs/*.html` (GitHub Pages compatible)

## GitHub Pages

Reports render to `docs/` folder. Images in `docs/*_files/` must NOT use Git LFS (configured in `.gitattributes`).
