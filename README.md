# Deadlift Study: RIR & Training Outcomes Meta-Analysis

A reproducible meta-analysis examining the relationship between **Repetitions in Reserve (RIR)** (proximity to failure) and resistance training outcomes (strength and hypertrophy).

## Overview

This project implements a multi-level meta-regression analysis to investigate whether training closer to muscular failure produces superior strength and hypertrophy adaptations compared to training with repetitions in reserve.

### Research Questions

1. Does training to failure (RIR = 0) produce greater strength gains than training with repetitions in reserve?
2. Does training to failure produce greater muscle hypertrophy than training with repetitions in reserve?
3. What is the dose-response relationship between RIR and training outcomes?

## Quick Start

```bash
# Prerequisites: R 4.4.2+, Quarto, pre-commit

# 1. Clone and setup
git clone <repository-url>
cd deadlift-study
make setup

# 2. Run analysis pipeline
make pipeline

# 3. Generate reports
make reports
```

## Project Structure

```
deadlift-study/
├── R/                           # Source code (R6 classes, box modules)
│   ├── __init__.R               # Module exports
│   ├── domain/                  # Business entities (DDD layer)
│   │   ├── study.R              # Study entity
│   │   ├── treatment_group.R    # TreatmentGroup entity
│   │   └── effect_size.R        # EffectSize value object
│   ├── calculators/             # Service layer
│   │   └── effect_size_calculator.R  # SMCR/ROMC calculators
│   ├── models/                  # Meta-analysis models
│   ├── validators/              # Data validation (pointblank)
│   └── utils/                   # Shared utilities
│       ├── config.R             # Configuration management
│       ├── logging.R            # Logging utilities
│       └── seed.R               # Reproducibility helpers
├── config/                      # YAML configuration
│   ├── default.yml              # Base configuration
│   └── studies/                 # Study-specific configs
│       ├── strength_rir.yml
│       └── hypertrophy_rir.yml
├── data/                        # Data files (DVC tracked)
│   ├── raw/                     # Original datasets
│   └── processed/               # Analysis-ready data (parquet)
├── analyses/                    # Quarto analysis documents
│   └── index.qmd                # Main analysis notebook
├── tests/                       # Test suite
│   └── testthat/
│       ├── setup.R              # Test fixtures
│       ├── test-treatment_group.R
│       └── test-effect_size_calculator.R
├── docs/                        # Documentation
│   ├── 01_CURRENT_EXPERIMENT_PLAN.md
│   ├── 02_METHODOLOGICAL_IMPROVEMENTS.md
│   ├── 03_PROJECT_ARCHITECTURE_PLAN.md
│   └── 04_RUNBOOKS.md
├── reports/                     # Generated reports
├── pelland/                     # Legacy analysis scripts
├── _targets.R                   # Pipeline definition
├── _quarto.yml                  # Quarto configuration
├── Makefile                     # Automation commands
├── .pre-commit-config.yaml      # Pre-commit hooks
└── .Renviron                    # Environment variables
```

## Architecture

### R6 Class Hierarchy

The project follows **Domain-Driven Design** with R6 classes:

```
Domain Layer (R/domain/)
├── Study               # Research study entity
│   ├── study_id, first_author, publication_year
│   ├── add_treatment_group()  → chainable
│   └── validate() → {is_valid, errors}
├── TreatmentGroup      # Experimental group within study
│   ├── group_id, rir, mean_pre, mean_post, sd_pre, n
│   ├── calculate_mean_change()
│   └── validate()
└── EffectSize          # Calculated effect size (value object)
    ├── effect_estimate, sampling_variance
    └── to_list()

Service Layer (R/calculators/)
├── SMCRCalculator      # Standardized Mean Change (Raw)
│   ├── calculate(TreatmentGroup) → EffectSize
│   └── calculate_batch(list) → list<EffectSize>
└── ROMCCalculator      # Response Ratio (log scale)
    ├── calculate(TreatmentGroup) → EffectSize
    └── calculate_batch(list) → list<EffectSize>
```

### Module System (box)

Code is organized using [box](https://klmr.me/box/) for Python-like explicit imports:

```r
# In _targets.R
box::use(
  R/utils/config[load_config, get_param],
  R/domain/study[Study],
  R/calculators/effect_size_calculator[SMCRCalculator]
)
```

### Pipeline Stages (targets)

```
┌─────────────────┐     ┌─────────────────┐     ┌─────────────────┐
│  Stage 1:       │────▶│  Stage 2:       │────▶│  Stage 3:       │
│  Configuration  │     │  Data Loading   │     │  Transformation │
│                 │     │                 │     │                 │
│  config_file    │     │  raw_data_file  │     │  treatment_grps │
│  analysis_cfg   │     │  raw_study_data │     │  effect_sizes   │
│  random_seed    │     │                 │     │  meta_data      │
└─────────────────┘     └─────────────────┘     └─────────────────┘
                                                        │
                        ┌─────────────────┐             │
                        │  Stage 5:       │◀────────────┘
                        │  Meta-Analysis  │
                        │                 │
                        │  (Models TBD)   │
                        └─────────────────┘
                                │
                        ┌─────────────────┐
                        │  Stage 6:       │
                        │  Outputs        │
                        │                 │
                        │  processed_data │
                        │  pipeline_done  │
                        └─────────────────┘
```

## Documentation

| Document | Description |
|----------|-------------|
| [01_CURRENT_EXPERIMENT_PLAN.md](docs/01_CURRENT_EXPERIMENT_PLAN.md) | Original analysis methodology from Pelland script |
| [02_METHODOLOGICAL_IMPROVEMENTS.md](docs/02_METHODOLOGICAL_IMPROVEMENTS.md) | Proposed enhancements: Bayesian, causal, conformal |
| [03_PROJECT_ARCHITECTURE_PLAN.md](docs/03_PROJECT_ARCHITECTURE_PLAN.md) | Technical architecture and design patterns |
| [04_RUNBOOKS.md](docs/04_RUNBOOKS.md) | SOPs, troubleshooting, and workflows |

## Technology Stack

| Category | Tool | Purpose |
|----------|------|---------|
| **Language** | R 4.4.2 | Statistical computing |
| **OOP** | R6 | Object-oriented design with reference semantics |
| **Modules** | box | Python-like explicit imports |
| **Pipeline** | targets | Workflow orchestration with caching |
| **Packages** | renv | Dependency management and lockfile |
| **Validation** | pointblank | Data quality checks |
| **Reports** | Quarto | Reproducible documents |
| **Formatting** | Air | Fast Rust-based code formatting |
| **Linting** | lintr/Jarl | Code quality checks |
| **Testing** | testthat + hedgehog | Unit + property-based tests |
| **Config** | config | YAML-based environment configs |
| **Logging** | logger | Production-grade logging |

## Configuration

### Environment Variables (.Renviron)

```bash
R_CONFIG_ACTIVE=default      # Options: default, debug, quick, full
DEADLIFT_SEED=42             # Reproducibility seed
DEADLIFT_LOG_LEVEL=INFO      # Logging: DEBUG, INFO, WARN, ERROR
```

### Configuration Profiles (config/default.yml)

| Profile | Use Case |
|---------|----------|
| `default` | Standard analysis settings |
| `debug` | Verbose logging, minimal iterations |
| `quick` | Fast development iterations |
| `full` | Production: 8 chains, 10k iterations |

## Key Commands

```bash
# Setup
make setup              # First-time setup (renv + hooks)
make install-r          # Install R packages only
make install-hooks      # Install pre-commit hooks

# Pipeline
make pipeline           # Run analysis pipeline
make pipeline-parallel  # Run with 4 parallel workers
make pipeline-viz       # Visualize pipeline DAG
make pipeline-status    # Check target status
make pipeline-clean     # Invalidate all targets

# Quality
make check              # Format check + lint + test
make format             # Format code (Air)
make lint               # Lint code (lintr/Jarl)
make test               # Run all tests
make test-coverage      # Run with coverage report
make check-coverage     # Verify 90%+ coverage

# Reports
make reports            # Render all Quarto reports

# Data
make validate-data      # Run pointblank validation
make dvc-pull           # Pull data from DVC remote
make dvc-push           # Push data to DVC remote

# Docker
make docker-build       # Build Docker image
make docker-run         # Run analysis in container

# Help
make help               # Show all commands
```

## Development Principles

This project adheres to:

### SOLID Principles

- **SRP**: Each R6 class has one actor (statistician, data manager, analyst)
- **OCP**: Calculators extensible via strategy pattern
- **LSP**: All TreatmentGroup subtypes substitutable
- **ISP**: Minimal public interfaces (`calculate`, `validate`)
- **DIP**: Services depend on abstractions (EffectSize), not implementations

### CUPID Principles

- **Composable**: Small classes (<10 public methods) that combine
- **Unix**: Each class does one thing well
- **Predictable**: Pure functions, parameterized seeds, immutable configs
- **Idiomatic**: R6 conventions, Pydantic-style validation
- **Domain-based**: Business terms (Study, EffectSize) not implementation (DataFrame, List)

### Scientific Validity

- **Reproducibility**: Parameterized random seeds (never global `set.seed()`)
- **Numerical stability**: Safe division, variance bounds checking
- **Statistical correctness**: Hedge's g correction, REML estimation
- **Data leakage prevention**: Separate fit/transform methods

### Code Standards

- **Coverage**: 90% minimum (enforced via `make check-coverage`)
- **Functions**: <30 lines (Extract Method refactoring)
- **Classes**: <10 public methods (Extract Class refactoring)
- **Testing**: TDD with classicist approach (real objects, no mocks)

## Adding a New Study

1. Add data file to `data/raw/`
2. Create study config in `config/studies/`:

```yaml
# config/studies/my_study.yml
inherits: default
data:
  primary_dataset: "my_study_data.xlsx"
models:
  frequentist:
    random_structure: "~1|study/group"
```

3. Set environment: `R_CONFIG_ACTIVE=my_study`
4. Run: `make pipeline`

## Contributing

See [04_RUNBOOKS.md](docs/04_RUNBOOKS.md) for:
- Development workflow
- Pre-commit hook configuration
- Testing guidelines
- Code review checklist

## License

[Add appropriate license]

## Citation

If you use this analysis, please cite:

```bibtex
@misc{deadlift-study-2024,
  title = {Proximity to Failure and Resistance Training Outcomes: A Meta-Analysis},
  author = {[Authors]},
  year = {2024},
  url = {[Repository URL]}
}
```
