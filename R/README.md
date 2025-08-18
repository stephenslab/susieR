# SusieR 2.0 Architecture

## Architecture Overview

```
Interface → Constructor → Workhorse → IBSS Core → Backend Methods
```

The package follows a clean separation of concerns with each layer having specific responsibilities.

## Data Objects

The package supports three types of data input:

- **`individual`**: Raw individual-level data with X and y
- **`ss`**: Sufficient statistics data with XtX, Xty, yty, and sample size (n)
- **`rss`**: Summary statistics data with z-scores, correlation matrix (R), and sample size (n)

### Key Files:
- `susie_data.R`: Creates and validates data objects. Handles standardization, intercept adjustment, and data-specific configurations. Automatically assigns appropriate S3 class based on input type.
- `generic_methods.R`: Defines generic S3 methods that operate on data objects (e.g., `single_effect_update()`, `update_variance_components()`, `check_convergence()`).

### Data Flow:
1. User provides raw data, summary statistics, or RSS inputs
2. Constructor validates inputs and creates appropriate data object with S3 class
3. For unmappable effects, individual data is automatically converted to sufficient statistics
4. Generic methods dispatch to appropriate backend based on data class

### Unmappable Effects:
Unmappable effects are specified as a field within data objects:
- `data$unmappable_effects = "none"` - Standard SuSiE (default)
- `data$unmappable_effects = "inf"` - Infinitesimal effects model
- `data$unmappable_effects = "ash"` - Adaptive shrinkage

## Model Components

The SuSiE model represents the response as a sum of L "single effects" plus noise. Each single effect is a sparse vector where at most one element is non-zero.

### Core Algorithm Files:

1. **`single_effect_regression.R`**: Single Effect Regression (SER) implementation
   - `single_effect_regression()`: Fits one sparse effect at a time
   - Computes posterior means, variances, and inclusion probabilities
   - Supports prior variance optimization via `optimize_V` parameter
   - Enhanced `optimize_prior_variance()` handles both standard and unmappable effects

2. **`iterative_bayesian_stepwise_selection.R`**: IBSS algorithm
   - `ibss_initialize()`: Creates initial model state with L effects
   - `ibss_fit()`: Main iteration loop that updates each effect sequentially
   - `ibss_finalize()`: Post-processing to compute credible sets and PIPs

3. **`susie_workhorse.R`**: Main orchestration layer
   - Manages the complete fitting pipeline: initialize → iterate → finalize
   - Coordinates variance component updates
   - Handles convergence checking based on specified method
   - Tracks fit history when `track_fit=TRUE`

## Backend Method Implementations

Each data type has a corresponding backend file implementing the S3 methods defined in `generic_methods.R`:

- `individual_data_methods.R` - Methods for class `individual`
- `sufficient_stats_methods.R` - Methods for class `ss`
- `rss_lambda_methods.R` - Methods for class `rss_lambda`

The backend system allows the same high-level algorithm to work with different data representations through S3 method dispatch. Each backend file contains a 1:1 correspondence with the generic methods, implementing data-specific computations.

Methods check the `data$unmappable_effects` field to apply appropriate algorithms when needed.

## User Interface

`susie_interface.R` provides the main user-facing functions:

- **`susie()`** - Fits SuSiE model with individual-level data (X, y)
- **`susie_ss()`** - Fits SuSiE model with summary statistics (XtX, Xty, yty, n)
- **`susie_rss()`** - Fits SuSiE model with RSS data (z-scores, R matrix)
  - Automatically handles both lambda = 0 (standard RSS) and lambda > 0 (correlated errors)
  - Lambda = 0 routes through sufficient statistics backend
  - Lambda > 0 uses specialized rss_lambda backend

All interface functions support the full range of algorithm options and call the same underlying `susie_workhorse()`.

## S3 Method Dispatch

This package uses R's S3 object-oriented system for method dispatch. When you call a generic function like `single_effect_update(data, model, l)`, R automatically dispatches to the appropriate method based on the class of the `data` object:

- `single_effect_update.individual()` for class `individual`
- `single_effect_update.ss()` for class `ss`
- `single_effect_update.rss_lambda()` for class `rss_lambda`

If a specific method doesn't exist, R falls back to the default method.

## Utility Functions

- `susie_utils.R`: Internal utility functions including:
  - `mom_unmappable()`: Method of moments for unmappable effects
  - `compute_eigen_decomposition()`: Eigen decomposition for unmappable effects
  - `compute_theta_blup()`: BLUP calculations for infinitesimal effects
  - Various helper functions for computations

- `susie_get_functions.R`: Exported accessor functions for extracting results:
  - `susie_get_cs()`: Extract credible sets
  - `susie_get_pip()`: Extract posterior inclusion probabilities
  - Other accessor functions for model components
