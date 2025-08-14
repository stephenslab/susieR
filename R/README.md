## Data

The package currently supports three types of data objects, each with its own S3 class:

- **`individual`**: Raw individual-level data with X and y
- **`ss`**: Sufficient statistics data with XtX, Xty, yty, and sample size n
- **`ss_inf`**: Sufficient statistics with infinitesimal effects (inherits from `ss`)

### Key Files:
- `susie_constructors.R`: Creates and validates data objects. Handles standardization, intercept adjustment, and null weight addition. Automatically assigns appropriate S3 class based on input type.
- `susie_generic_methods.R`: Defines generic methods that operate on data objects (e.g., `get_var_y()`, `configure_data()`, `add_eigen_decomposition()`).

### Data Flow:
1. User provides raw data or summary statistics
2. Constructor validates inputs and creates appropriate data object with S3 class
3. For non-sparse methods, individual data is automatically converted to sufficient statistics
4. Generic methods dispatch to appropriate backend based on data class

## Model

The SuSiE model represents the response as a sum of L "single effects" plus noise. Each single effect is a sparse vector where at most one element is non-zero (though the package also supports non-sparse extensions for infinitesimal effects).

### Core Algorithm Components:

1. **`ser.R`**: Single Effect Regression (SER) implementation
   - `single_effect_regression()`: Fits one sparse effect at a time
   - Computes posterior means, variances, and inclusion probabilities
   - Supports prior variance optimization via `optimize_V` parameter

2. **`ibss.R`**: Iterative Bayesian Stepwise Selection (IBSS) algorithm
   - `ibss_initialize()`: Creates initial model state with L effects, handles `model_init` for warm starts
   - `ibss_fit()`: Main iteration loop that updates each effect sequentially via `single_effect_update()`
   - `ibss_finalize()`: Post-processing to compute credible sets, PIPs, and format output

3. **`susie_workhorse.R`**: Main orchestration layer
   - Manages the complete fitting pipeline: initialize → iterate → finalize
   - Handles convergence checking (ELBO-based for standard, PIP-based for non-sparse)
   - Updates variance components between iterations
   - Tracks fit history when `track_fit=TRUE`

## Backend Implementations

Each data type has a corresponding backend file that implements the S3 methods defined in `susie_generic_methods.R`:

- `susie_individual_backend.R` - Methods for class `individual`
- `susie_ss_backend.R` - Methods for class `ss`  
- `susie_non_sparse_backend.R` - Methods for class `ss_inf` (and future `ss_ash`)

The backend system allows the same high-level algorithm to work with different data representations through S3 method dispatch. Non-sparse methods inherit most functionality from the `ss` backend, only overriding methods where the algorithm differs (e.g., variance estimation, convergence checking).

## Interface

`susie_interface.R` currently provides the two main user-facing functions:

- **`susie()`** - Fits SuSiE model with individual-level data (X, y)
- **`susie_ss()`** - Fits SuSiE model with summary statistics (XtX, Xty, yty, n)

Both functions support all algorithm options and call the same underlying `susie_workhorse()`.

## A note on S3 classes

This package uses R's S3 object-oriented system for method dispatch. When you call a generic function like `single_effect_update(data, ...)`, R automatically dispatches to the appropriate method based on the class of the `data` object:

- `single_effect_update.individual()` for class `individual`
- `single_effect_update.ss()` for class `ss`
- `single_effect_update.ss_inf()` for class `ss_inf`

### Class Inheritance

The non-sparse classes use multiple inheritance:
```r
class(data) <- c("ss_inf", "ss")
```

This means that:
1. R first looks for `method.ss_inf()`
2. If not found, it falls back to `method.ss()`
3. If neither exists, it calls `method.default()`

This inheritance allows non-sparse methods to reuse most of the `ss` implementation, only overriding specific methods where the algorithm differs (e.g., variance updates, convergence checking).
