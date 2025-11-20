# susieR 2.0 Architecture

## Overview

susieR 2.0 implements a unified architecture incorporating various extensions to the Sum of Single Effects model for Bayesian variable selection regression. 
The package supports multiple data types (individual-level, sufficient statistics, regression summary statistics) through a single algorithmic pipeline using S3 method dispatch.

## Architecture Diagram

```
Interface → Constructor → Workhorse → IBSS Core → Backend Methods
               ↓            ↓
          (data, params) → model
```

## Core Object Definitions

The architecture revolves around three key objects:

### **Data Object**
- **Purpose**: Contains input data in processed, algorithm-ready form
- **S3 Classes**: `individual`, `ss`, `rss_lambda` (determines method dispatch)
- **Mutability**: Immutable - never modified after creation
- **Contents**: 
  - Input matrices: X/y (individual), XtX/Xty/yty (ss), z/R (rss_lambda)
  - Metadata: n, p
  - Scaling attributes: For compute_Xb() compatibility
  - Specialized fields: Eigen decomposition for unmappable effects/rss_lambda


### **Params Object**
- **Purpose**: Contains ALL algorithm parameters and user settings
- **Mutability**: Immutable - never modified after validation
- **Contents**:
  - Algorithm parameters: L, max_iter, tol, convergence_method
  - Estimation settings: estimate_prior_method, estimate_residual_method
  - Model options: unmappable_effects, refine, standardize, intercept

### **Model Object**
- **Purpose**: Contains fitted SuSiE model state, results, and algorithm outputs
- **Mutability**: Mutable - updated throughout fitting process
- **Contents**:
  - Model matrices: alpha, mu, mu2, V, sigma2
  - Fitted values: Xr (individual), XtXr (ss), Rz (rss_lambda)
  - Algorithm outputs: ELBO, niter, converged
  - Final results: credible sets, PIPs, intercept, z-scores

## Constructor Pattern

### **Constructor Workflow**:
1. **Interface functions** (`susie()`, `susie_ss()`, and `susie_rss()`) take user inputs and call constructors functions
2. **Constructors** create validated (data, params) objects
3. **Workhorse** Validated (data, params) objects are directly forwarded to the workhorse function for the main SuSiE algorithm

### **Constructor Functions** (`susie_constructors.R`):
- `individual_data_constructor()` → Processes X, y matrices → (data, params)
- `sufficient_stats_constructor()`→ Processes XtX, Xty, yty → (data, params)  
- `summary_stats_constructor()`: Routes RSS inputs based on lambda parameter
   - If `lambda = 0` → Converts RSS data to SS → `sufficient_stats_constructor()` → (data, params)
   - If `lambda > 0` → `rss_lambda_constructor()`→ Processes z, R for regularized LD → (data, params)

### **Data Type Support**:

Each data object receives an S3 class to automatically route to the appropriate backend function based on the data object's S3 class.

- **`individual`**: Individual-level data (X, y matrices)
- **`ss`**: Sufficient statistics (XtX, Xty, yty, n)
- **`rss_lambda`**: RSS with regularized LD matrix (z, R, lambda > 0)

## Model Components

### Core Algorithm Files:

1. **`susie_workhorse.R`**: Main orchestration layer
   - Manages the complete fitting pipeline: initialize → iterate → finalize
   - Coordinates variance component updates
   - Handles convergence checking based on specified method
   - Tracks fit history when `track_fit=TRUE`

2. **`iterative_bayesian_stepwise_selection.R`**: IBSS algorithm
   - `ibss_initialize()`: Creates initial model state with L effects
   - `ibss_fit()`: Main iteration loop that updates each effect sequentially
   - `ibss_finalize()`: Post-processing to compute credible sets and PIPs

3. **`single_effect_regression.R`**: Single Effect Regression (SER) implementation
   - `single_effect_regression()`: Fits one sparse effect at a time
   - `optimize_prior_variance()`: Optimizes the prior variance for the lth effect
   - `single_effect_update()`: Implements the complete SER update pipeline

## Backend Method Implementations

Each data type has a corresponding backend file implementing the S3 methods defined in `generic_methods.R`:

- `individual_data_methods.R` - Methods for class `individual`
- `sufficient_stats_methods.R` - Methods for class `ss`
- `rss_lambda_methods.R` - Methods for class `rss_lambda`

The backend system allows the same high-level algorithm to work with different data representations through S3 method dispatch. Each backend file contains a 1:1 correspondence with the generic methods, implementing data-specific computations.

## Utility Functions

- **`susie_utils.R`**: Internal utility functions organized by:
  - **Fundamental Building Blocks** (general purpose-helpers, matrix operations)
  - **Data Processing & Validation** (input validation, data conversion)
  - **Model Initialization** (set up model state)
  - **Core Algorithm Components** (posterior mean calculation, lbf calculation)
  - **Variance Esimation** (residual variance and unmappable effects variance esimation)
  - **Convergence & Optimization** (Convergence checking, ELBO computation)
  - **Credible Sets & Post-processing** (Generate credible sets, pips)

- **`susie_rss_utils.R`**: Internal utility functions specific to RSS data with lambda > 0, organized by:
  - **Fundamental Computations** (core RSS computations)
  - **RSS Model Methods** (lambda estimation, precomputations)
  - **Diagnostic & Quality Control** (detect allele switch)

- **`susie_get_functions.R`**: Exported accessor functions for extracting results:
  - `susie_get_cs()`: Extract credible sets
  - `susie_get_pip()`: Extract posterior inclusion probabilities
  - Other accessor functions for model components
