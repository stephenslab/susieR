## Data

- `susie_constructors.R` implements data object construction and validation for different input types.
- `susie_generic_methods.R` defines S3 generic methods for operations on data objects.

## Model

1. `ser.R` implements Single Effect Regression (SER) models using different backends.
2. `ibss.R` implements the Iterative Bayesian Stepwise Selection (IBSS) algorithm core.
    - `ibss_initialize()` sets up initial model state.
    - `ibss_fit()` performs the main fitting iterations.
    - `ibss_finalize()` computes credible sets and final outputs.
3. `susie_engine.R` orchestrates the complete SuSiE model fitting pipeline.

## Backend Implementations

- `susie_individual_backend.R` implements methods for individual-level data (class `individual`).
- `susie_ss_backend.R` implements methods for sufficient statistics (class `ss`).
- `susie_non_sparse_backend.R` implements non-sparse extensions (classes `ss_inf` and TODO: `ss_ash`).

## Interface

- `susie_interface.R` implements the main user-facing functions `susie()` and `susie_ss()`.
- These functions provide a unified interface regardless of data type or method choice.

## A note on S3 classes

TODO
