# susieR Testing Framework

This directory contains the comprehensive test suite for the susieR 2.0 package, with **>1,000 total tests** ensuring code correctness, stability, and consistency with the reference implementation.

## File Organization

- `testthat/`: Directory containing test files
    - `helper_*.R`: Helper functions used in testing for simulating data, assigning attributes, and more
    - `test_*.R`: Unit tests validating correctness and stability for all functions in susieR 2.0
    - `reference/`: Reference tests ensuring consistency with original susieR 1.0 implementation
