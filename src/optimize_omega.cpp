// -*- mode: C++; c-indent-level: 2; c-basic-offset: 2; indent-tabs-mode: nil; -*-
//
// Rcpp/Armadillo utilities for multi-panel LD mixture.
//
// The RSS-lambda model is:  z ~ N(R b, sigma^2 R + lambda I)
// With K reference panels:  R(omega) = sum_k omega_k R_k
// Omega is optimized on the simplex to maximize the profile Eloglik.
//
// Key insight: the log-likelihood is CONCAVE in omega (from matrix analysis:
// log det is concave on PD matrices, quadratic form is linear in R),
// so the optimization problem is convex with a unique global optimum.
//
// These C++ functions handle the O(p^3) heavy lifting:
//   compute_eloglik_from_eigen — evaluate Eloglik given eigendecomposition
//   compute_null_loglik         — null marginal loglik for initialization
//   eigen_R_omega               — form R(omega) and eigendecompose
//
// The R side handles optimization logic:
//   K=2: Brent's method via optimize()
//   K>2: Frank-Wolfe (conditional gradient) on simplex

#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

// ============================================================================
// Compute the expected log-likelihood for a given eigendecomposition of R.
//
// Eloglik = -p/2 log(2pi) - 1/2 sum log(sigma2*d_j + lambda)
//           - 1/2 * ER2
// where ER2 = z'S^{-1}z - 2*zbar'RS^{-1}z + zbar'RS^{-1}Rzbar
//           - tr(RS^{-1}R Z'Z) + diag(RS^{-1}R)' postb2
// ============================================================================

// [[Rcpp::export]]
double compute_eloglik_from_eigen(
    const arma::vec& eigenvalues,  // p-vector: eigenvalues of R(omega)
    const arma::mat& eigenvectors, // p x p eigenvectors
    const arma::vec& z,            // p-vector of z-scores
    const arma::vec& zbar,         // p-vector: colSums(alpha * mu)
    const arma::vec& diag_postb2,  // p-vector: colSums(alpha * mu2)
    const arma::mat& Z,            // L x p matrix: alpha * mu
    double sigma2,
    double lambda,
    double z_null_norm2) {

  int p = z.n_elem;

  // S_diag = sigma2 * d + lambda;  Dinv = 1 / S_diag
  vec S_diag = sigma2 * eigenvalues + lambda;
  vec Dinv = 1.0 / S_diag;
  for (int j = 0; j < p; j++) {
    if (!std::isfinite(Dinv(j))) Dinv(j) = 0.0;
  }

  vec Vtz = eigenvectors.t() * z;
  vec D = eigenvalues;
  vec DinvD2 = Dinv % (D % D);

  // -1/2 sum log(S_diag)
  double logdet_term = -0.5 * accu(log(S_diag));

  // z'S^{-1}z
  double zSinvz = dot(Dinv % Vtz, Vtz);
  if (lambda > 0) zSinvz += z_null_norm2 / lambda;

  // -2 * zbar' * V * (Dinv .* D .* Vtz)
  vec RSinvz = eigenvectors * (Dinv % D % Vtz);
  double term2 = -2.0 * dot(zbar, RSinvz);

  // (V'zbar)' * DinvD2 * (V'zbar)
  vec Vtzbar = eigenvectors.t() * zbar;
  double term3 = dot(Vtzbar % Vtzbar, DinvD2);

  // -sum_{l} (Z_l * V)^2 * DinvD2
  mat ZV = Z * eigenvectors;
  double term4 = 0.0;
  for (unsigned int l = 0; l < ZV.n_rows; l++) {
    term4 += dot(ZV.row(l) % ZV.row(l), DinvD2.t());
  }
  term4 = -term4;

  // diag(RS^{-1}R)' * postb2
  vec diag_RSinvR = eigenvectors % eigenvectors * DinvD2;
  double term5 = dot(diag_RSinvR, diag_postb2);

  double ER2 = zSinvz + term2 + term3 + term4 + term5;
  return -(double)p / 2.0 * log(2.0 * M_PI) + logdet_term - 0.5 * ER2;
}

// ============================================================================
// Null marginal log-likelihood for panel initialization (no signal model):
// log p(z | R_k, sigma2, lambda) =
//   -1/2 sum log(sigma2*d_k + lambda) - 1/2 [sum(Vtz_k^2 / S_k) + z_null/lambda]
// ============================================================================

// [[Rcpp::export]]
double compute_null_loglik(
    const arma::vec& eigenvalues,
    const arma::vec& Vtz,
    double sigma2,
    double lambda,
    double z_null_norm2) {

  vec S_diag = sigma2 * eigenvalues + lambda;
  double logdet = 0.0;
  double quad = 0.0;
  int p = eigenvalues.n_elem;
  for (int j = 0; j < p; j++) {
    if (S_diag(j) > 0) {
      logdet += log(S_diag(j));
      quad += Vtz(j) * Vtz(j) / S_diag(j);
    }
  }
  if (lambda > 0) quad += z_null_norm2 / lambda;
  return -0.5 * logdet - 0.5 * quad;
}

// ============================================================================
// Form R(omega) = sum_k omega_k R_k and eigendecompose.
// Returns list(values = p-vector descending, vectors = p x p matrix).
// ============================================================================

// [[Rcpp::export]]
Rcpp::List eigen_R_omega(
    const Rcpp::List& panel_R_list,
    const arma::vec& omega,
    int K, int p) {

  mat R_omega = zeros<mat>(p, p);
  for (int k = 0; k < K; k++) {
    R_omega += omega(k) * as<mat>(panel_R_list[k]);
  }
  R_omega = 0.5 * (R_omega + R_omega.t());

  vec eigval;
  mat eigvec;
  eig_sym(eigval, eigvec, R_omega);

  // Sort descending, clip negatives
  vec eigval_desc = reverse(eigval);
  mat eigvec_desc = fliplr(eigvec);
  for (int j = 0; j < p; j++) {
    if (eigval_desc(j) < 0) eigval_desc(j) = 0;
  }

  // Return as plain vector (not Nx1 matrix)
  NumericVector vals = wrap(eigval_desc);
  vals.attr("dim") = R_NilValue;

  return Rcpp::List::create(
    Named("values") = vals,
    Named("vectors") = wrap(eigvec_desc)
  );
}

// ============================================================================
// Combined: form R(omega), eigendecompose, and evaluate Eloglik.
//
// This avoids returning the p x p eigenvector matrix to R and back,
// which is the main bottleneck when R calls eigen_R_omega + compute_eloglik
// separately. Used by the omega optimizer (Brent / Frank-Wolfe).
// ============================================================================

// [[Rcpp::export]]
double eval_omega_eloglik(
    const Rcpp::List& panel_R_list,
    const arma::vec& omega,
    const arma::vec& z,
    const arma::vec& zbar,
    const arma::vec& diag_postb2,
    const arma::mat& Z,
    double sigma2,
    double lambda,
    int K, int p) {

  // Form R(omega) = sum_k omega_k R_k
  mat R_omega = zeros<mat>(p, p);
  for (int k = 0; k < K; k++) {
    R_omega += omega(k) * as<mat>(panel_R_list[k]);
  }
  R_omega = 0.5 * (R_omega + R_omega.t());

  // Eigendecompose
  vec eigval;
  mat V;
  eig_sym(eigval, V, R_omega);

  // Sort descending, clip negatives
  vec D = reverse(eigval);
  V = fliplr(V);
  for (int j = 0; j < p; j++) {
    if (D(j) < 0) D(j) = 0;
  }

  // Compute z_null_norm2 = ||z||^2 - ||V'z||^2
  vec Vtz = V.t() * z;
  double z_null_norm2 = dot(z, z) - dot(Vtz, Vtz);
  if (z_null_norm2 < 0) z_null_norm2 = 0;

  // Now compute Eloglik (same math as compute_eloglik_from_eigen)
  vec S_diag = sigma2 * D + lambda;
  vec Dinv = 1.0 / S_diag;
  for (int j = 0; j < p; j++) {
    if (!std::isfinite(Dinv(j))) Dinv(j) = 0.0;
  }

  vec DinvD2 = Dinv % (D % D);

  double logdet_term = -0.5 * accu(log(S_diag));

  double zSinvz = dot(Dinv % Vtz, Vtz);
  if (lambda > 0) zSinvz += z_null_norm2 / lambda;

  vec RSinvz = V * (Dinv % D % Vtz);
  double term2 = -2.0 * dot(zbar, RSinvz);

  vec Vtzbar = V.t() * zbar;
  double term3 = dot(Vtzbar % Vtzbar, DinvD2);

  mat ZV = Z * V;
  double term4 = 0.0;
  for (unsigned int l = 0; l < ZV.n_rows; l++) {
    term4 += dot(ZV.row(l) % ZV.row(l), DinvD2.t());
  }
  term4 = -term4;

  vec diag_RSinvR = V % V * DinvD2;
  double term5 = dot(diag_RSinvR, diag_postb2);

  double ER2 = zSinvz + term2 + term3 + term4 + term5;
  return -(double)p / 2.0 * log(2.0 * M_PI) + logdet_term - 0.5 * ER2;
}
