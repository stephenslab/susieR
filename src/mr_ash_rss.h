#ifndef MR_ASH_RSS_H
#define MR_ASH_RSS_H

#include <RcppArmadillo.h>
#include <cmath>
#include <vector>
#include <string>
#include <algorithm>
#include <unordered_map>
#include <omp.h>

using namespace arma;
using namespace std;

/**
 * Softmax function
 *
 * @param x Input vector
 * @return Softmax output vector
 */
vec softmax_rss(const vec& x) {
	vec y = exp(x - max(x));
	return y / sum(y);
}

/**
 * Bayesian regression with Normal prior from sufficient statistics
 *
 * @param xTx X'X (scalar)
 * @param xTy X'y (scalar)
 * @param sigma2_e Error variance
 * @param sigma2_0 Prior variance
 * @return An unordered_map containing the least-squares estimate (bhat, s2), the posterior mean and standard deviation (mu1, sigma2_1), and the log-Bayes factor (logbf)
 */
unordered_map<string, double> bayes_ridge_sufficient(double xTx, double xTy, double sigma2_e, double sigma2_0) {
	// Compute the least-squares estimate and its variance
	double bhat = xTy / xTx;
	double s2 = sigma2_e / xTx;

	// Compute the posterior mean and variance assuming a normal prior with zero mean and variance sigma2_0
	double sigma2_1 = 1 / (1 / s2 + 1 / sigma2_0);
	double mu1 = sigma2_1 / s2 * bhat;

	// Compute the log-Bayes factor
	double logbf = log(s2 / (sigma2_0 + s2)) / 2 + (pow(bhat, 2) / s2 - pow(bhat, 2) / (sigma2_0 + s2)) / 2;

	// Return the least-squares estimate (bhat, s2), the posterior mean and standard deviation (mu1, sigma2_1), and the log-Bayes factor (logbf)
	return {{"bhat", bhat}, {"s2", s2}, {"mu1", mu1}, {"sigma2_1", sigma2_1}, {"logbf", logbf}};
}

/**
 * Bayesian regression with mixture-of-normals prior from sufficient statistics
 *
 * @param xTx X'X (scalar)
 * @param xTy X'y (scalar)
 * @param sigma2_e Error variance
 * @param w0 Mixture weights
 * @param sigma2_0 Mixture variances
 * @return An unordered_map containing the log-Bayes factor (logbf), the posterior assignment probabilities (w1), the posterior mean (mu1) and variance (sigma2_1) of the coefficients, and the posterior mean (mu1_k) and variance (sigma2_1_k) for each mixture component
 */
unordered_map<string, vec> bayes_mix_sufficient(double xTx, double xTy, double sigma2_e, const vec& w0, const vec& sigma2_0) {
	// Get the number of mixture components (K)
	int K = sigma2_0.n_elem;

	// Compute the Bayes factors and posterior statistics separately for each mixture component
	mat out(K, 5);
	for (int i = 0; i < K; i++) {
		unordered_map<string, double> ridge_out = bayes_ridge_sufficient(xTx, xTy, sigma2_e, sigma2_0[i]);
		out(i, 0) = ridge_out["bhat"];
		out(i, 1) = ridge_out["s2"];
		out(i, 2) = ridge_out["mu1"];
		out(i, 3) = ridge_out["sigma2_1"];
		out(i, 4) = ridge_out["logbf"];
	}

	// Compute the posterior assignment probabilities for the latent indicator variable
	vec w1 = softmax_rss(out.col(4) + log(w0));

	// Compute the posterior mean (mu1) and variance (sigma2_1) of the regression coefficients
	vec mu1_k_vec = out.col(2);
	vec sigma2_1_k_vec = out.col(3);

	double mu1 = sum(w1 % mu1_k_vec);
	double sigma2_1 = sum(w1 % (square(mu1_k_vec) + sigma2_1_k_vec)) - pow(mu1, 2);

	// Compute the log-Bayes factor as a linear combination of the individual BFs for each mixture component
	double u = max(out.col(4));
	double logbf = u + log(sum(w0 % exp(out.col(4) - u)));

	// Return the posterior assignment probabilities (w1), the posterior mean (mu1) and variance (sigma2_1) of the coefficients,
	// and the posterior mean (mu1_k) and variance (sigma2_1_k) for each mixture component, and the log-Bayes factor (logbf)
	return {{"w1", w1}, {"mu1", vec(1, fill::value(mu1))}, {"sigma2_1", vec(1, fill::value(sigma2_1))},
		{"mu1_k", mu1_k_vec}, {"sigma2_1_k", sigma2_1_k_vec}, {"logbf", vec(1, fill::value(logbf))}};
}

/**
 * Bayesian multiple regression with mixture-of-normals prior from sufficient statistics
 *
 * @param XTy X'y vector
 * @param XTX X'X matrix
 * @param yTy y'y scalar
 * @param n Sample size
 * @param sigma2_e Error variance
 * @param sigma2_0 Mixture variances
 * @param w0 Mixture weights
 * @param mu1_init Initial value for mu1
 * @param tol Convergence tolerance
 * @param max_iter Maximum number of iterations
 * @param update_w0 Whether to update w0
 * @param update_sigma Whether to update sigma2_e
 * @param compute_ELBO Whether to compute the Evidence Lower Bound (ELBO)
 * @param ncpus Number of CPUs to use for parallel processing
 * @return An unordered_map containing the posterior assignment probabilities (w1), the posterior mean (mu1) and variance (sigma2_1) of the coefficients, the error variance (sigma2_e), the mixture weights (w0), and optionally the ELBO
 */
unordered_map<string, mat> mr_ash_sufficient(const vec& XTy, const mat& XTX, double yTy, int n, double& sigma2_e,
                                             const vec& sigma2_0, vec& w0, const vec& mu1_init, double tol = 1e-8,
                                             int max_iter = 1e5, bool update_w0 = true, bool update_sigma = true,
                                             bool compute_ELBO = true, bool verbose = false, int ncpus = 1) {
	// Set the number of threads for OpenMP
	int nProcessors = omp_get_max_threads();
	if (ncpus < nProcessors) nProcessors = ncpus;
	omp_set_num_threads(nProcessors);

	// Initialize parameters
	int p = XTX.n_cols;
	int K = sigma2_0.n_elem;
	vec mu1_t = mu1_init;
	vec sigma2_1_t(p, fill::zeros);
	mat w1_t(p, K, fill::zeros);
	mat mu1_k_t(p, K, fill::zeros);
	mat sigma2_1_k_t(p, K, fill::zeros);
	vec err(p, fill::value(datum::inf));
	int t = 0;
	double ELBO = 0;

	// Iterate until convergence
	while (any(err > tol)) {
		double ELBO0 = ELBO;
		double var_part_ERSS = 0;
		double neg_KL = 0;

		// Update iterator
		t++;

		// Exit loop if maximum number of iterations is reached
		if (t > max_iter) {
			cerr << "Max number of iterations reached. Try increasing max_iter." << endl;
			break;
		}

		// Save current estimates
		vec mu1_tminus1 = mu1_t;

		vec XTrbar = XTy - XTX * mu1_t;

		// Loop through the variables
		for (int j = 0; j < p; j++) {
			// Remove j-th effect from expected residuals
			vec XTrbar_j = XTrbar + XTX.col(j) * mu1_t[j];

			double xTrbar_j = XTrbar_j[j];
			double xTx = XTX(j, j);

			// Run Bayesian SLR
			unordered_map<string, vec> bfit = bayes_mix_sufficient(xTx, xTrbar_j, sigma2_e, w0, sigma2_0);

			// Update variational parameters
			mu1_t[j] = bfit["mu1"][0];
			sigma2_1_t[j] = bfit["sigma2_1"][0];
			w1_t.row(j) = bfit["w1"].t();
			mu1_k_t.row(j) = bfit["mu1_k"].t();
			sigma2_1_k_t.row(j) = bfit["sigma2_1_k"].t();

			// Compute ELBO parameters
			if (compute_ELBO) {
				var_part_ERSS += sigma2_1_t[j] * xTx;
				neg_KL += bfit["logbf"][0] + (1 / (2 * sigma2_e)) * (-2 * xTrbar_j * mu1_t[j] + (xTx * (sigma2_1_t[j] + pow(mu1_t[j], 2))));
			}

			// Update expected residuals
			XTrbar = XTrbar_j - XTX.col(j) * mu1_t[j];
		}

		// Update w0 if requested
		if (update_w0) {
			w0 = sum(w1_t, 0).t() / p;
		}

		// Compute distance in mu1 between two successive iterations
		err = abs(mu1_t - mu1_tminus1);

		// Compute ERSS and ELBO
		double ERSS = yTy - 2 * dot(XTy, mu1_t) + as_scalar(mu1_t.t() * XTX * mu1_t) + var_part_ERSS;
		if (compute_ELBO) {
			ELBO = -0.5 * log(n) - 0.5 * n * log(2 * datum::pi * sigma2_e) - (1 / (2 * sigma2_e)) * ERSS + neg_KL;
			if (verbose) {
				// Print out useful info
				cout << "Iteration: " << t << ", Max beta diff: " << max(err) << ", ELBO diff: " << ELBO - ELBO0 << ", ELBO: " << ELBO << endl;
			}
		} else {
			if (verbose) {
				// Print out useful info
				cout << "Iteration: " << t << ", Max beta diff: " << max(err) << endl;
			}
		}

		// Update residual variance if requested
		if (update_sigma) {
			sigma2_e = ERSS / n;
		}
	}

	// Return the posterior assignment probabilities (w1), the posterior mean (mu1) and variance (sigma2_1) of the coefficients,
	// the error variance (sigma2_e), the mixture weights (w0), and optionally the ELBO
	return {{"mu1", mat(mu1_t)}, {"sigma2_1", mat(sigma2_1_t)}, {"w1", w1_t},
		{"sigma2_e", mat(1, 1, fill::value(sigma2_e))}, {"w0", mat(w0)}, {"ELBO", mat(1, 1, fill::value(ELBO))}};
}

/**
 * Rescale posterior mean and covariance
 *
 * @param mu1 Posterior mean vector
 * @param sigma2_1 Posterior covariance matrix
 * @param sx Scaling vector
 * @return An unordered_map containing the rescaled posterior mean (mu1_orig) and covariance (sigma2_1_orig)
 */
unordered_map<string, mat> rescale_post_mean_covar(const vec& mu1, const mat& sigma2_1, const vec& sx) {
	vec mu1_orig = mu1 / sx;
	mat sigma2_1_orig = diagmat(1 / sx) * sigma2_1 * diagmat(1 / sx);
	return {{"mu1_orig", mat(mu1_orig)}, {"sigma2_1_orig", sigma2_1_orig}};
}

/**
 * Bayesian multiple regression with mixture-of-normals prior
 *
 * @param bhat Observed effect sizes (standardized)
 * @param shat Standard errors of effect sizes
 * @param z Z-scores
 * @param R Correlation matrix
 * @param var_y Variance of the outcome
 * @param n Sample size
 * @param sigma2_e Error variance
 * @param s0 Prior variances for the mixture components
 * @param w0 Prior weights for the mixture components
 * @param mu1_init Initial value for the posterior mean of the coefficients
 * @param tol Convergence tolerance
 * @param max_iter Maximum number of iterations
 * @param update_w0 Whether to update the mixture weights
 * @param update_sigma Whether to update the error variance
 * @param compute_ELBO Whether to compute the Evidence Lower Bound (ELBO)
 * @param standardize Whether to standardize the input data
 * @param ncpus Number of CPUs to use for parallel processing
 * @return An unordered_map containing the posterior mean (mu1) and covariance (sigma2_1) of the coefficients, the posterior assignment probabilities (w1), the error variance (sigma2_e), the mixture weights (w0), and optionally the ELBO
 */
unordered_map<string, mat> mr_ash_rss(const vec& bhat, const vec& shat, const vec& z, const mat& R, double var_y, int n,
                                      double sigma2_e, const vec& s0, vec& w0, const vec& mu1_init, double tol = 1e-8,
                                      int max_iter = 1e5, bool update_w0 = true, bool update_sigma = true, bool compute_ELBO = true,
                                      bool standardize = false, int ncpus = 1) {
	// Get number of variables
	int p = z.n_elem;

	// Initialize regression coefficients to 0 if not provided
	vec mu1_init_use = mu1_init;
	if (mu1_init.is_empty()) {
		mu1_init_use = vec(p, fill::zeros);
	}

	// Compute Z-scores if not provided
	vec z_use = z;
	if (z.is_empty()) {
		z_use = bhat / shat;
	}

	// Compute PVE-adjusted Z-scores if sample size is provided
	vec adj(p, fill::ones);
	if (std::isfinite(n)) {
		adj = (n - 1) / (square(z_use) + n - 2);
		z_use %= sqrt(adj);
	}
	// Compute X'X and X'y
	mat XtX;
	vec Xty;
	if (std::isfinite(var_y) && !shat.is_empty()) {
		vec XtXdiag = var_y * adj / square(shat);
		XtX = diagmat(sqrt(XtXdiag)) * R * diagmat(sqrt(XtXdiag));
		XtX = 0.5 * (XtX + XtX.t());
		Xty = z_use % sqrt(adj) % (var_y / shat);
	} else {
		// The effects are on the standardized X, y scale
		XtX = (n - 1) * R;
		Xty = z_use * sqrt(n - 1);
		var_y = 1.0;
	}

	// Adjust X'X and X'y if X is standardized
	vec sx(p, fill::ones);
	if (standardize) {
		vec dXtX = XtX.diag();
		sx = sqrt(dXtX / (n - 1));
		sx.replace(0, 1);
		XtX = diagmat(1 / sx) * XtX * diagmat(1 / sx);
		Xty /= sx;
		mu1_init_use %= sx;
	}

	// Run variational inference
	unordered_map<string, mat> result = mr_ash_sufficient(Xty, XtX, var_y * (n - 1), n, sigma2_e, s0, w0, mu1_init_use,
	                                                      tol, max_iter, update_w0, update_sigma, compute_ELBO, standardize, ncpus);

	// Rescale posterior mean and covariance if X was standardized
	if (standardize) {
		unordered_map<string, mat> out_adj = rescale_post_mean_covar(vectorise(result["mu1"]), result["sigma2_1"], sx);
		result["mu1"] = out_adj["mu1_orig"];
		result["sigma2_1"] = out_adj["sigma2_1_orig"];
	}

	return {{"mu1", result["mu1"]}, {"sigma2_1", result["sigma2_1"]}, {"w1", result["w1"]},
		{"sigma2_e", result["sigma2_e"]}, {"w0", result["w0"]}, {"ELBO", result["ELBO"]}};
};

#endif
