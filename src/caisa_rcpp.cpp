#include "mr_ash.h"

// FUNCTION DEFINITIONS
// --------------------
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

Rcpp::List caisa_rcpp       (const arma::mat& X, const arma::vec& y,
                             const arma::vec& w, const arma::vec& sa2,
                             arma::vec& pi, arma::vec& beta,
                             arma::vec& r, double sigma2, const arma::uvec& o,
                             int maxiter, int miniter,
                             double convtol, double epstol, std::string method_q,
                             bool updatepi, bool updatesigma,
                             bool verbose) {
  
  // See "mr_ash.h"
  
  // ---------------------------------------------------------------------
  // DEFINE SIZES
  // ---------------------------------------------------------------------
  int n                   = X.n_rows;
  int p                   = X.n_cols;
  int K                   = sa2.n_elem;
  
  // ---------------------------------------------------------------------
  // PREDEFINE LOCAL VARIABLES
  // ---------------------------------------------------------------------
  arma::vec varobj(maxiter);
  int iter               = 0;
  int i                  = 0;
  int j;
  
  double a1;
  double a2;
  arma::vec piold;
  arma::vec betaold;
  
  // ---------------------------------------------------------------------
  // PRECALCULATE
  // ---------------------------------------------------------------------
  arma::mat S2inv        = 1 / outerAddition(1/sa2, w);
  S2inv.row(0).fill(epstol);
  
  // ---------------------------------------------------------------------
  // START LOOP : CYCLE THROUGH COORDINATE ASCENT UPDATES
  // ---------------------------------------------------------------------
  for (iter = 0; iter < maxiter; iter++) {
    
    // reset parameters
    a1                   = 0;
    a2                   = 0;
    piold                = pi;
    pi.fill(0);
    betaold              = beta;
    
    // ---------------------------------------------------------------------
    // RUN COORDINATE ASCENT UPDATES : INDEX 1 - INDEX P
    // ---------------------------------------------------------------------
    for (j = 0; j < p; j++){
      
      updatebetaj(X.col(o(i)), w(o(i)), beta(o(i)), r, piold, pi, sigma2, sa2, S2inv.col(o(i)), a1, a2, o(i), p, epstol);
      i++;
      // updatebetaj(X.col(j), w(j), beta(j), r, piold, pi, sigma2, sa2, S2inv.col(j), a1, a2, j, p, epstol);
      
    }
    
    // ---------------------------------------------------------------------
    // CALCULATE VARIATIONAL OBJECTIVE 1
    // ---------------------------------------------------------------------
    varobj(iter)          = arma::dot(r,r) - arma::dot(square(beta), w) + a1;
    
    // ---------------------------------------------------------------------
    // UPDATE SIGMA2 IF REQUESTED
    // ---------------------------------------------------------------------
    if (updatesigma) {
      if (method_q == std::string("sigma_indep_q")) {
        sigma2            = varobj(iter) + p * (1.0 - pi(0)) * sigma2;
        sigma2            = sigma2 / (n + p * (1.0 - pi(0)));
      } else if (method_q == std::string("sigma_dep_q")) {
        sigma2            = varobj(iter) / n;
      }
    }
    
    if (updatepi) {
      // if updatepi == true, we update pi
      piold               = pi;
    }
    // ---------------------------------------------------------------------
    // CALCULATE VARIATIONAL OBJECTIVE 2
    // ---------------------------------------------------------------------
    varobj(iter)          = varobj(iter) / sigma2 / 2.0 +
                            log(2.0 * M_PI * sigma2) / 2.0 * n -
                            dot(pi, log(piold + epstol)) * p + a2;
    
    for (j = 1; j < K; j++){
      varobj(iter)       += pi(j) * log(sa2(j)) * p / 2;
    }
    
    if (!updatepi) {
      // if updatepi == false, we do not update pi
      pi                  = piold;
    }
    
    // ---------------------------------------------------------------------
    // CHECK CONVERGENCE
    // ---------------------------------------------------------------------
    if (iter >= miniter - 1) {
      double beta_norm = arma::norm(beta, 2);
      if (arma::norm(betaold - beta, 2) < convtol * std::max(1.0, beta_norm)) {
        iter++;
        break;
      }
      
      if (iter > 0) {
        if (varobj(iter) > varobj(iter - 1)){
          break;
        }
      }
    }
  }
  
  if (verbose) {
    Rprintf("Mr.ASH terminated at iteration %d: max|beta|=%.4e, sigma2=%.4e, pi0=%.4f\n", 
              iter, arma::max(arma::abs(beta)), sigma2, pi(0));
  }
  // ---------------------------------------------------------------------
  // RETURN VALUES
  // ---------------------------------------------------------------------
  return Rcpp::List::create(Rcpp::Named("beta")    = beta,
                            Rcpp::Named("sigma2")  = sigma2,
                            Rcpp::Named("pi")      = pi,
                            Rcpp::Named("iter")    = iter,
                            Rcpp::Named("varobj")  = varobj.subvec(0,iter-1));
}
