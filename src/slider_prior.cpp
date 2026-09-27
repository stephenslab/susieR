#include <cpp11.hpp>
#include <R.h>
#include <Rinternals.h>
#include <algorithm>
#include <cmath>
#include <climits>
#include <vector>

// Integrate a fixed finite slider prior after the two residual cross-products.
// Every SNP uses O(K) scalar work, independently of the sample size.
[[cpp11::register]]
SEXP slide_prior_ser_native(SEXP xx_, SEXP xh_, SEXP hh_, SEXP xy_, SEXP hy_,
                           SEXP V_, SEXP sigma2_, SEXP grid_, SEXP prior_,
                           SEXP forced_) {
  R_xlen_t p=XLENGTH(xx_), K=XLENGTH(grid_);
  if (p<1 || p>INT_MAX || K<1 || K>INT_MAX ||
      TYPEOF(xx_)!=REALSXP || TYPEOF(xh_)!=REALSXP ||
      TYPEOF(hh_)!=REALSXP || TYPEOF(xy_)!=REALSXP || TYPEOF(hy_)!=REALSXP ||
      TYPEOF(grid_)!=REALSXP || TYPEOF(prior_)!=REALSXP ||
      TYPEOF(forced_)!=LGLSXP || XLENGTH(xh_)!=p || XLENGTH(hh_)!=p ||
      XLENGTH(xy_)!=p || XLENGTH(hy_)!=p || XLENGTH(forced_)!=p ||
      XLENGTH(prior_)!=K)
    Rf_error("Invalid finite-slider sufficient statistics or prior.");
  long double V=Rf_asReal(V_), sigma2=Rf_asReal(sigma2_), total=0;
  if (!std::isfinite(V) || V<0 || !std::isfinite(sigma2) || sigma2<=0)
    Rf_error("Invalid slider variance.");
  int zero=-1;
  for (R_xlen_t k=0;k<K;++k) {
    if (!R_FINITE(REAL(grid_)[k]) || std::abs(REAL(grid_)[k])>1 ||
        !R_FINITE(REAL(prior_)[k]) || REAL(prior_)[k]<0)
      Rf_error("Invalid finite-slider grid or probabilities.");
    if (REAL(grid_)[k]==0) zero=(int)k;
    total+=REAL(prior_)[k];
  }
  if (total<=0 || zero<0)
    Rf_error("The slider prior needs positive total mass and a grid containing zero.");
  SEXP ans=PROTECT(Rf_allocVector(VECSXP,4));
  SEXP names=PROTECT(Rf_allocVector(STRSXP,4));
  const char* labels[]={"summary","weights","mu","mu2"};
  for (int i=0;i<4;++i) SET_STRING_ELT(names,i,Rf_mkChar(labels[i]));
  Rf_setAttrib(ans,R_NamesSymbol,names);
  SEXP summary=PROTECT(Rf_allocMatrix(REALSXP,(int)p,7));
  SEXP weights=PROTECT(Rf_allocMatrix(REALSXP,(int)p,(int)K));
  SEXP mu=PROTECT(Rf_allocMatrix(REALSXP,(int)p,(int)K));
  SEXP mu2=PROTECT(Rf_allocMatrix(REALSXP,(int)p,(int)K));
  SET_VECTOR_ELT(ans,0,summary); SET_VECTOR_ELT(ans,1,weights);
  SET_VECTOR_ELT(ans,2,mu); SET_VECTOR_ELT(ans,3,mu2);
  std::vector<long double> logmass(K), means(K), seconds(K);
  for (R_xlen_t j=0;j<p;++j) {
    if (j%4096==0) R_CheckUserInterrupt();
    long double xx=REAL(xx_)[j], xh=REAL(xh_)[j], hh=REAL(hh_)[j];
    long double xy=REAL(xy_)[j], hy=REAL(hy_)[j];
    if (!std::isfinite(xx) || !std::isfinite(xh) || !std::isfinite(hh) ||
        !std::isfinite(xy) || !std::isfinite(hy) || xx<0 || hh<0 ||
        LOGICAL(forced_)[j]==NA_LOGICAL)
      Rf_error("Invalid finite-slider sufficient statistics.");
    long double peak=-INFINITY;
    for (R_xlen_t k=0;k<K;++k) {
      long double d=REAL(grid_)[k];
      long double wk=LOGICAL(forced_)[j] ? (k==zero ? 1.L : 0.L) :
        REAL(prior_)[k]/total;
      means[k]=seconds[k]=0;
      logmass[k]=-INFINITY;
      if (wk>0) {
        long double s=std::max(0.L,xx+2*d*xh+d*d*hh);
        long double t=s==0 ? 0.L : xy+d*hy;
        long double q=1+V*s/sigma2;
        long double variance=V/q, mean=variance*t/sigma2;
        long double bf=(V*t*t/(sigma2*sigma2*q)-std::log1p(V*s/sigma2))/2;
        means[k]=mean; seconds[k]=variance+mean*mean;
        logmass[k]=std::log(wk)+bf;
        peak=std::max(peak,logmass[k]);
      }
    }
    long double norm=0;
    for (R_xlen_t k=0;k<K;++k) norm+=std::exp(logmass[k]-peak);
    long double sums[7]={0,peak+std::log(norm),0,0,0,0,0};
    for (R_xlen_t k=0;k<K;++k) {
      long double prob=std::exp(logmass[k]-peak)/norm, d=REAL(grid_)[k];
      REAL(weights)[j+p*k]=prob;
      REAL(mu)[j+p*k]=means[k]; REAL(mu2)[j+p*k]=seconds[k];
      sums[0]+=prob*d;
      sums[2]+=prob*means[k]; sums[3]+=prob*seconds[k];
      sums[4]+=prob*d*means[k]; sums[5]+=prob*d*seconds[k];
      sums[6]+=prob*d*d*seconds[k];
    }
    for (int c=0;c<7;++c) REAL(summary)[j+p*c]=sums[c];
  }
  UNPROTECT(6);
  return ans;
}
