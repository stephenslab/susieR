#include <R.h>
#include <Rinternals.h>
#include <algorithm>
#include <cmath>
#include <climits>
#include <vector>

// The derivative numerator is a cubic. Its derivative is quadratic, so
// splitting at those turning points brackets EVERY stationary root in [-1,1].
// Bisection avoids unstable divisions in the closed-form cubic formula.
static long double poly(const long double *a, long double x) {
  return ((a[3]*x+a[2])*x+a[1])*x+a[0];
}
static std::vector<long double> candidates(long double xx, long double xh,
    long double hh, long double xy, long double hy) {
  long double A=1+xx, B=xh, D=hh, C=xy, E=hy;
  long double a[4] = {A*C*E-B*(C*C+A), A*E*E-D*(C*C+A)-2*B*B,
                      B*E*E-C*D*E-3*B*D, -D*D};
  std::vector<long double> result{0,-1,1}, breaks{-1,1};
  long double scale=0;
  for (int k=0; k<4; ++k) scale=std::max(scale,std::abs(a[k]));
  if (scale==0) return result;
  for (int k=0; k<4; ++k) a[k]/=scale;
  long double qa=3*a[3], qb=2*a[2], qc=a[1];
  if (qa==0) {
    if (qb!=0) { long double t=-qc/qb; if(t>-1 && t<1) breaks.push_back(t); }
  } else {
    long double disc=qb*qb-4*qa*qc;
    if (disc>=0) {
      long double q=-0.5L*(qb+std::copysign(std::sqrt(disc),qb));
      if(q==0) {
        long double t=-qb/(2*qa); if(t>-1 && t<1) breaks.push_back(t);
      } else {
        long double t=q/qa; if(t>-1 && t<1) breaks.push_back(t);
        t=qc/q; if(t>-1 && t<1) breaks.push_back(t);
      }
    }
  }
  std::sort(breaks.begin(),breaks.end());
  // Include turning points themselves: this catches a repeated root, and
  // adding other interior candidates cannot exceed the actual global maximum.
  for (auto t:breaks) result.push_back(t);
  for (size_t k=1;k<breaks.size();++k) {
    long double lo=breaks[k-1], hi=breaks[k], fl=poly(a,lo), fh=poly(a,hi);
    if (fl==0 || fh==0 || (fl>0)==(fh>0)) continue;
    for (int iteration=0;iteration<60;++iteration) {
      long double mid=(lo+hi)/2, fm=poly(a,mid);
      if(fm==0) { lo=hi=mid; break; }
      if((fl>0)==(fm>0)) {lo=mid; fl=fm;} else hi=mid;
    }
    result.push_back((lo+hi)/2);
  }
  return result;
}

extern "C" SEXP slide_ser(SEXP xx_, SEXP xh_, SEXP hh_, SEXP xy_, SEXP hy_,
                            SEXP V_, SEXP sigma2_, SEXP fixed_) {
  R_xlen_t n=XLENGTH(xx_);
  if(TYPEOF(xx_)!=REALSXP || TYPEOF(xh_)!=REALSXP || TYPEOF(hh_)!=REALSXP ||
     TYPEOF(xy_)!=REALSXP || TYPEOF(hy_)!=REALSXP || TYPEOF(fixed_)!=REALSXP ||
     XLENGTH(xh_)!=n || XLENGTH(hh_)!=n || XLENGTH(xy_)!=n || XLENGTH(hy_)!=n ||
     XLENGTH(fixed_)!=n || n>INT_MAX) Rf_error("Invalid slider sufficient statistics.");
  long double V=Rf_asReal(V_), sigma2=Rf_asReal(sigma2_);
  if(!std::isfinite(V) || V<0 || !std::isfinite(sigma2) || sigma2<=0)
    Rf_error("Invalid slider variance.");
  SEXP ans=PROTECT(Rf_allocMatrix(REALSXP,(int)n,7));
  for(R_xlen_t j=0;j<n;++j) {
    if(j%16384==0) R_CheckUserInterrupt();
    long double xx=REAL(xx_)[j], xh=REAL(xh_)[j], hh=REAL(hh_)[j];
    long double xy=REAL(xy_)[j], hy=REAL(hy_)[j];
    double fixed=REAL(fixed_)[j];
    long double best=R_FINITE(fixed)?fixed:0, best_lbf=-INFINITY;
    std::vector<long double> ds;
    if(R_FINITE(fixed) || V==0) ds.push_back(best);
    else ds=candidates(V*xx/sigma2,V*xh/sigma2,V*hh/sigma2,
                       std::sqrt(V)*xy/sigma2,std::sqrt(V)*hy/sigma2);
    for(auto d:ds) {
      long double s=std::max(0.L,xx+2*d*xh+d*d*hh), t=xy+d*hy;
      long double Q=1+V*s/sigma2;
      long double bf=(V*t*t/(sigma2*sigma2*Q)-std::log1p(V*s/sigma2))/2;
      if(bf>best_lbf) {best=d;best_lbf=bf;}
    }
    long double s=std::max(0.L,xx+2*best*xh+best*best*hh), t=xy+best*hy;
    long double variance=V/(1+V*s/sigma2), mean=variance*t/sigma2;
    // Match SuSiE's no-information convention on zero-design columns.
    if(s==0) {mean=0;variance=0;best_lbf=0;}
    REAL(ans)[j]=best;
    REAL(ans)[j+n]=best_lbf;
    REAL(ans)[j+2*n]=mean;
    REAL(ans)[j+3*n]=variance+mean*mean;
    REAL(ans)[j+4*n]=s;
    REAL(ans)[j+5*n]=t;
    REAL(ans)[j+6*n]=variance;
  }
  UNPROTECT(1);
  return ans;
}
