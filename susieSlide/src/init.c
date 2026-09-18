#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
extern SEXP slide_ser(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
static const R_CallMethodDef callMethods[] = {
  {"slide_ser", (DL_FUNC) &slide_ser, 8},
  {NULL, NULL, 0}
};
void R_init_susieSlide(DllInfo *dll) {
  R_registerRoutines(dll, NULL, callMethods, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
