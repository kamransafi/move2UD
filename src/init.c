#include <R.h>
#include <Rinternals.h>
#include <stdlib.h>
#include <R_ext/Rdynload.h>

extern SEXP bgb(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP dbbmm2(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP llBGBvar(SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"bgb",      (DL_FUNC) &bgb,      10},
    {"dbbmm2",   (DL_FUNC) &dbbmm2,   10},
    {"llBGBvar", (DL_FUNC) &llBGBvar,  2},
    {NULL, NULL, 0}
};

void R_init_move2UD(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
