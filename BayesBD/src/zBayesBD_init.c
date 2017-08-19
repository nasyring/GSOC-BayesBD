#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* .Call calls */
extern SEXP BayesBD_BayesBDbinary(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP BayesBD_BayesBDnormal(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP BayesBD_eigenfun(SEXP, SEXP);
extern SEXP BayesBD_unisliceL(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"BayesBD_BayesBDbinary", (DL_FUNC) &BayesBD_BayesBDbinary, 8},
    {"BayesBD_BayesBDnormal", (DL_FUNC) &BayesBD_BayesBDnormal, 9},
    {"BayesBD_eigenfun",      (DL_FUNC) &BayesBD_eigenfun,      2},
    {"BayesBD_unisliceL",     (DL_FUNC) &BayesBD_unisliceL,     8},
    {NULL, NULL, 0}
};

void R_init_BayesBD(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}