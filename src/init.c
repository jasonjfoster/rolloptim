#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP _rollport_roll_max_ratio(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _rollport_roll_max_return(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _rollport_roll_min_risk(SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"_rollport_roll_max_ratio",  (DL_FUNC) &_rollport_roll_max_ratio,  6},
    {"_rollport_roll_max_return", (DL_FUNC) &_rollport_roll_max_return, 5},
    {"_rollport_roll_min_risk",   (DL_FUNC) &_rollport_roll_min_risk,   5},
    {NULL, NULL, 0}
};

void R_init_rollport(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}