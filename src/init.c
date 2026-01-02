#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP _rolloptim_roll_max_mean(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _rolloptim_roll_max_utility(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _rolloptim_roll_min_rss(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _rolloptim_roll_min_var(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"_rolloptim_roll_max_mean",    (DL_FUNC) &_rolloptim_roll_max_mean,    6},
    {"_rolloptim_roll_max_utility", (DL_FUNC) &_rolloptim_roll_max_utility, 6},
    {"_rolloptim_roll_min_rss",     (DL_FUNC) &_rolloptim_roll_min_rss,     5},
    {"_rolloptim_roll_min_var",     (DL_FUNC) &_rolloptim_roll_min_var,     6},
    {NULL, NULL, 0}
};

void R_init_rolloptim(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}