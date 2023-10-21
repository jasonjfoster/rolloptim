#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
 Check these declarations against the C/Fortran source code.
 */

/* .Call calls */
extern SEXP _rollport_roll_max_return(void *, void *, void *, void *, void *);
extern SEXP _rollport_roll_max_utility(void *, void *, void *, void *, void *, void *);
extern SEXP _rollport_roll_min_var(void *, void *, void *, void *);

static const R_CallMethodDef CallEntries[] = {
  {"_rollport_roll_max_return",  (DL_FUNC) &_rollport_roll_max_return,  5},
  {"_rollport_roll_max_utility", (DL_FUNC) &_rollport_roll_max_utility, 6},
  {"_rollport_roll_min_var",     (DL_FUNC) &_rollport_roll_min_var,     4},
  {NULL, NULL, 0}
};

void R_init_rollport(DllInfo *dll)
{
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}