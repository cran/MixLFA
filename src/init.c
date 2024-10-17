#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

// Déclaration de la fonction Rcpp exportée
SEXP _MixLFA_MLFA(SEXP nb_composantes, SEXP nb_facteur, SEXP X, SEXP Y, SEXP Z, SEXP id, SEXP Lam, SEXP beta, SEXP tau, SEXP s_e, SEXP s_b, SEXP max_it, SEXP pro2, SEXP nb_sujet, SEXP J, SEXP p, SEXP fixed_factor);

// Enregistrement des routines
static const R_CallMethodDef CallEntries[] = {
  {"_MixLFA_MLFA", (DL_FUNC) &_MixLFA_MLFA, 17},
  {NULL, NULL, 0}
};

void R_init_MixLFA(DllInfo *dll) {
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
