#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* .C calls */
extern void ruben(double *, int *, double *, int *, double *, double *, int *, double *, double *, int *, double *);

static const R_CMethodDef CEntries[] = {
  {"ruben", (DL_FUNC) &ruben, 11},
  {NULL, NULL, 0}
};

void R_init_manta(DllInfo *dll) {
  R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}

void R_unload_manta(DllInfo *info) { // #nocov start
    /* Release resources. */
}                                    // #nocov end
