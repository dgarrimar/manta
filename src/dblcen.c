#include <R.h>
#include <Rinternals.h>

/* Double Centering for Classical Multidimensional Scaling */

/* NB: this does not duplicate A */
SEXP dblcen(SEXP A)
{
  int n = nrows(A);
  double *a = REAL(A);
  size_t N = n; /* avoid integer overflow with long vectors */

for(int i = 0; i < n; i++) {
  double sum = 0;
  for(int j = 0; j < n; j++) sum += a[i+j*N];
  sum /= n;
  for(int j = 0; j < n; j++) a[i+j*N] -= sum;
}
for(int j = 0; j < n; j++) {
  double sum = 0;
  for(int i = 0; i < n; i++) sum += a[i+j*N];
  sum /= n;
  for(int i = 0; i < n; i++) a[i+j*N] -= sum;
}
return A;
}
