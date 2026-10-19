/* Truncated multivariate series product: out[Kk[p]] += a[Ii[p]] * b[Jj[p]] over the pair list,
   with the same product on absolute-value arrays (rounding bookkeeping of rseries.py). */
#include <complex.h>
#include <string.h>
void tsmul(long np, const int *Ii, const int *Jj, const int *Kk, long n,
           const double complex *a, const double complex *b, double complex *out,
           const double *aa, const double *ab, double *outa)
{
    memset(out, 0, n * sizeof(double complex));
    if (outa) memset(outa, 0, n * sizeof(double));
    for (long p = 0; p < np; p++) {
        out[Kk[p]] += a[Ii[p]] * b[Jj[p]];
        if (outa) outa[Kk[p]] += aa[Ii[p]] * ab[Jj[p]];
    }
}
