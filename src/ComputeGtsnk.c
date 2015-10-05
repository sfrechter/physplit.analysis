#include <R.h>

void ComputeGtsnk(double *X, double *a, int *dim, double *G){
  int T = dim[0], S = dim[1], N = dim[2], K = dim[3];
  int t,s,n,k, tp;
  int iXs, iXn = 0;
  double *g, *a0, ar;
  g = (double *)R_alloc(T, sizeof(double));
  g[0] = 0;
  a0 = a;
  for (k = 0; k<K; k++){
    iXn = 0;
    for (n = 0; n<N; n++){
      iXs = 0;
      a = a0; // Reset a to the top for this k-index.
      for (s = 0; s<S; s++){
	*(G++) = 0;
	g[1] = *a;
	ar   = *(a+1);
	// Compute the filter shape
	for (t = 2; t<T; t++) g[t] = g[t-1]*ar;
	for (t = 2; t<T; t++) g[t] = g[t]*((double)t); // Compute the convolution
	for (t = 1; t<T; t++){
	  for (tp = 1; tp<=t; tp++)
	    (*G) += X[(t - tp) + iXs + iXn]*g[tp];
	  G++;
	}
	iXs += T;
	a += 2;
      }
      iXn += T*S;
    }
    a0 = a;
  }
}
