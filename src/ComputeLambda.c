#include <R.h>

void ComputeLambda(double *X, double *a, int *dim,
		   double *alpha, double *v0, double *l0,
		   double *U, double *U1, double *V, double *L){
  // X is in tsn format
  // a is in 2sk format
  // U,L are in tsnk format.
  int t, s, n, k;
  int T = dim[0]; // # of time bins
  int S = dim[1]; // # of odors
  int N = dim[2]; // # of cells
  int K = dim[3]; // # of clusters

  for (k = 0; k<K; k++){
    for (n = 0; n<N; n++){
      for (s = 0; s<S; s++){
	*U = (*X)*(*a);
	*U1= (*X); // The case for a = 1.
	*V = (*U - *v0);
	*L = *l0 + (*V > 0 ? (*V)* (*alpha) : 0);
	U++; U1++; V++; L++; X++;
	for (t = 1; t<T; t++){
	  *U = (*X)*(*a) + (*(U  - 1))*(*(a+1));
	  *U1= (*X)      + (*(U1 - 1))*(*(a+1));
	  *V = (*U - *v0);
	  *L = *l0 + (*V > 0 ? (*V)*(*alpha) : 0);
	  U++; U1++; V++; L++; X++;
	}
	a += 2;
      }
      v0++;
      l0++;
      alpha++;
      a-=2*S;
    }
    X     -= T*S*N;
    alpha -= N;
    v0    -= N;
    l0    -= N;
    a     += 2*S;
  }
}
