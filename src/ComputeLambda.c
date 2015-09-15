#include <R.h>

void ComputeLambda(double *X, double *a, int *dim,
		   double *alpha, double *v0, double *l0,
		   double *U, double *U1, double *V, double *L){
  // X is in tsn format
  // a is in 2sk format
  // U,L are in tskn format.
  int t, s, n, k;
  int T = dim[0]; // # of time bins
  int S = dim[1]; // # of odors
  int N = dim[2]; // # of cells
  int K = dim[3]; // # of clusters
  double *X0 = X; // Store the location of the top of X
  double *a0 = a;

  for (n = 0; n<N; n++){
    a = a0; // Reset to the top of the cell.
    for (k = 0; k<K; k++){
      X = X0; // Reset X for each cluster.
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
	a+=2;
      }
    }
    X0 += T*S; // Move to the next cell.
    alpha++;
    v0++;
    l0++;
  }
}
