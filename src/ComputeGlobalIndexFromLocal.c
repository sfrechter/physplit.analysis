#include <R.h>
int ComputeGlobalIndexFromLocal(int il, int n, int *dl, int *dg){
  // Computes a global index from the local one provided by first
  // converting the index into local coordinates, and then converting
  // those coordiantes into a global index.
  int i, ig = 0;
  int *cl = (int *)R_alloc(n, sizeof(int)); // The local coordinates
  int d;
  if (n>1){
    for (i = 0; i<n-1; i++){
      d     =  dl[i+1]/dl[i];
      cl[i] = (il % d);
      il   -=  cl[i];
      il   /=  d;
      ig   +=  cl[i]*dg[i];
    }
  }
  cl[n-1] = il;
  ig += cl[n-1]*dg[n-1];
  return ig;
}
