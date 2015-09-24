#include <R.h>
int ComputeGlobalIndexFromLocal(int il, int n, int *dl, int *dg);

void BroadcastMultiply(double *X, double *y, int *dims, int *pndims, int *margins, int *pnmargins, double *Xy, int *computeInds, int *ii, int *pni, int *ir, int *pnr){
  // Multiplies X by y along the specified margins, broadcasting out to the remaining dimensions.
  //
  // Algorithm:
  // Indices into X can be decomposed as
  // i = i0 + d1 i1 + d1 d2 i2 + ... + prod(d_i,1,k) ik
  //
  // To broadcast multiply, we first split the coordinates into those
  // specified in the margins (I), and those in the remainder (R). Let
  // #In be the size of that dimension, and similarly Rn.
  //
  // What we want to do is, for a given (i1,i2,...), to multiply all
  // elements of X with coordinates (i1,i2,..|r1,r2,..) with y(i1,i2,...).
  //
  // We can do this by:
  // FOR i in 1:length(y)
  // 1. Map i to coordinates i => (i1,i2...) => OFFSET  
  // 2. FO (r in 1:#elelments in other dims)
  // 3.     Map r to coordinates (r1,r2,...)
  //        Map (r1,r2,...) to a global index => INDEX
  //        Multiply X[OFFSET + INDEX] by y[i];
  // We can compute the offsets and indices before hand.

  int i, j, nrdims, ni, nr, ji, jr, ndims = *pndims, nmargins = *pnmargins;
  int *dg, *di, *dli, *dgi, *dr, *dlr, *dgr, *imargin;

  if (*computeInds>0){
    nrdims = ndims - nmargins;
    dg  = (int *)R_alloc(ndims,    sizeof(int)); // The global coefficients overall (1, #(1), #(1)#(2), ...)

    di  = (int *)R_alloc(nmargins, sizeof(int)); // The dimensions of the margins:    (#I1, #I2, ...)
    dli = (int *)R_alloc(nmargins, sizeof(int)); // The local coefs for the margins: (1, #I1, #I1 #I2, ...)
    dgi = (int *)R_alloc(nmargins, sizeof(int)); // The global coefs for the margins: (#(...I1), #(...I2), ...)

    dr  = (int *)R_alloc(nrdims, sizeof(int)); // The dimensions of the remaining:    (#R1, #R2, ...)
    dlr = (int *)R_alloc(nrdims, sizeof(int)); // The local coefs for the remaining:  (1, #R1, #R1 #R2, ...)
    dgr = (int *)R_alloc(nrdims, sizeof(int)); // The global coefs for the remaining: (#(...R1), #(...R2), ...)

    // Margin indicator
    imargin = (int *)R_alloc(ndims, sizeof(int));
    for (i = 0; i<ndims; i++)    imargin[i] = 0;
    for (i = 0; i<nmargins; i++) imargin[margins[i]] = 1;

    dg[0] = 1;
    for (i = 1; i<ndims; i++) dg[i] = dg[i-1]*dims[i-1];

    ji = 0; jr = 0;
    ni = 1; nr = 1;
    for (i = 0; i<ndims; i++)
      if (imargin[i]){
	di[ji]     = dims[i];
	ni        *= dims[i];
	dgi[ji]    = dg[i];
	ji++;
      }
      else{
	dr[jr]     = dims[i];
	nr        *= dims[i];
	dgr[jr]    = dg[i];
	jr++;
      }
  
    dli[0] = 1;
    for (i=1; i<nmargins; i++)
      dli[i] = dli[i-1]*di[i-1];
  
    dlr[0] = 1;
    for (i=1; i<nrdims; i++)
      dlr[i] = dlr[i-1]*dr[i-1];
  
    // Precompute the OFFSETS
    for (i = 0; i<ni; i++)
      ii[i] = ComputeGlobalIndexFromLocal(i, nmargins, dli, dgi);
  
    // Precompute the INDICES
    for (i = 0; i<nr; i++)
      ir[i] = ComputeGlobalIndexFromLocal(i, nrdims, dlr, dgr);

    *pni = ni;
    *pnr = nr;
  }
  else{
    ni = *pni;
    nr = *pnr;
  }
    
  // Actually perform the multiplication
  for (i = 0; i<ni; i++)
    for (j = 0; j<nr; j++)
      Xy[ii[i] + ir[j]] = X[ii[i] +ir[j]]*y[i];

}
