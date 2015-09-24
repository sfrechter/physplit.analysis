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

void ApplySum(double *X, int *dims, int *pndims,  int *margins, int *pnmargins, double *Y){
  // Computes the sum of the elements of X over all non-margin dimensions.
  //
  // Algorithm:
  // Indices into X can be decomposed as
  // i = i0 + d1 i1 + d1 d2 i2 + ... + prod(d_i,1,k) ik
  //
  // To marginalize, we hold indices in the margin dimensions fixed,
  // while we vary those in the remaining dimensions.
  //
  // Call the set of margins to keep I and the set to marginalize over
  // R (for the remaining dims). Let #In be the size of that dimension,
  // and similarly Rn.
  //
  // What we want to do is, for a given (i1,i2,...), to sum over (r1,r2,...).
  //
  // The index (i1,i2,...) into the output array is easy to compute:
  // iy = i1 + #I1 i2 + #I1 #I2 i2 + ...
  //
  // The index into the input matrix for (i1,i2,... | r1, r2,...) can
  // be written as the sum of an offset (due to i) and a running
  // index (due to r). The offset is simply
  //
  // OFFSET = #(...I1) i1 + #(...I2) i2 + ...
  //
  // where #(...Ik) = prod(#Ii, 1, k)
  //
  // INDEX = #(...R1) r1 + #(...R2) r2 + ...
  //
  // We'll be indexing linearly through i and r, so we'll need to have
  // functions that map i to OFFSET and r to index. This can be done
  // in two steps: a mapping into a 'local' set of coordinates, and
  // then the mapping into a 'global', linear index:
  //
  // i -> (i1,i2,...) -> OFFSET
  //
  // To perform the first mapping, we start with i, and iteratively
  // 1. i  = i/Ik
  // 2. ik = mod(i, Ik);
  // 3. i = (i - ik);
  // Until i = 0.
  //
  // Once we have the ik, mapping to th offset is easy.  But we need
  // two sets of numbers: (1, #I1, #I1 #I2, ...) to compute the first
  // mapping, and (#(...I1), #(...I2),...) to compute the second.
  int i, j, nrdims, ni, nr, ji, jr, ndims = *pndims, nmargins = *pnmargins;
  int *dg, *di, *dli, *dgi, *dr, *dlr, *dgr, *ii, *ir, *imargin;

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
  ii = (int *)R_alloc(ni, sizeof(int));
  for (i = 0; i<ni; i++)
    ii[i] = ComputeGlobalIndexFromLocal(i, nmargins, dli, dgi);
  
  // Precompute the INDICES
  ir = (int *)R_alloc(nr, sizeof(int));
  for (i = 0; i<nr; i++)
    ir[i] = ComputeGlobalIndexFromLocal(i, nrdims, dlr, dgr);

  // Actually perform the sum
  for (i = 0; i<ni; i++){
    Y[i] = 0;
    for (j = 0; j<nr; j++)
      Y[i] += X[ii[i] + ir[j]];    
  }
}
