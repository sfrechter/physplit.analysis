#include <R.h>

void ComputeQ(double *LL, int *d,
	      double *qnk, double *qtsnk){  
  int T = d[0], S = d[1], N = d[2], K = d[3];
  int i, in, ik;
  double *lpnk = (double *)R_alloc(N*K, sizeof(double));
  double *lmax = (double *)R_alloc(N, sizeof(double));
  double *qsum = (double *)R_alloc(N, sizeof(double));

  // Initialize
  for (in = 0; in<N; in++){
    *lmax++ = -INFINITY;
    *qsum++ = 0;
  }
  lmax -= N;
  qsum -= N;
    
  // Compute lpnk and the maximum across k
  for (ik = 0; ik<K; ik++){
    for (in = 0; in<N; in++){
      *lpnk = 0;
      for (i = 0; i<T*S;i++)
	*lpnk += *LL++;
      if (*lpnk>*lmax) *lmax = *lpnk;
      lpnk++;
      lmax++;
    }
    lmax -= N;
  }
  
  // Subtract out lmax
  lpnk -= N*K;
  for (ik = 0; ik<K; ik++){
    for (in = 0; in<N; in++)
      *lpnk++ -= *lmax++;
    lmax -= N;
  }
  
  // Compute qnk = exp(lpnk) and figure out the sums for each n.
  lpnk -= N*K;
  for (ik = 0; ik<K; ik++){
    for (in = 0; in<N; in++){
      *qnk  = exp(*lpnk++);
      *qsum++ += *qnk++;
    }
    qsum -= N;
  }

  // Normalize qnk so that sums over k all equal 1.
  qnk -= N*K;
  for (ik = 0; ik<K; ik++){
    for (in = 0; in<N; in++)
      *qnk++/=*qsum++;
    qsum-=N;
  }  

  // Now expand qnk out to qtsnk
  qnk -= N*K;
  for (ik = 0; ik<K; ik++)
    for (in = 0; in<N; in++){
      for (i = 0; i<T*S; i++)
	*qtsnk++ = *qnk;
      qnk++;
    }
  
}
