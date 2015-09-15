#' Title
#'
#' @param X dfgkdf;g asdfsadf
#' @param Y
#' @param Fits
#' @param numClusters
#' @param kalpha
#' @param thalpha
#' @param tauth
#' @param taua
#' @param tauu0
#' @param filterLength
#' @param be
#' @param numIters
#' @param fvalTol
#' @param dt
#' @param seed
#' @param emptyAction
#' @param initMode
#' @param verbose
#' @param timer
#' @param slopeRatioToStop
#' @param numSlopePoints
#' @param checkToStopEvery
#' @export
#'
ClusterLhnData <- function(X, Y, Fits, numClusters=3, kalpha=10, thalpha=3/20, tauth = 0.1, taua=1, tauu0=0.5, filterLength=20, be=5, numIters=10000, fvalTol=1e-2, dt=1e-5, seed=0, emptyAction="singleton", initMode="random", verbose=TRUE, timer="OFF", slopeRatioToStop=100, numSlopePoints=20, checkToStopEvery=100){
  set.seed(seed);
  Timers = TIMER_INIT(status=timer);

  T = dim(X)[[1]]; # Number of bins
  S = dim(X)[[2]]; # Number of odors
  N = dim(X)[[3]]; # Number of cells
  K = numClusters; # Number of clusters

  SCALE_GRAD <- function(g){
    ng = sqrt(sum(g^2));
    return(g/max(1,ng));
  }

  DFUN1 <- function(a1,a2) sum((a1[[1]]*(a1[[2]]^(0:(T-1))) - a2[[1]]*(a2[[2]]^(0:(T-1))))^2);
  DFUN  <- function(A1,A2) sum(sapply(1:K, function(i) DFUN1(A1[,i],A2[,i])));

  Y = array(rep(Y,K), dim=c(dim(Y), K));
  lfY= log(factorial(Y));

  l0     = matrix(rexp(N,rate=1/tauu0),        nrow=N);
  v0     = matrix(rexp(N,rate=1/tauth),        nrow=N);
  al     = matrix(rgamma(N, shape=kalpha, scale=thalpha), nrow=N);
  F      = matrix(NaN, nrow=numIters, ncol=1);

  ## Initialize the clusters
  ainit = array(unlist(Fits[,"aest"]), c(2,S,N), dimnames = list(list("aa","ar"), dimnames(X[[2]]), dimnames(X[[3]])));
  a = array(0,c(2,S,K));

  if (initMode=="random"){
    a = array(runif(n = 2*K*S), dim=c(2,K,S), dimnames = list(list("aa","ar"), dimnames(X[[3]]), sapply(1:numClusters, function (i) {sprintf("cluster%d", i)})));
  }else if (initMode=="kmeans"){
    isample = sample(N,size=S,replace=FALSE);
    a      = ainit[,,isample];
  }else if (initMode =="kmeans++"){
    ## Pick the first cluster center at random from the data
    iclust = sample(N,1);
    a[, ,1] = ainit[,,iclust];
    if (K>1)
      for (i in 2:K){
        iavail = setdiff(1:N,iclust);
        ## Compute the distance from each (non-cluster center) datapoint to the cluster centers
        d = c();
        for (j in 1:length(iavail)){
          dk = c();
          a0  = ainit[,,iavail[[j]]];
          for (k in 1:length(iclust))
            dk[[k]] = DFUN(a0,ainit[,,iclust[[k]]]);
          d[[j]] = min(dk);
        }

        p = d/sum(d);
        if (length(p)==1){
          iclust[[i]] = iavail;
        }else{
          iclust[[i]] = sample(iavail, 1, prob=p);
        }
        a[,,i] = ainit[,,iclust[[i]]];
      }
  }else{
    stop(sprintf("Unknown initMode '%s'.\n", initMode));
  }

  ## Begin the iterations
  if (verbose) pb = txtProgressBar(min=1,max=numIters,initial=1,style=3);
  Timers <- TIMER_TIC("ALL_ITERS", Timers);
  for (t in 1:numIters){
    Drive = ComputeLambda(X,a,al,v0,l0);

    L = Drive$L; L[L==0] = 1e-12;
    U = Drive$U;
    U1= Drive$U1;
    V = Drive$V;

    LL = -L + Y*log(L) - lfY;
    ## E-STEP: Compute the cluster responsibilites
    Timers <- TIMER_TIC("E_STEP", Timers);
    lpnk = apply(LL, c(3,4), FUN="sum");
    lpnk = sweep(lpnk, 1, apply(lpnk,1,FUN="max"));
    qnk  = exp(lpnk);
    qnk  = sweep(qnk, 1, apply(qnk,1,FUN="sum"), FUN="/");

    qtsnk = aperm(array(rep(qnk,S*T),dim=c(N,K,T,S)), c(3,4,1,2));
    ## Compute the objective function
    H = -qnk*log(qnk); H = sum(H[!is.nan(H)]); # Entropy
    Eqll = sum(qtsnk * LL); # Expected log likelihood
    lprAl = sum((kalpha-1)*log(al) - al/thalpha); ## Prior on alpha
    F[[t]] = H + Eqll + lprAl;
    Timers <- TIMER_TOC("E_STEP", Timers);

    ## Break out before updating if you're on the last iteration.
    ## To keep the parameters and results current.
    if (t == numIters){
      exitMode = "ITERS";
      break;
    }else if ((t %% checkToStopEvery)==0){
      x = 1:numSlopePoints;
      y0 = F[x];
      mstart = lm(y0 ~ x);
      y1 = F[t + ((-numSlopePoints+1):0)];
      mend = lm(y1 ~ x);
      ratio = abs(mstart$coefficients[2]/mend$coefficients[2]);
      if (ratio > slopeRatioToStop){
        exitMode= "SLOPE_RATIO";
        if (verbose)
          cat(sprintf("Breaking at t = %d due to SLOPE_RATIO %8.3f > %8.3f\n", t, ratio, slopeRatioToStop));
        break;
      }
    }

    ## M-STEP: Update the parameters
    ## 1. Compute the gradients
    Timers <- TIMER_TIC("M_GRADIENTS", Timers);
    Z   = Y/L - 1;
    Zq  = Z*qtsnk;
    ZqV = Zq*(V>0);

    gradL0 =  apply(Zq,     3, FUN="sum");
    gradAl =  apply(ZqV*V,  3, FUN="sum") + (kalpha - 1)/al - 1/thalpha; ## Add a prior on Al
    gradV0 = -apply(ZqV,    3, FUN="sum")*al;

    grada = apply(Zq*U1, c(2,4),FUN="sum");

    G     = ComputeGtsnk(X,a);
    gradr = apply(Zq*G, c(2,4),FUN="sum");
    Timers <- TIMER_TOC("M_GRADIENTS", Timers);

    ## 2. Update the parameters
    Timers <- TIMER_TIC("M_UPDATES", Timers);
    l0 = l0 + SCALE_GRAD(gradL0)*dt; l0[l0<0] = 0;
    al = al + SCALE_GRAD(gradAl)*dt; al[al<0] = 0;
    v0 = v0 + SCALE_GRAD(gradV0)*dt;

    aa = a[1,,] + SCALE_GRAD(grada)*dt; aa[aa<0] = 0; a[1,,] = aa;
    ar = a[2,,] + SCALE_GRAD(gradr)*dt; ar[ar<0] = 0; ar[ar>1] = 1; a[2,,] = ar;
    Timers <- TIMER_TOC("M_UPDATES", Timers);

    if (verbose) setTxtProgressBar(pb, t);
  }
  Timers <- TIMER_TOC("ALL_ITERS", Timers);
  TIMER_SUMMARY(Timers);
  if (verbose) close(pb);

  ## Extract a clustering from the probabilities
  clust =  apply(qnk, 1, "which.max");
  pclust =  apply(qnk, 1, "max");
  dclust = -diag(apply(LL[,,,clust],c(3,4),FUN="sum"));
  results= list(a = a, al = al, v0 = v0, l0 = l0, qnk = qnk, F = F[1:t], L = L,clust=clust, pclust=pclust, dclust=dclust, numIters = t, seed = seed, exitMode = exitMode);
}
