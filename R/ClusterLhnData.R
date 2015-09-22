#' A function for fitting Shahar's LHN data to the linear-nonlinear-poisson
#' model.
#'
#' In the following, T: the number of time bins, S: the number of odors, N: the
#' number of cells, K: the number of clusters.
#'
#' @param X A binary T x S x N array indicating the odor window for each odor
#'   and cell.
#' @param Y An integer valued T x S x N array containing the spike counts in
#'   each bin.
#' @param ainit A 2 x S x N containing the parameters of the fitted drives for
#'   each cell and odor.
#' @param numClusters The number of clusters to use.
#' @param kalpha The shape parameter for the gamma prior on alpha.
#' @param thalpha The scale parameter for the gamma prior on alpha
#' @param tauv0 The rate constant for the exponential prior on v0.
#' @param taua The rate constant for the exponential prior on the drive
#'   parameters
#' @param taul0 The rate constant for the exponential prior on l0.
#' @param numIters The maximum number of iterations.
#' @param dt The time constant of the updates.
#' @param seed The random seed to use.
#' @param initMode The initialization mode for the clustering. Can be "random",
#'   "kmeans", or "kmeans++".
#' @param verbose If TRUE will print out the progress of the algorithm and other
#'   diagonstic information.
#' @param timer If "ON" will time different blocks of the code.
#' @param slopeRatioToStop The ratio of the rate of change of the objective at
#'   the end to the start above which to terminate.
#' @param numSlopePoints How many points to take to compute the slope of the
#'   objective
#' @param checkToStopEvery How often to compute the stopping ratio.
#' @return A list consisting of
#'
#'   \item{seed}{The random seed used.}
#'
#'   \item{a}{A 2 x S x K array containing the learned drive parameters}
#'
#'   \item{al}{A N x 1 vector containing the learned alpha values}
#'
#'   \item{v0}{A N x 1 vector containing the learned v0 parameters}
#'
#'   \item{l0}{A N x 1 vector containing the learned l0 parameters}
#'
#'   \item{qnk}{A N x K matrix of the cluster responsibilities for each data
#'   point.}
#'
#'   \item{L}{A T x S x N x K array containing the final lambda values}
#'
#'   \item{numIters}{The actual number of iterations that ran.}
#'
#'   \item{F}{A numIters x 1 array containing the objective function as function
#'   of the number of iterations.}
#'
#'   \item{clust}{A N x 1 vector of cluster assignments.}
#'
#'   \item{pclust}{A N x 1 vector of the probabilities of the cluster chosen.}
#'
#'   \item{dclust}{A N x 1 vector of distances to its cluster center.}
#'
#'   \item{exitMode}{A string with the exit mode of the algorithm: "ITERS" if it
#'   hit the maximum number of iterations, "SLOPE_RATIO" if it exited early due
#'   to the slope ratio.}
#' @export
ClusterLhnData <- function(X, Y, ainit = NULL, numClusters=3, kalpha=10, thalpha=3/20, tauv0 = 0.1, taua=1, taul0=0.5, minIters = 0, numIters=10000, dt=1e-5, seed=0, initMode="random", verbose=TRUE, timer="OFF", slopeRatioToStop=100, numSlopePoints=20, checkToStopEvery=100, keepHistory = NULL){
  set.seed(seed);
  Timers = TIMER_INIT(status=timer);

  ## If X and Y have data for just 1 cell, add a singleton dimension so the rest of the code doesn't complain.
  if (length(dim(X))==2) X = array(X,dim = c(dim(X),1));
  if (length(dim(Y))==2) Y = array(Y,dim = c(dim(Y),1));

  T = dim(X)[[1]]; # Number of bins
  S = dim(X)[[2]]; # Number of odors
  N = dim(X)[[3]]; # Number of cells
  K = numClusters; # Number of clusters

  if (K>N)
      stop(sprintf("Number of clusters %d is greater than number of data points %d.", K, N));

  SCALE_GRAD <- function(g){
    ng = sqrt(sum(g^2));
    return(g/max(1,ng));
  }

  DFUN1 <- function(a1,a2) sum((a1[[1]]*(a1[[2]]^(0:(T-1))) - a2[[1]]*(a2[[2]]^(0:(T-1))))^2);
  DFUN  <- function(A1,A2) sum(sapply(1:K, function(i) DFUN1(A1[,i],A2[,i])));

  Y = array(rep(Y,K), dim=c(dim(Y), K));
  lfY= log(factorial(Y));

  l0     = matrix(rexp(N,rate=1/taul0),        nrow=N);
  v0     = matrix(rexp(N,rate=1/tauv0),        nrow=N);
  al     = matrix(rgamma(N, shape=kalpha, scale=thalpha), nrow=N);
  F      = matrix(NaN, nrow=numIters, ncol=1);

  ## Initialize the clusters
  a = array(0,c(2,S,K));
  if (initMode=="single"){
      if (K != 1)
          stop("Number of clusters must equal 1 when fitting a single data point.");
      a[1,,] = 1;
      a[2,,] = 0.1;
  }else if (initMode=="random" || is.null(ainit)){
    a = array(runif(n = 2*K*S), dim=c(2,S,K));
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
  ## Initialize the history object
  history = NULL;

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

    args = c(list(history, !is.null(keepHistory), t, numIters), keepHistory);
    history = do.call("HISTORY", args);

    ## Break out before updating if you're on the last iteration.
    ## To keep the parameters and results current.
    if (t>=minIters){
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
            cat(sprintf("\nBreaking at t = %d due to SLOPE_RATIO %8.3f > %8.3f\n", t, ratio, slopeRatioToStop));
          break;
        }
      }
    }

    ## M-STEP: Update the parameters
    ## 1. Compute the gradients
    Timers <- TIMER_TIC("M_GRADIENTS", Timers);
    Z   = Y/L - 1;
    Zq  = Z*qtsnk;
    ZqV = Zq*(V>0);
    alb = aperm(array(rep(al,T*S*K),dim=c(N,T,S,K)),c(2,3,1,4));
    ZqVa = ZqV*alb;

    gradL0 =  apply(Zq,     3, FUN="sum");
    gradAl =  apply(ZqV*V,  3, FUN="sum") + ((kalpha - 1)/al - 1/thalpha); ## Add a prior on Al
    gradV0 = -apply(ZqV,    3, FUN="sum")*al;

    grada = apply(ZqVa*U1, c(2,4),FUN="sum");
    G     = ComputeGtsnk(X,a);
    gradr = apply(ZqVa*G, c(2,4),FUN="sum");
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

  if (length(clust)>1){
      dclust = -diag(apply(LL[,,,clust],c(3,4),FUN="sum"));
  }
  else{
      dclust = -sum(LL);
  }

  Lclust = array(data=0, dim=c(T,S,N));
  for (i in 1:N)
      Lclust[,,i] = L[,,i,clust[[i]]];
  results= list(a = a, al = al, v0 = v0, l0 = l0, qnk = qnk, F = F[1:t], L = L, Lclust = Lclust, clust=clust, pclust=pclust, dclust=dclust, numIters = t, seed = seed, exitMode = exitMode, history = history);
}
