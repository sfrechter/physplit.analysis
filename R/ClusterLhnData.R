#' A function for fitting Shahar's LHN data to the linear-nonlinear-poisson
#' model.
#'
#' In the following, T: the number of time bins, S: the number of odors, N: the
#' number of cells, K: the number of clusters.
#' @param Data FIXME - what should the input data look like?!
#' @param numClusters The number of clusters to use.
#' @param kalpha The shape parameter for the gamma prior on alpha.
#' @param thalpha The scale parameter for the gamma prior on alpha
#' @param sdv0 The standard deviation of the gaussian prior on membrane potential offset.
#' @param taua The rate constant for the exponential prior on the drive
#'   parameters
#' @param taul0 The rate constant for the exponential prior on l0.
#' @param minIters The minumum number of iterations.
#' @param numIters The maximum number of iterations.
#' @param dt The time constant of the updates.
#' @param seed The random seed to use.
#' @param initMode The initialization mode for the clustering. Can be "random",
#'   "kmeans", or "kmeans++".
#' @param iclust An initial clustering assignment, if any.
#' @param verbose If TRUE will print out the progress of the algorithm and other
#'   diagonstic information.
#' @param timer If "ON" will time different blocks of the code.
#' @param slopeRatioToStop The ratio of the rate of change of the objective at
#'   the end to the start above which to terminate.
#' @param numSlopePoints How many points to take to compute the slope of the
#'   objective
#' @param checkToStopEvery How often to compute the stopping ratio.
#' @param keepHistory A least of strings containing the variables to track.
#' @param keepHistoryAt A list of iterations at which to record history. If NULL defaults to all.
#' @param maxPreFitIters The maximum number of iterations to pre fit the cell-specific parameters to the clusters. If set to 0 will not prefit the parameters.
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
#'   \item{Lclust}{A T x S x N array containing the final lambda value for the most likely cluster for each fit.}
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
#'
#'   \item{history}{A list containing the values of the tracked variables for the specified iterations.}
#'
#'   \item{misc}{A miscellaneous list to hold other variables, used mainly for debugging.}
#' @export
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @importFrom stats runif
ClusterLhnData <- function(Data, numClusters=3, kalpha=10, thalpha=3/20, sdv0 = 0.1, taua=1, taul0=0.5, minIters = 0, numIters=10000, dt=1e-5,
                           seed=0, initMode="random", iclust = NULL,verbose=TRUE, timer="OFF",
                           slopeRatioToStop=100, numSlopePoints=20, checkToStopEvery=100,
                           keepHistory = NULL, keepHistoryAt = NULL, maxPreFitIters = 1){
    Timers = TIMER_INIT(status=timer);

    X = Data$X;
    Y = Data$Y;

    ## If X and Y have data for just 1 cell, add a singleton dimension so the rest of the code doesn't complain.
    if (length(dim(X))==2) X = array(X,dim = c(dim(X),1));
    if (length(dim(Y))==2) Y = array(Y,dim = c(dim(Y),1));

    T = dim(X)[[1]]; # Number of bins
    S = dim(X)[[2]]; # Number of odors
    N = dim(X)[[3]]; # Number of cells
    K = numClusters; # Number of clusters

    if (K>N)
        stop(sprintf("Number of clusters %d is greater than number of data points %d.", K, N));

    if (length(seed)==1){
        set.seed(seed);
    }else if (length(seed) == N){
         initMode = "fixed";
     }else{
          stop("Seed must either be an integer or a vector with NUM_CELLS elements.");
      }

    SCALE_GRAD <- function(g){
        ng = sqrt(sum(g^2));
        return(g/max(1,ng));
    }

    DFUN1 <- function(a1,a2) sum((a1[[1]]*(a1[[2]]^(0:(T-1))) - a2[[1]]*(a2[[2]]^(0:(T-1))))^2);
    DFUN  <- function(A1,A2) sum(sapply(1:K, function(i) DFUN1(A1[,i],A2[,i])));

    Y = array(rep(Y,K), dim=c(dim(Y), K));
    # lfY= log(factorial(Y));

    if (initMode=="single"){
        l0     = matrix(apply(Y,MARGIN=3,FUN=mean),     nrow=N);
        v0     = matrix(0, nrow=N, ncol = 1);
        al     = matrix(1, nrow=N, ncol=1);
    }else{
         ainit  = Data$afit;
         l0     = matrix(sapply(Data$Fits, FUN=function(x) x$l0), nrow=N, ncol=1);
         v0     = matrix(sapply(Data$Fits, FUN=function(x) x$v0), nrow=N, ncol=1);
         al     = matrix(sapply(Data$Fits, FUN=function(x) x$al), nrow=N, ncol=1);
     }
    F      = matrix(NaN, nrow=numIters, ncol=1);

    misc = NULL; ## List to contain miscellaneous information.
    ## Initialize the clusters
    a = array(0,c(2,S,K));
    if (initMode=="single"){
        if (K != 1)
            stop("Number of clusters must equal 1 when fitting a single data point.");
        a[1,,] = 1;
        a[2,,] = 0.1;
    }else if (initMode=="fixed"){
         cat(sprintf("initMode = specifed.\n"));
         a[1,,] = 1;
         a[2,,] = 0.1;
         qnk = matrix(0, nrow = N, ncol = K);
         qtsnk = array(0, dim = c(T,S,N,K));
         for (i in 1:N){
             for (j in 1:K){
                 qnk[[i,j]]   = as.double(seed[[i]] == j);
                 qtsnk[,,i,j] = as.double(seed[[i]] == j);
             }
         }
     }else if (initMode=="random"){
         a = array(runif(n = 2*K*S), dim=c(2,S,K));
     }else if (initMode=="kmeans"){
          if (is.null(iclust)){
              isample = sample(N,size=K,replace=FALSE);
          }else{
               isample = iclust;
           }
          a = ainit[,,isample];
          misc$iclust = isample;
          if (maxPreFitIters > 0){
              ifit = setdiff(1:N, isample);
              minFitIters=1000;
              numFitIters=10000;
              dtFit = 1e-2;
              if (verbose)
                  cat("Finding parameters for non-cluster-center cells...\n");
              l0f = NULL;
              alf = NULL;
              v0f = NULL;
              Fbest = NULL;
              converged = FALSE;
              fitIters = 1;
              while (!converged && fitIters <= maxPreFitIters){
                  if (verbose)
                      cat(sprintf("Fit iteration %d\n", fitIters));
                  res = FitCellSpecificParameters(Data$X[,,ifit], Data$Y[,,ifit], a, l0 = l0f, al = alf, v0 = v0f, kalpha=kalpha, thalpha=thalpha, sdv0 = tauv0, taua=taua, taul0=taul0, minIters = minFitIters, numIters=numFitIters, dt=dtFit, seed = seed, verbose=verbose, timer="OFF", slopeRatioToStop=500, numSlopePoints=20, checkToStopEvery=100, keepHistory = NULL, keepHistoryAt = NULL);
                  if (!is.null(Fbest)){
                      improved = which((res$Fbest - Fbest)/abs(Fbest)>0.001);
                      if (length(improved)>0){
                          Fbest[improved] = res$Fbest[improved];
                          if (verbose){
                              cat(sprintf("%d cells improved.\n", length(improved)));
                              print(improved);
                          }
                          l0f[improved,]=res$l0.best[improved];
                          v0f[improved,]=res$v0.best[improved];
                          alf[improved,]=res$al.best[improved];
                      }else{
                           converged = TRUE;
                       }
                  }else{
                       Fbest = res$Fbest;
                       l0f = matrix(res$l0.best, nrow = N - K, ncol = K);
                       v0f = matrix(res$v0.best, nrow = N - K, ncol = K);
                       alf = matrix(res$al.best, nrow = N - K, ncol = K);
                   }

                  fitIters = fitIters + 1;
              }
              misc$ifit  = ifit;
              misc$res = res;
              l0[ifit] = res$l0.best;
              al[ifit] = res$al.best;
              v0[ifit] = res$v0.best;
          }
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

    ## Indices for marginalizing out all dims but the third.
    ii3 = NULL; ir3 = NULL;
    ii24= NULL; ir24= NULL;

    ## Initialize the clustering
    clust = matrix(NA,nrow=N,ncol=1);
    ## Begin the iterations
    if (verbose) pb = txtProgressBar(min=0,max=numIters,initial=1,style=3);
    Timers <- TIMER_TIC("ALL_ITERS", Timers);
    startTime = proc.time();
    for (t in 1:numIters){
        Timers <- TIMER_TIC("COMPUTE_LAMBDA", Timers);
        Drive = ComputeLambda(X,a,al,v0,l0);
        L = Drive$L; L[L==0] = 1e-12;
        U = Drive$U;
        U1= Drive$U1;
        V = Drive$V;
        Timers <- TIMER_TOC("COMPUTE_LAMBDA", Timers);
        LL = -L + Y*log(L);
        ## E-STEP: Compute the cluster responsibilites
        Timers <- TIMER_TIC("E_STEP_Q", Timers);
        if (initMode != "fixed"){
            Q     = ComputeQ(LL);
            qnk   = Q$qnk;
            qtsnk = Q$qtsnk;
        }
        Timers <- TIMER_TOC("E_STEP_Q", Timers);

        ## Compute a clustering
        clust[,1]  =  apply(qnk, 1, "which.max");

        ## Compute the objective function
        Timers <- TIMER_TIC("E_STEP_SCALARS", Timers);
        H = -qnk*log(qnk); H = sum(H[!is.nan(H)]); # Entropy
        Eqll = sum(qtsnk * LL); # Expected log likelihood
        lprAl = sum((kalpha-1)*log(al) - al/thalpha); ## Prior on alpha
        F[[t]] = H + Eqll + lprAl;
        Timers <- TIMER_TOC("E_STEP_SCALARS", Timers);

        Timers <- TIMER_TIC("HISTORY", Timers);
        args    = c(list(history, !is.null(keepHistory), t, numIters, at=keepHistoryAt), keepHistory);
        history = do.call("HISTORY", args);
        Timers <- TIMER_TOC("HISTORY", Timers);

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
        Timers <- TIMER_TIC("M_STEP_Z", Timers);
        Z   = Y/L - 1;
        Zq  = Z*qtsnk;
        ZqV = Zq*(V>0);
        if (is.null(ii3)){
            res = BroadcastMultiply(ZqV, al, 3, ii=ii3, ir=ir3);
            ii3 = res$ii; ir3 = res$ir;
        }
        ZqVa = BroadcastMultiply(ZqV, al, 3, ii=ii3, ir=ir3)$Xy;   ## alb = aperm(array(rep(al,T*S*K),dim=c(N,T,S,K)),c(2,3,1,4));
        Timers <- TIMER_TOC("M_STEP_Z", Timers);

        Timers <- TIMER_TIC("M_STEP_GRADS", Timers);
        gradL0 =  ApplySum(Zq,    3, ii=ii3, ir=ir3)$Y;
        gradAl =  ApplySum(ZqV*V, 3, ii=ii3, ir=ir3)$Y + ((kalpha - 1)/al - 1/thalpha); ## Add a prior on Al.
        gradV0 = -ApplySum(ZqV,   3, ii=ii3, ir=ir3)$Y*al;

        if (is.null(ii24)){
            res  = ApplySum(ZqVa*U1, c(2,4), ii=ii24, ir=ir24);
            ii24 = res$ii; ir24 = res$ir;
        }

        grada = ApplySum(ZqVa*U1, c(2,4), ii=ii24, ir=ir24)$Y;
        G     = ComputeGtsnk(X,a);
        gradr = ApplySum(ZqVa*G,  c(2,4), ii=ii24, ir=ir24)$Y;
        Timers <- TIMER_TOC("M_STEP_GRADS", Timers);

        ## 2. Update the parameters
        Timers <- TIMER_TIC("M_STEP_UPDATES", Timers);
        l0 = l0 + SCALE_GRAD(gradL0)*dt; l0[l0<1e-6] = 1e-6;
        if (initMode != "single"){
            al = al + SCALE_GRAD(gradAl)*dt; al[al<0] = 0;
            v0 = v0 + SCALE_GRAD(gradV0)*dt;
        }

        aa = a[1,,] + SCALE_GRAD(grada)*dt; aa[aa<0] = 0; a[1,,] = aa;
        ar = a[2,,] + SCALE_GRAD(gradr)*dt; ar[ar<0] = 0; ar[ar>1] = 1; a[2,,] = ar;
        Timers <- TIMER_TOC("M_STEP_UPDATES", Timers);

        if (verbose) setTxtProgressBar(pb, t);
    }
    endTime <- proc.time();
    Timers <- TIMER_TOC("ALL_ITERS", Timers);
    TIMER_SUMMARY(Timers);
    if (verbose) close(pb);

    elapsedTime = endTime[[3]] - startTime[[3]];
    cat(sprintf("Completed %d iterations in %1.1f seconds\n%1.3f seconds / iteration\n%1.3e seconds / TSNK / iteration.\n", t, elapsedTime, elapsedTime/t, elapsedTime/t/T/S/N/K));

    ## Extract a clustering from the probabilities
    clust  =  apply(qnk, 1, "which.max");
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
    results= list(a = a, al = al, v0 = v0, l0 = l0, qnk = qnk, F = F[1:t], L = L, Lclust = Lclust, clust=clust, pclust=pclust, dclust=dclust, numIters = t, seed = seed, exitMode = exitMode, history = history, misc = misc);
}
