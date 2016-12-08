#' @export
#' @importFrom stats rexp rnorm rgamma lm
FitCellSpecificParameters <- function(X, Y, a, l0=NULL, v0 = NULL, al=NULL, kalpha=10, thalpha=3/20, sdv0 = 0.1, taua=1, taul0=0.5, minIters = 0, numIters=10000, dt=1e-5, seed = 0, verbose=TRUE, timer="OFF", slopeRatioToStop=500, numSlopePoints=20, checkToStopEvery=100, keepHistory = NULL, keepHistoryAt = NULL){
    Timers = TIMER_INIT(status=timer);

    ## If X and Y have data for just 1 cell, add a singleton dimension so the rest of the code doesn't complain.
    if (length(dim(X))==2) X = array(X,dim = c(dim(X),1));
    if (length(dim(Y))==2) Y = array(Y,dim = c(dim(Y),1));

    T = dim(X)[[1]]; # Number of bins
    S = dim(X)[[2]]; # Number of odors
    N = dim(X)[[3]]; # Number of cells
    K = dim(a)[[3]]; # Number of clusters

    SCALE_GRAD <- function(g){
        ng = sqrt(sum(g^2));
        return(g/max(1,ng));
    }

    Y = array(rep(Y,K), dim=c(dim(Y), K));

    set.seed(seed);

    if (is.null(l0))
        l0     = 0*array(rexp(N*K, taul0),    dim=c(N,K));

    if (is.null(v0))
        v0     = 0*array(rnorm(N*K, sd=sdv0), dim=c(N,K));

    if (is.null(al))
        al     = 1 + 0*array(rgamma(N*K,shape=kalpha,rate=thalpha), dim=c(N,K));

    F      = array(NaN, dim = c(numIters,N,K));
    slopeRatio = array(NaN, dim = c(N,K));

    ## Initialize the history object
    history = NULL;

    ## Indices for marginalizing out dims
    ii34 = NULL; ir34 = NULL;

    ## Begin the iterations
    if (verbose) pb = txtProgressBar(min=1,max=numIters,initial=1,style=3);
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

        ## Compute the objective function
        Timers <- TIMER_TIC("E_STEP_SCALARS", Timers);
        if (is.null(ii34)){
            res = ApplySum(LL, c(3,4), ii=ii34, ir=ir34);
            ii34 = res$ii; ir34 = res$ir;
        }
        Eqll = ApplySum(LL, c(3,4), ii=ii34, ir=ir34)$Y;
        lprAl = (kalpha-1)*log(al) - al/thalpha; ## Prior on alpha
        F[t,,] = Eqll + lprAl;
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
                 for (n in 1:N){
                     for (k in 1:K){
                         y0 = F[x,n,k];
                         mstart = lm(y0 ~ x);
                         y1 = F[t + ((-numSlopePoints+1):0),n,k];
                         mend = lm(y1 ~ x);
                         ratio = abs(mstart$coefficients[2]/mend$coefficients[2]);
                         slopeRatio[n,k] = ratio;
                     }
                 }
                 if (min(slopeRatio[!is.nan(slopeRatio)]) > slopeRatioToStop){ ## We'll get NaN's if the slope is zero in both cases.
                     exitMode= "SLOPE_RATIO";
                     if (verbose)
                         cat(sprintf("\nBreaking at t = %d due to min SLOPE_RATIO %8.3f > %8.3f\n", t, min(slopeRatio), slopeRatioToStop));
                     break;
                 }
             }
        }


        ## M-STEP: Update the parameters
        ## 1. Compute the gradients
        Timers <- TIMER_TIC("M_STEP_Z", Timers);
        Z   = Y/L - 1;
        ZV  = Z*(V>0);
        Timers <- TIMER_TOC("M_STEP_Z", Timers);

        Timers <- TIMER_TIC("M_STEP_GRADS", Timers);
        gradL0 =  ApplySum(Z,    c(3,4), ii=ii34, ir=ir34)$Y;
        gradAl =  ApplySum(ZV*V, c(3,4), ii=ii34, ir=ir34)$Y + ((kalpha - 1)/al - 1/thalpha); ## Add a prior on Al.
        gradV0 = -ApplySum(ZV,   c(3,4), ii=ii34, ir=ir34)$Y*al;
        Timers <- TIMER_TOC("M_STEP_GRADS", Timers);

        ## 2. Update the parameters
        Timers <- TIMER_TIC("M_STEP_UPDATES", Timers);
        l0 = l0 + SCALE_GRAD(gradL0)*dt; l0[l0<1e-6] = 1e-6;
        al = al + SCALE_GRAD(gradAl)*dt; al[al<0] = 0;
        v0 = v0 + SCALE_GRAD(gradV0)*dt;
        Timers <- TIMER_TOC("M_STEP_UPDATES", Timers);

        if (verbose) setTxtProgressBar(pb, t);
    }
    endTime <- proc.time();
    Timers <- TIMER_TOC("ALL_ITERS", Timers);
    TIMER_SUMMARY(Timers);
    if (verbose) close(pb);

    elapsedTime = endTime[[3]] - startTime[[3]];
    cat(sprintf("Completed %d iterations in %1.1f seconds\n%1.3f seconds / iteration\n%1.3e seconds / TSNK / iteration.\n", t, elapsedTime, elapsedTime/t, elapsedTime/t/T/S/N/K));

    ibest   = drop(apply(F[t,,],MARGIN=1,FUN="which.max"));
    Fbest   = diag(F[t,,ibest]);
    l0.best = diag(l0[,ibest]);
    v0.best = diag(v0[,ibest]);
    al.best = diag(al[,ibest]);

    results= list(a = a, al = al, v0 = v0, l0 = l0, al.best = al.best, l0.best = l0.best, v0.best = v0.best, ibest = ibest, Fbest = Fbest, F = F[1:t,,], slopeRatio = slopeRatio, L = L, numIters = t, seed = seed, exitMode = exitMode, history = history);
}
