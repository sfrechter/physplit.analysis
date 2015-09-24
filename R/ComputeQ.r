ComputeQ <- function(LL){
    ## L is in tsnk format.
    T = dim(LL)[[1]];
    S = dim(LL)[[2]];
    N = dim(LL)[[3]];
    K = dim(LL)[[4]];
    d = c(T,S,N,K);
    results = .C("ComputeQ",
        LL,as.integer(d),
        as.double(vector("double", N*K)),
        as.double(vector("double", prod(d))),
        PACKAGE = "physplit.analysis");
    qnk   = array(results[[3]], dim=c(N,K));
    qtsnk = array(results[[4]], dim=d);

    list(qnk = qnk, qtsnk = qtsnk);
}

