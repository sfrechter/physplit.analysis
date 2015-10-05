ComputeGtsnk <- function(X,a){
    ## X is in tsn format.
    ## a is in 2sk format.

    T = dim(X)[[1]];
    S = dim(X)[[2]];
    N = dim(X)[[3]];
    K = dim(a)[[3]];
    d = c(dim(X),K);
    results = .C("ComputeGtsnk",
        X,a,as.integer(d),
        as.double(vector("double", prod(d))), PACKAGE = "physplit.analysis");

    G = array(results[[4]], dim=c(T,S,N,K));
}

