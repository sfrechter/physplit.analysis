#' @export
ComputeLambda <- function(X,a,al,v0,l0){
    ## X is in tsn format.
    ## a is in 2sk format.
    ## U,L are in tskn format.
    T = dim(X)[[1]];
    S = dim(X)[[2]];
    N = dim(X)[[3]];
    K = dim(a)[[3]];    
    d = c(dim(X),K);

    ## If the parameters are specified as 1D matrices, then we assume
    ## one parameter per cell, otherwise one per cell AND
    ## cluster. This latter mode is useful when pre fitting parameters
    ## to all clusters.
    
    if ((length(v0) != N) && any(dim(v0) != c(N,K))) stop("Expected len(v0) = N or dim(v0) = (N,K).");
    if ((length(l0) != N) && any(dim(l0) != c(N,K))) stop("Expected len(l0) = N or dim(l0) = (N,K).");
    if ((length(al) != N) && any(dim(al) != c(N,K))) stop("Expected len(al) = N or dim(al) = (N,K).");

    clusterSpecificParams = (length(al) != N);

    results = .C("ComputeLambda",
        X,a,as.integer(d),
        al, v0, l0,
        as.integer(clusterSpecificParams),
        as.double(vector("double", prod(d))),
        as.double(vector("double", prod(d))),
        as.double(vector("double", prod(d))),
        as.double(vector("double", prod(d))),
        PACKAGE = "physplit.analysis");

    U  = array(results[[8]],  dim=d);
    U1 = array(results[[9]],  dim=d);
    V  = array(results[[10]],  dim=d);
    L  = array(results[[11]], dim=d);

    list(U = U, U1 = U1, V = V, L = L);
}

