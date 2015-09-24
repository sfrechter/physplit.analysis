ComputeLambda <- function(X,a,al,v0,l0){
    ## X is in tsn format.
    ## a is in 2sk format.
    ## U,L are in tskn format.
    T = dim(X)[[1]];
    S = dim(X)[[2]];
    N = dim(X)[[3]];
    K = dim(a)[[3]];
    d = c(dim(X),K);
    results = .C("ComputeLambda",
        X,a,as.integer(d),
        al, v0, l0,
        as.double(vector("double", prod(d))),
        as.double(vector("double", prod(d))),
        as.double(vector("double", prod(d))),
        as.double(vector("double", prod(d))),
        PACKAGE = "physplit.analysis");

    U  = array(results[[7]],  dim=d);
    U1 = array(results[[8]],  dim=d);
    V  = array(results[[9]],  dim=d);
    L  = array(results[[10]], dim=d);

    list(U = U, U1 = U1, V = V, L = L);
}

