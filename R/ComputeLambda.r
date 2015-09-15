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

    U  = aperm(array(results[[7]], dim=c(T,S,K,N)),c(1,2,4,3));
    U1 = aperm(array(results[[8]], dim=c(T,S,K,N)),c(1,2,4,3));
    V  = aperm(array(results[[9]], dim=c(T,S,K,N)),c(1,2,4,3));
    L  = aperm(array(results[[10]],dim=c(T,S,K,N)),c(1,2,4,3));

    list(U = U, U1 = U1, V = V, L = L);
}

