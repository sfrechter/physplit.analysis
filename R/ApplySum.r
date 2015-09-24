#' @export
ApplySum <- function(X, margins){
    dims = dim(X);
    nd = length(dims);

    if (length(margins)==0)
        return(sum(X));
    if (length(margins)==nd)
        return(X);

    ny = prod(dims[margins]);
    results = .C("ApplySum",
        as.double(X), as.integer(dims), as.integer(nd), as.integer(margins-1), as.integer(length(margins)),
        as.double(vector("double", ny)), PACKAGE = "physplit.analysis");
    if (length(margins)==1){
        Y = array(results[[6]], dim=c(dims[margins],1));
    }
    else{
        Y = array(results[[6]], dim=dims[margins]);
    }
}

