#' @export
ApplySum <- function(X, margins, ii = NULL, ir = NULL){
    dims = dim(X);
    nd = length(dims);

    if (length(margins)==0)
        return(list(Y = sum(X), ii = c(), ir=1:nd));
    if (length(margins)==nd)
        return(list(Y = X, ii = 1:nd, ir = c()));

    ny = prod(dims[margins]);

    computeIndices = is.null(ii);
    if (!computeIndices){
        ## Don't compute the indices
        results = .C("ApplySum",
            as.double(X), as.integer(dims), as.integer(nd), as.integer(margins-1), as.integer(length(margins)),
            as.double(vector("double", ny)),
            as.integer(computeIndices), 
            ii, length(ii), ir, length(ir), PACKAGE = "physplit.analysis");
    }else{
         ## We need the indices to be computed.
         nx = prod(dims); ni = prod(dims[margins]); nr = nx/ni;
         ii = as.integer(vector("integer", ni));
         ir = as.integer(vector("integer", nr));         
         results = .C("ApplySum",
             as.double(X), as.integer(dims), as.integer(nd), as.integer(margins-1), as.integer(length(margins)),
             as.double(vector("double", ny)),
             as.integer(computeIndices), 
             ii, ni, ir, ni, PACKAGE = "physplit.analysis");
         ii = results[[8]];
         ir = results[[10]];
     }
    if (length(margins)==1){
        Y = array(results[[6]], dim=c(dims[margins],1));
    }
    else{
        Y = array(results[[6]], dim=dims[margins]);
    }
    list(Y = Y, ii = ii, ir = ir);    
}

