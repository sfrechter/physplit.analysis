#' @export
BroadcastMultiply <- function(X, y, margins, ii = NULL, ir = NULL){
    dx = dim(X);
    dy = drop(dim(y));

    if (prod(dx[margins])!=prod(dy))
        stop("Number of elements in margin dims does not equal number of elements of y.");

    ndx = length(dx);
    nx  =   prod(dx);
    
    if (length(margins)==0)
        return(list(Xy = X*y, ii = 0,  ir=1:nx));

    if (length(margins)==ndx)
        return(list(Xy = X*y, ii = 1:nx, ir=0));
   
    computeIndices = is.null(ii);
    if (!computeIndices){
        ## Don't compute the indices
        results = .C("BroadcastMultiply",
            as.double(X), as.double(y), as.integer(dx), as.integer(ndx), as.integer(margins-1), as.integer(length(margins)),
            as.double(vector("double", nx)),
            as.integer(computeIndices), 
            ii, length(ii), ir, length(ir), PACKAGE = "physplit.analysis");
    }else{
         ## We need the indices to be computed.
         ni = prod(dx[margins]); nr = nx/ni;
         ii = as.integer(vector("integer", ni));
         ir = as.integer(vector("integer", nr));         
         results = .C("BroadcastMultiply",
             as.double(X), as.double(y), as.integer(dx), as.integer(ndx), as.integer(margins-1), as.integer(length(margins)),
             as.double(vector("double", nx)),
             as.integer(computeIndices), 
             ii, ni, ir, nr, PACKAGE = "physplit.analysis");
         ii = results[[9]];
         ir = results[[11]];
     }
    Xy = array(results[[7]], dim=dx);
    list(Xy = Xy, ii = ii, ir = ir);    
}

