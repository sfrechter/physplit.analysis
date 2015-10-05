#' @export
HISTORY <- function(history, keepHistory, iter, numIters, at=NULL, ...){
    if (keepHistory){
        if (is.null(at))
            at = 1:numIters;
        fields = list(...);
        values = lapply(fields, FUN=function(x) eval.parent(as.symbol(x),n=3));

        ihist = which(at==iter);
        if (length(ihist)>0){
            for (ih in ihist){
                for (i in 1:length(fields)){
                    fld = fields[[i]];
                    d = dim(values[[i]]);
                    if (is.null(history[[fld]]))
                        history[[fld]] = drop(array(dim=c(d,length(at))));
                    history[[fld]][prod(d)*(ih-1) + (1:prod(d))] = values[[i]];
                }
            }
        }
    }
    history;
}

