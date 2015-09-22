HISTORY <- function(history, keepHistory, iter, numIters, ...){
    if (keepHistory){
        fields = list(...);
        values = lapply(fields, FUN=function(x) eval.parent(as.symbol(x),n=3));

        for (i in 1:length(fields)){
            fld = fields[[i]];
            d = dim(values[[i]]);
            if (is.null(history[[fld]]))
                history[[fld]] = drop(array(dim=c(d,numIters)));
            history[[fld]][prod(d)*(iter-1) + (1:prod(d))] = values[[i]];
        }
    }
    history;
}

