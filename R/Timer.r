TIMER_INIT <- function(status="ON"){
    return(list(status=status));
}
TIMER_TIC <- function(whichTimer, T){
    if (T$status=="OFF")
        return(T);

    if (is.null(T[[whichTimer]])){
        T[[whichTimer]] = list(startTime=proc.time(), stopTime=Inf, dt = list());
    }else{
         T[[whichTimer]]$startTime = proc.time();
         T[[whichTimer]]$stopTime = Inf;
     }

    return(T);
}

TIMER_TOC <- function(whichTimer, T){
    if (T$status=="OFF")
        return(T);
    if (is.null(T[[whichTimer]]))
        stop(sprintf("No timer called '%s' found.", whichTimer));
    T[[whichTimer]]$stopTime = proc.time();
    T[[whichTimer]]$dt[[length(T[[whichTimer]]$dt)+1]] = (T[[whichTimer]]$stopTime - T[[whichTimer]]$startTime)[["elapsed"]];
    return(T);
}

TIMER_SUMMARY <- function(T){
    if (T$status == "OFF")
        return(NULL);
    fields = names(T);
    cat("Timer summary:\n");
        
    for (f in fields){
        if (f=="status")
            next;
        tt = unlist(T[[f]]$dt);
        cat(sprintf("\t%20s: %8d calls, total time %8.3f, %8.3f +/- %8.3f / call\n",f,length(tt),sum(tt),mean(tt), sd(tt)));
    }
}
