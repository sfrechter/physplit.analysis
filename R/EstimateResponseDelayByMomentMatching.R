#' Estimates the response delay by cross correlating the observed response with the odor profile convolved with an exponential filter.
#'
#' A function that detects the delay of a response relative to the
#' odor onset by finding the delay with maximum cross correlation of
#' the given response rate timecourse, with that of the odor onset
#' convolved with a decaying exponential. The parameters of the
#' exponential are found by moment matching.
#'
#' @param y A length T vector of the mean firing rate of the neuron in a series of bins.
#' @param x0 A binary length T vector indicating the odor window.
#' @param verbose Whether to print out the results
#' @param plot Whether to plot the delayed exponential on top of the firing rate.
#' @param minDelay The minimum delay value that can be used.
#' @param minDelay The minimum delay value that can be used.
#' @param fitOdorProfileOnError Whether to fit the odor profile itself if learning the filter parameters failed.
#' @return A list consisting of
#'
#'   \item{a0}{The amplitude of the exponential.}
#'
#'   \item{a1}{The autoregressive coefficient of the exponential}
#'
#'   \item{tau}{The time constant of the exponential. a1 = exp(-1/tau).}
#'
#'   \item{l0}{The baseline firing rate.}
#'
#'   \item{delay}{The delay with the highest correlation.}
#'
#'   \item{status}{A string. "OK" if the parameters could be learned, "ERROR" if not.}
#' @export
EstimateResponseDelayByMomentMatching <- function(y, x0, verbose = FALSE, plot=FALSE, minDelay = 0, fitOdorProfileOnError = TRUE){
    t0 = min(which(x0>0));
    l0 = mean(y[1:(t0-1)]);
    Tp = sum(x0);

    z  = (y - l0)*((y-l0)>=0);
    M  = sum(z);
    M2 = sum(z^2);

    b = M/Tp;
    ## Try to figure out tau by solving M2/b^2 - Tp  = x - x exp(-Tp/x)
    res = tryCatch({
                       uniroot(function(x) M2/b^2 - Tp - x * exp(-Tp/x) + x, c(0,100))
                   }, error = function(e) {
                          return(NULL);
                      }
                   );
    ## If there's too much noise this sometimes doesn't work.
    ## so just assume 0 delay.
    a1 = NULL;
    a0 = NULL;
    tau = NULL;
    delay = 0;

    if (is.null(res)){
        status = "ERROR";
        if (fitOdorProfileOnError){
            cc = ccf(y,x0,plot=FALSE,lag.max=length(y));
            delay = max(cc$lag[which.max(cc$acf)], minDelay);
        }
    }else{
         status = "OK";
         tau = res$root;
         a1 = exp(-1/tau);
         a0 = b*(1 - a1);
         f = stats::filter(x0, a1, method="recursive")*a0;
         cc = ccf(y,f,plot=FALSE,lag.max=length(y));
         delay = max(cc$lag[which.max(cc$acf)], minDelay);
     }

    if (verbose){
        cat(sprintf("b : %1.3f\n", b));
        cat(sprintf("a0: %1.3f\n", a0));
        cat(sprintf("a1: %1.3f\n", a1));
    }

    if (plot){
        pdf(file="estimateResponseOnsetByMomentMatching.pdf", width=8,height=8);
        plot(y);
        lines(binhf::shift(f,places=delay),col="red");
        dev.off();
    }
    return(list(a1=a1,a0=a0,tau=tau,delay = delay, status=status, l0 = l0));
}
