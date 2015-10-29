#' A function that uses model selection to detect variability in the
#' baseline firing rate of a cell. The algorithm is to assume Poisson
#' spiking in the baseline period and compare the posterior on a
#' single rate explaining the data to having two rates. This is an
#' approximation for the ideal situation of comparing a single rate to
#' the possibility of having more than one rate, i.e. two,three,four,
#' etc. That is, given the vector of baseline rates r, we're saying
#'
#' relatve prob. of single rate = p(M_1|r)/(p(M_2|r) + p(M_3|r) + ...) ~ p(M_1|r)/p(M_2|r).
#'
#' Hence the procedure is theoretically biased towards the single rate
#' model, and the user can supply a counter bias to this in the bias
#' argument.
#'
#' Comparing the two models requires a prior probability on each:
#'
#' p(M_1|r) / p(M_2|r) = p(r|M_1)p(M_1)/p(r|M_2)p(M_2)
#'
#' We actually compare the logarithm of this
#'
#' log [p(M_1|r)/p(M_2|r)] = log [p(r|M_1)/p(r|M_2)] + log [p(M_1)/p(M_2)]
#' 
#' We assume p(M_1) = p(M_2) so that the default comparison is that of
#' the likelihoods, and the latter term is zero, but the user can
#' provide the logarithm of the prior probability ration in the bias
#' term. The final decision of the model is:
#'
#' M1 if log [p(r|M_1)/p(r|M_2)] + bias > 0, otherwise M2.
#'
#' The likelihoods are computed by marginalizing over the rates. For
#' example, for the two rate case, we have
#'
#' p(r|M_2) = int(l1, lmin, lmax) int(l0, lmin, lmax) p(r|l0,l1,M_2) p(l0)p(l1)
#'
#' The prior on the rates are assumed uniform, i.e. p(l0) = p(l1) = 1/(lmax - lmin).
#'
#' The integration is performed numerically by splitting the range
#' (lmin-lmax) into a fixed pitch grid, computing the probability at
#' each grid point, and summing the result scaled by the grid pitch.
#' 
#' @param r A vector of length N contain the baseline spike counts for N odors.
#' @param lrange The range of baseline rates to integrate over. If lrange is NULL, will compute this from the data.
#' @param nsd The number of standard deviations around the observed range of rates to integrate over.
#' @param nlvals The number of points along each dimension to use in perform the numerical integration.
#' @param bias The bias term in favour of M1.
#' @return A list consisting of
#'
#'   \item{bestModel}{1 or 2, indicating the better model.}
#'
#'   \item{lpM1}{The log likelihood of the single rate model.}
#'
#'   \item{lpM2}{The log likelihood of the two rate model.}
#' @export
DetectVariableBaselineUsingBayesianModelSelection <- function(r, lrange = NULL, nsd = 5, nlvals = 1001, bias = 0){    
    n  = length(r);

    if (is.null(lrange)){
        lmin = min(r);
        lmax = max(r);
        
        if (lmin == lmax)
            lmax = lmin + 1e-6;
        
        ## To determine the range of the integration, we'll pick a range [l0,
        ## l1] around [lmin - lmax] such that the observed values of the min
        ## and the max are three standard deviations above, and below,
        ## respectively.
        
        ll  = lmin;
        if (ll>0)
            ll = uniroot(f=function(l) l + nsd*sqrt(l) - lmin, c(0,lmin))$root;
        uu = 1;
        if (lmax>1)
            tryCatch({
                         uu = uniroot(f = function(l) l - nsd*sqrt(l) - lmax, c(lmax,10*lmax))$root;
                     },error= function(e){
                           uu = lmax + nsd*sqrt(lmax);
                       });
        
    }else{
         ll = min(lrange);
         uu = max(lrange);
     }
    if (ll<1e-12)
        ll = 1e-12;    
    lvals = seq(ll,uu,length.out=nlvals);
    dl = lvals[[2]]-lvals[[1]];

    ## Compute the probability of the 1 rate model.
    ## p(r|M1) = 1/lrange int_l  l^sum(r) exp(-len(r)* l);
    
    lpi = outer(log(lvals),r) - outer(lvals,rep(1,n))
    lp  = apply(lpi, MARGIN=1,FUN="sum");
    lpMax = max(lp);
    lps = lp - lpMax;
    Q = sum(exp(lps))*dl;
    lpM1 = lpMax + log(Q) - log(diff(range(lvals)));
    
    ## Compute the probability of the 2 rate model.
    ## p(r|M2) = 1/2 * 1/lrange^2 int(l1,lmin,lmax) int(l0,lmin,lmax) prod(i=1,N) (0.5 * p(r_i|l0) + 0.5*p(r_i|l1))
    
    ll = expand.grid(lvals,lvals);
    l0 = ll[,1]; l1 = ll[,2];
    nl = nrow(l0);
    lpi  = log(exp(outer(log(l0),r) - outer(l0,rep(1,n))) + exp(outer(log(l1),r) - outer(l1,rep(1,n))));
    lp   = apply(lpi, MARGIN=1,FUN="sum");
    maxLp = max(lp);
    lps = lp - maxLp;

    Q = sum(exp(lps)*(l0<=l1))*dl*dl;

    lpM2 = maxLp + log(Q) - n*log(2)-2*log(diff(range(lvals)))

    return(list(bestModel = 1 + (lpM2>(lpM1+bias)), lpM1 = lpM1, lpM2 = lpM2));
}
