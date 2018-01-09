#' Sample from finite population with known number of true positives
#'
#' @details Imagine we have integers \code{1:N} (the population size). We say
#'   that the first \code{npositive} integers are the true positives in the
#'   population (npositive may be 0).
#'
#'   Now we draw a sample of size \code{n} from the integers \code{1:N}. For
#'   this sample we say that the observed positives are those integers
#'   \code{<=npositive}.
#'
#' @param n sample size
#' @param N Population size
#' @param p positive probability (e.g. of being an LHN)
#' @param npositive Number of true positives (e.g. actual number LHNs in a
#'   tract)
#' @param replicates Number of samples to draw
#' @return A vector of length \code{replicates} with counts of sample positives.
#' @export
#' @examples
#' # Draw a random sample of size 24 tracings for a population of 96 profiles in
#' # a tract known to have 36 LHNs
#' sample_finite_population(24, N=96, npositive=36)
#' # Draw 1000 random samples, plot the distribution of observed sample positives
#' rand.samples=sample_finite_population(24, N=96, npositive=36, replicates=1000)
#' mean(rand.samples)
#' quants=quantile(rand.samples, c(0.05,0.95))
#' \donttest{
#' library(ggplot2)
#' qplot(rand.samples, binwidth=1, xlab='Observed Sample Positives') +
#'   geom_vline(xintercept = quants, colour='red')
#' }
#'
#' \donttest{
#' # Compare with binomial distribution
#' resdf=data.frame(x=popsample(10,50,p=.5, replicates = 100000), type='popsample')
#' resdf=rbind(resdf, data.frame(x=rbinom(100000, size=10,p=.5), type='rbinom'))
#' library(ggplot2)
#' qplot(x, col=type, data=resdf, geom='density')
#' }
#' @seealso \code{\link{truepos_given_sample}} for estimating the true number of
#'   positives in a finite population given a sample.
#'
sample_finite_population <- function(n, N, p=NULL, npositive=NULL, replicates=1) {
  if(is.null(npositive))
    npositive=round(p*N)
  # imagine we have integers 1:N (the population size)
  # 1:npositive (are the true positives)
  # Now we draw a sample of size n from the integers 1:N
  # the number of sample positives is the number of sampled values in the range
  # 1:npositive

  # This looks a bit like a binomial distribution but in fact has lower
  # variance since the extreme values are slightly less likley to occur
  replicate(sum(sample.int(n=N, size = n)<=npositive), n = replicates)
}

#' Estimate distribution of true positives given sampling resuts
#'
#' @details The idea is to generate random realisations for all possible numbers
#'   of true positives, choose only those cases that resulted in the observed
#'   number of sample positives, and then use that empirical distribution of
#'   simulated true positives to estimate the unknown true positive.
#'
#'   We can use this approach to estimate
#'
#' @param samplepos Number of positives observed in sample
#' @param n Sample size
#' @param N Population size
#' @param replicates Number of replicates per tested true pos number
#'
#' @return a vector containing population true positive counts that could have
#'   generated the observed number of sample positives.
#' @export
#'
#' @examples
#' # Imagine we have sampled 10 profiles from a tract of 48 and found 2 LHNs
#' tps=truepos_given_sample(samplepos = 2, n=10, N=48)
#'
#' hist(tps, breaks=0:49-.5, col='red')
#' plot(ecdf(tps))
#' # calculate median and 90% range
#' quantile(tps, c(.05,.5,.95))
#' # calculate mode i.e. most likely true positive number
#' pmode <- function(x) as.integer(names(which.max(table(x))))
#' pmode(tps)
#' @importFrom tidyr gather
#' @seealso \code{\link{sample_finite_population}}
truepos_given_sample <- function(samplepos, n, N, replicates=1000) {
  m=mapply(sample_finite_population, n=n, N=N, npositive=0:N, replicates=replicates)
  dm=as.data.frame(m)
  colnames(dm)=1:ncol(dm)
  gm=gather(dm, key = 'truepos', value = 'n')
  gm$truepos=as.integer(gm$truepos)
  gm$truepos[gm$n==samplepos]
}
