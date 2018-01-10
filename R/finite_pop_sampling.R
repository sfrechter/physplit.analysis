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
#'   NB this is the hypergeometric distribution.
#'
#' @param n sample size
#' @param N Population size
#' @param p positive probability (e.g. of being an LHN)
#' @param npositive Number of true positives (e.g. actual number LHNs in a
#'   tract)
#' @param replicates Number of samples to draw
#' @return A vector of length \code{replicates} with counts of sample positives.
#' @export
#' @importFrom stats rhyper
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
#' @family population-sampling
sample_finite_population <- function(n, N, p=NULL, npositive=NULL, replicates=1) {
  if(is.null(npositive))
    npositive=round(p*N)
  # imagine we have integers 1:N (the population size)
  # 1:npositive (are the true positives)
  # Now we draw a sample of size n from the integers 1:N
  # the number of sample positives is the number of sampled values in the range
  # 1:npositive
  rhyper(nn=replicates, m=npositive, n=N-npositive, k=n)
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
#'   generated the observed number of sample positives. It has class
#'   \code{truepos}.
#' @export
#'
#' @examples
#' # Imagine we have sampled 10 profiles from a tract of 48 and found 2 LHNs
#' tps=truepos_given_sample(samplepos = 2, n=10, N=48)
#'
#' hist(tps, breaks=0:49-.5, col='red')
#' plot(ecdf(tps))
#'
#' # the mode should be the Maximum Likelihood Estimate
#' # (if enough replicates were used)
#' summary(tps)
#' # 95% confidence interval
#' summary(tps, alpha=.05)
#'
#' @importFrom tidyr gather
#' @family population-sampling
truepos_given_sample <- function(samplepos, n, N, replicates=1000) {
  m=mapply(rhyper, m=0:N, n=N:0, k=n, nn=replicates)
  dm=as.data.frame(m)
  colnames(dm)=1:ncol(dm)
  gm=gather(dm, key = 'truepos', value = 'n')
  gm$truepos=as.integer(gm$truepos)
  res=gm$truepos[gm$n==samplepos]
  class(res)='truepos'
  res
}

#' @export
#' @rdname truepos_given_sample
#' @param object Sample counts to summarise
#' @param alpha The confidence interval is (1-alpha)*100\% (i.e. alpha=0.1 => 90\% CI)
#' @importFrom stats quantile
summary.truepos <- function(object, alpha=0.1, ...) {
  quantile.levels=c(alpha/2, 0.5, 1-alpha/2)
  qs=quantile(object, quantile.levels)
  res=c(qs[1], Median=unname(qs[2]), Mode=pmode(object), qs[3])
  class(res)=c("summaryDefault", "table")
  res
}

pmode <- function(x) as.integer(names(which.max(table(x))))

#' Approximate (1-alpha)100\% confidence interval for proportion of a population
#'
#' @details Note that this value is generally similar to that obtained with the
#'   sampling based approach in \code{\link{truepos_given_sample}} when
#'   \code{pest=0.5} but becomes an increasingly bad approximation as the
#'   proportion tends to 0 or 1.
#'
#' @param n Sample size
#' @param pest Estimated positive proportion
#' @param N Population size (default implies infinite population)
#' @param alpha The confidence interval is (1-alpha)*100\%
#' @family population-sampling
#' @references See https://onlinecourses.science.psu.edu/stat414/node/264
#' @export
#' @examples
#' # 95% confidence interval for population size 43, sample size 10 and estimated
#' # proportion of 0.1 ...
#' # expressed as a proportion
#' prop.ci(10, pest=0.1, N=43)
#' # as a number of positives
#' prop.ci(10, pest=0.1, N=43)*43
#'
#' ## Compare with sampling based calculation
#' prop.ci(10, pest=0.5, N=48, alpha=.1)*48
#' summary(truepos_given_sample(samplepos = 5, n=10, N=48))
#' # more different
#' prop.ci(10, pest=0.2, N=48, alpha=.1)*48
#' summary(truepos_given_sample(samplepos = 2, n=10, N=48))
#' @importFrom stats qnorm
prop.ci <- function(n, pest=0.5, N=Inf, alpha=0.05) {
  # normal quantile for given alpha
  z.a=qnorm(p=1-alpha/2)

  # finite population correction
  correction=if(is.finite(N)) (N-n)/(N-1) else 1

  half.interval=z.a * sqrt(pest*(1-pest)/n * correction)
  c(pest-half.interval, pest+half.interval)
}

#' Estimate sample size to find population proportion with given tolerance
#'
#' @details Estimated proportion should have (1-alpha)*100\% confidence of being
#'   within p +/- epsilon of the true proportions
#'
#' @param epsilon Tolerange of proportion estimate (see details)
#' @param pest The estimated proportion (leave at 0.5 if unknown)
#' @param N Population size (default implies infinite population)
#' @inheritParams prop.ci
#' @references See https://onlinecourses.science.psu.edu/stat414/node/264
#' @return numeric
#' @family population-sampling
#' @export
#'
#' @examples
#' required.sample.size(0.04, alpha=0.1)
#' required.sample.size(0.04, alpha=0.1, N=40)
#' # epsilon here is equivalent to +/- 4
#' required.sample.size(10/40, alpha=0.1, N=40)
required.sample.size <- function(epsilon, pest=0.5, alpha=0.05, N=Inf) {
  z.a=qnorm(p=1-alpha/2)
  m=z.a^2*pest*(1-pest)/epsilon^2
  if(is.finite(N)) {
    # finite population correction
    m / (1 + (m-1) / N)
  } else {
    m
  }
}
