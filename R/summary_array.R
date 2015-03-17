#' Create the raw summary array for all spikes
#'
#' @param x A list of smoothed psth objects for cell/odours combinations.
#'   Defaults to \code{\link{smSpikes}}
#' @export
#' @import physplitdata
#' @examples
#' summary_array=create_raw_summary_array()
#' # histogram of baseline firing for all cells / odours
#' hist(summary_array[,,'baseline'])
#' # collapse all data for different odours for the same cell
#' # i.e. average baseline firing rate for each cell
#' baseline_cell=rowMeans(summary_array[,,'baseline'], na.rm=TRUE)
#' hist(baseline_cell, xlab='Firing freq /Hz')
#'
#' # Plot density distributions by cell group
#' pns=subset(PhySplitDB, Group=='PN' & cell %in% names(smSpikes))$cell
#' physplit=PhySplitDB[match(names(smSpikes), PhySplitDB$cell), ]
#' rownames(physplit)=physplit$cell
#' physplit$baseline=baseline_cell[physplit$cell]
#' library(ggplot2)
#' qplot(baseline,col=Group, data=subset(physplit, Group%in%c("L","O","PN")), geom='density')
create_raw_summary_array<-function(x=physplitdata::smSpikes) {

  allfreqs=lapply(x,function(psthsforcell) sapply(psthsforcell,function(psth) psth$freq))
  allodours=unique(unlist(sapply(x,names)))
  allfreqs_allodours=lapply(allfreqs,addnacols,allodours)

  num_cells <- length(x)
  num_odours <- length(allodours)
  num_stats <- 7
  stat_names <- c("baseline", "max1", "max2", "max3", "max4", "max5", "max6")
  summary_array <- array(dim=c(num_cells, num_odours, num_stats), dimnames=list(names(x), allodours, stat_names))

  summary_array[, , 'baseline'] <- t(sapply(allfreqs_allodours, function(x) colMeans(x[2:9, ])))
  colMax=function(x) apply(x, 2, max)
  summary_array[, , 'max1'] <- t(sapply(allfreqs_allodours, function(x) colMax(x[12:19, ])))
  summary_array[, , 'max2'] <- t(sapply(allfreqs_allodours, function(x) colMax(x[16:23, ])))
  summary_array[, , 'max3'] <- t(sapply(allfreqs_allodours, function(x) colMax(x[20:27, ])))
  summary_array[, , 'max4'] <- t(sapply(allfreqs_allodours, function(x) colMax(x[24:31, ])))
  summary_array[, , 'max5'] <- t(sapply(allfreqs_allodours, function(x) colMax(x[28:35, ])))
  summary_array[, , 'max6'] <- t(sapply(allfreqs_allodours, function(x) colMax(x[32:39, ])))
  summary_array
}
