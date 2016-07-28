#' Convert any of Shahar's ids (or a filename) to a short id
#'
#' @details shortids look like 120726c3
#' @param x A stack, cell or short id or a filename that begins with one of these.
#'
#' @return A vector of shortids, named by the original id when there were more than 1
#' @export
#'
#' @examples
#' anyid2shortid("120726BJK1742SF274LC.lsm")
anyid2shortid <- function(x) {
  x=basename(x)
  if(length(x)>1) return(sapply(x, anyid2shortid))
  # catch anything that already looks like a short id
  if(grepl("^1[0-9]{5}c[0-9]$", x)) return(x)

  if(isTRUE(substr(x,1,4) == "nm20")) return(substr(x, 5, 12))
  stack=substr(x,1,7)
  x=PhySplitDB$cell[match(stack, physplitdata::PhySplitDB$stack)]
  return(substr(x, 5, 12))
}
