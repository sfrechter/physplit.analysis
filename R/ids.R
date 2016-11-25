#' Convert any identifier / file path to Shahar's short id or stack id
#'
#' @description \code{anyid2shortid} converts any of Shahar's ids (or a
#'   filename) to a short id
#'
#' @details shortids look like 120726c3
#' @param x A stack, cell or short id or a filename that begins with one of
#'   these.
#'
#' @return A vector of short ids or stack ids, named by the original id when
#'   there were more than 1
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

#' \code{anyid2stack} converts any of Shahar's ids (or a filename) to a stack id
#' @rdname anyid2shortid
#' @export
#' @examples
#' anyid2stack("120726BJK1742SF274LC.lsm")
#' anyid2stack("120726c3")
anyid2stack <- function(x) {
  m=match(anyid2shortid(x), physplitdata::PhySplitDB$shortid)
  physplitdata::PhySplitDB$stack[m]
}
