#' Reorder odour response matrix adding nas as necessary
#'
#' given an odour response matrix with columns for different odours and rows for
#' different sweep timepoints, reorder columns in desired order adding NAs for
#' missing odours
#'
#' @param mat odour response matrix where odours are coloumns and cells are rows
#' @param odours character vector of odour names
#' @return gives back a tidyied up matrix with coloumns aorderd by
#' \code{odours} adding \code{NA}s for missing odour response data
#' @export
addnacols<-function(mat,odours){
  odours_we_have=colnames(mat)
  missing_odours=setdiff(odours,odours_we_have)
  # make an empty matrix to put our data
  finalmat=matrix(ncol=length(odours),nrow=nrow(mat))
  colnames(finalmat)=odours
  # assign the data that we do have to the final matrix
  finalmat[,odours_we_have]=mat
  finalmat
}
