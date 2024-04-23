
#' @export
dimnames.FBM <- function(x){
  attr(x,"dim_names")
}

#' @export
"dimnames<-.FBM" <- function(x, value){
  if (any(!is.list(value),length(value)!=2)){
    stop("value for dimnames should be a list of two elements: a vector of row names and a vector of column names")
  }
  dimnames_dim <- unlist(lapply(value, length))
  if (any(dimnames_dim[1]!=dim(x)[1], dimnames_dim[2]!=dim(x)[2])){
    stop("row and column names should be the same length as their respective dimensions in the FBM")
  }
  attr(x,"dim_names")<-lapply(value, as.character)
  return(x)
}
