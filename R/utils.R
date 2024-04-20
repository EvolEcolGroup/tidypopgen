# utility internal functions

# simple function to extract the extension of a file
file_ext <- function(x){
  utils::tail(unlist(strsplit(x,".",fixed = TRUE)),n=1)
}

# stop if not diploid
stopifnot_diploid <- function(x){
  if (attr(x,"ploidy")!=2){
    stop("this function only works on diploid data")
  }
}


