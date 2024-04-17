# utility internal functions

# simple function to extract the extension of a file
file_ext <- function(x){
  utils::tail(unlist(strsplit(x,".",fixed = TRUE)),n=1)
}

