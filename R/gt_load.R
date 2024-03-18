#' Load a gen_tibble
#'
#' Load a `gen_tibble` previously saved with [gt_save()].
#'
#' @param file_name the file name, including the full path. If it does not end
#' with '.gt', the extension will be added.
#' @returns a [`gen_tibble`]
#' @export

gt_load <-function(file_name = NULL) {
  if (file_ext(file_name)!="gt"){
    file_name <- paste0(file_name,".gt")
  }
  if (!file.exists(file_name)){
    stop(file_name, " does not exist")
  }
  x <- readRDS(file_name)
  backingfiles <- gt_get_file_names(x)
  if (!file.exists(backingfiles[1])){
    stop(backingfiles[1], " does not exist")
  }
  if (!file.exists(backingfiles[2])){
    stop(backingfiles[2], " does not exist")
  }
  return(x)
}

gt_get_bigsnp<-function(x){
  attr(x$genotypes,"bigsnp")
}
