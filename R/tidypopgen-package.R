#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @import dplyr
#' @import tibble
#' @importFrom Rcpp sourceCpp
#' @importFrom rlang :=
#' @useDynLib tidypopgen, .registration = TRUE
#' @importFrom ggplot2 autoplot
## usethis namespace: end
NULL

#' @importFrom generics tidy
#' @export
generics::tidy

#' @importFrom generics augment
#' @export
generics::augment

#' @importFrom ggplot2 autoplot
#' @export
ggplot2::autoplot

# Unexported functions from bigsnpr
getCode <- utils::getFromNamespace("getCode", "bigsnpr")
getInverseCode <- utils::getFromNamespace("getInverseCode", "bigsnpr")
readbina <- utils::getFromNamespace("readbina", "bigsnpr")
write.table2 <- utils::getFromNamespace("write.table2", "bigsnpr")
writebina <- utils::getFromNamespace("writebina", "bigsnpr")
