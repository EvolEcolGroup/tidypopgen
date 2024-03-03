
# Look into eval_select, which can work with a vector with names
sub_sel <- function(.data, FUN) {
  # defuse FUN
  defused_FUN <- enquo(FUN)
  tidyselect::eval_select(expr=defused_FUN, data = .data)
}

my_data <- 1:5
names(my_data) <- c("rs123","rs143","rs67","ps4343","rs125")
sub_sel (my_data, ends_with("3"))
sub_sel (my_data, c(3,1,4))


## we could parse the defused FUN and eval_select only if it contains a selection
## helper
## or better select_loci_if
## and select_loci
## we can then add a flip_select_loci_if
## and flip_select_loci

