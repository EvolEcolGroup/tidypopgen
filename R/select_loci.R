#' The `select` verb for `loci`
#'
#' An equivalent to [dplyr::select()] that works on the `genotype` column
#' of a `gen_tibble`, using the mini-grammar available for `tidyselect`. The
#' `select`-like evaluation only has access to the names of the loci (i.e. it
#' can select only based on names, not summary statistics of those loci; look
#' at [select_loci_if()] for that feature.
#'
#' Note that the `select_loci` verb does not modify the backing FBM files,
#' but rather it subsets the list of loci to be used stored in the `gen_tibble`.
#'
#' For large datasets (millions of loci), select loci with [dplyr::all_of()]
#' or [dplyr::any_of()], with helpers that match names (such as
#' [dplyr::starts_with()] or [dplyr::everything()]), with numeric positions
#' written as numbers (e.g. `1:100` or `-c(1, 5)`), with a function call
#' that gives positions or names (e.g. `which(...)`, `sample(...)` or
#' `sort(x)`), or with a variable that holds positions or names (e.g.
#' `select_loci(x, idx)`), each optionally negated with `-` or `!`. These
#' selections are resolved directly, which is much faster and uses much less
#' memory than the general `tidyselect` machinery, which has to create an R
#' symbol for every locus name (and R never frees symbols). As in
#' `tidyselect`, a locus with the same name as the variable is selected
#' instead of the value of the variable, and selecting with a variable gives
#' the deprecation warning of `tidyselect`: prefer `all_of(idx)`, which can
#' not be confused with a locus named `idx`. Other expressions (e.g.
#' combinations such as `c(all_of(x), starts_with("rs"))`, or bare loci
#' names) still work, but are slow and memory hungry on very large datasets.
#' To always use `tidyselect` (e.g. if a selection is not resolved as
#' `tidyselect` would resolve it), set `.engine = "tidyselect"`.
#' @param .data a `gen_tibble`
#' @param .sel_arg one unquoted expression, using the mini-grammar of
#'  [dplyr::select()] to select loci. Variable names
#' can be used as if they were positions in the data frame, so expressions
#' like x:y can be used to select a range of variables.
#' @param .engine how the selection is resolved: "auto" (the default)
#'   resolves the common selections described in Details directly, and leaves
#'   the others to `tidyselect`; "tidyselect" always uses `tidyselect`, which
#'   is slow and memory hungry with millions of loci.
#' @returns a `gen_tibble` with a subset of the loci.
#' @seealso [dplyr::select()]
#' @export
#' @examples
#' example_gt <- load_example_gt("gen_tbl")
#'
#' # Select loci by name
#' example_gt_subset <- example_gt %>%
#'   select_loci(all_of(c("rs1", "rs2", "rs3")))
#' show_loci(example_gt_subset)
#'
#' # Select loci by index
#' example_gt_subset <- example_gt %>% select_loci(all_of(c(4, 2, 1)))
#' show_loci(example_gt_subset)
#'
#' # Always use tidyselect to resolve the selection
#' example_gt_subset <- example_gt %>%
#'   select_loci(all_of(c(4, 2, 1)), .engine = "tidyselect")
#'
select_loci <- function(.data, .sel_arg, .engine = c("auto", "tidyselect")) {
  .engine <- match.arg(.engine)
  # defuse the selection criteria
  sel_defused <- rlang::enquo(.sel_arg)

  loci_sel <- NULL
  if (.engine == "auto") {
    # evaluate an ordinary function call (e.g. which() or sample()), as
    # tidyselect would, and select with its value
    sel_defused <- select_loci_eval_call(sel_defused)
    # resolve the most common selections (all_of(), any_of(), numeric
    # positions, values) directly, without building a named vector of all
    # loci
    loci_sel <- select_loci_fast(sel_defused, .data)
  }
  if (is.null(loci_sel)) {
    # general case: create a named vector of all loci
    loci_names_vec <- stats::setNames(
      seq_len(count_loci(.data)),
      nm = loci_names(.data)
    )
    # now evaluate the selection criteria, giving it the loci names to work
    # with
    loci_sel <- tidyselect::eval_select(
      expr = sel_defused,
      data = loci_names_vec
    )
  }

  # modify the table (vctrs::vec_slice() gives the same result as `[`, but
  # with less overhead)
  attr(.data$genotypes, "loci") <-
    vctrs::vec_slice(attr(.data$genotypes, "loci"), unname(loci_sel))
  # prioritise "gen_tbl" class over "sf"
  obj_class <- class(.data)
  if ("sf" %in% obj_class) {
    obj_class <-
      c("gen_tbl", "sf", obj_class[!obj_class %in% c("gen_tbl", "sf")])
    class(.data) <- obj_class
  }
  .data
}

#' Resolve simple loci selections without tidyselect
#'
#' `tidyselect::eval_select()` needs a named vector of all loci and, for any
#' expression that is a call (including `all_of(x)`), it builds a data mask
#' binding every locus name to a symbol. With millions of loci this takes
#' seconds and gigabytes of memory, and the symbols are never released by R.
#' This function resolves the most common selections directly, with the same
#' semantics as `tidyselect`:
#' - `all_of(x)` and `any_of(x)` (also as `tidyselect::` or `dplyr::` calls),
#'   optionally negated with `-` or `!`, where `x` (evaluated in the caller
#'   environment, as `tidyselect` does) is a numeric vector of positions, a
#'   character vector (or factor) of loci names, or `NULL`;
#' - selection helpers that only use the names of the loci (`starts_with()`,
#'   `ends_with()`, `contains()`, `matches()`, `num_range()`, `everything()`
#'   and `last_col()`), optionally negated, which are evaluated with
#'   [tidyselect::with_vars()] rather than with a data mask;
#' - expressions made only of numeric literals combined with `c()`, `:`,
#'   `(` and unary `-` (e.g. `c(1, 5, 2)`, `1:100`, `-5`, `-(1:10)`);
#' - values, optionally negated: literals (e.g. `5` or `"rs1"`), vectors
#'   injected with `!!`, and the values of function calls (e.g. `which(...)`),
#'   which [select_loci_eval_call()] puts in place of the calls. They are
#'   used as the argument of `all_of()`;
#' - bare variables (e.g. `idx` in `select_loci(x, idx)`), optionally
#'   negated, that are not loci names and hold a vector of positions or
#'   names (see [select_loci_symbol_value()]): their value is used as the
#'   argument of `all_of()`, with the deprecation warning of `tidyselect`
#'   (see [warn_external_vector()]).
#'
#' Anything else, or anything unusual (missing values, zeros in literals,
#' positions out of range, mixed positive and negative positions, names not
#' found by `all_of()`), returns `NULL`, so that [select_loci()] falls back
#' on `tidyselect`, which handles it or generates the appropriate error.
#' @param sel_quo the defused selection (a quosure)
#' @param .data the `gen_tibble`
#' @returns an integer vector of positions in the loci table, or NULL if the
#'   selection could not be resolved
#' @keywords internal
#' @noRd
select_loci_fast <- function(sel_quo, .data) {
  if (rlang::quo_is_missing(sel_quo)) {
    return(NULL)
  }
  expr <- rlang::quo_get_expr(sel_quo)
  n_loci <- count_loci(.data)
  # a negated selection (-x or !x) keeps all loci except the selected ones
  negate <- rlang::is_call(expr, c("-", "!"), n = 1)
  if (negate) {
    expr <- expr[[2]]
  }

  tidy_ns <- c("", "tidyselect", "dplyr")
  if (rlang::is_call(expr, c("all_of", "any_of"), n = 1, ns = tidy_ns)) {
    # any_of() silently drops loci that do not exist, all_of() errors
    strict <- rlang::is_call(expr, "all_of", ns = tidy_ns)
    # tidyselect evaluates the argument of all_of()/any_of() in the caller
    # environment (not in the data mask), so we do the same; if evaluation
    # fails, let tidyselect generate the error
    sel <- tryCatch(
      list(eval(expr[[2]], rlang::quo_get_env(sel_quo))),
      error = function(e) NULL
    )
    if (is.null(sel)) {
      return(NULL)
    }
    pos <- selection_value_positions(sel[[1]], .data, strict = strict)
    if (is.null(pos)) {
      return(NULL)
    }
  } else if (is.null(expr) || is.atomic(expr)) {
    # a value: a literal (e.g. 5 or "rs1"), a vector injected with !!, or the
    # value of a call (see select_loci_eval_call()), which tidyselect uses as
    # it uses the argument of all_of()
    pos <- selection_value_positions(expr, .data, strict = TRUE)
    if (is.null(pos)) {
      return(NULL)
    }
  } else if (is.symbol(expr)) {
    # a bare variable (e.g. select_loci(x, idx)): unless a locus has the same
    # name, tidyselect uses its value as it uses the argument of all_of()
    # (see select_loci_symbol_value())
    sel <- select_loci_symbol_value(expr, sel_quo, .data)
    if (is.null(sel)) {
      return(NULL)
    }
    pos <- selection_value_positions(sel[[1]], .data, strict = TRUE)
    if (is.null(pos)) {
      return(NULL)
    }
    # warn only now that the selection is resolved: when we fall back on
    # tidyselect, it gives the warning itself
    warn_external_vector(as.character(expr), rlang::quo_get_env(sel_quo))
  } else if (rlang::is_call(expr, tidyselect_name_helpers, ns = tidy_ns)) {
    # selection helpers that only need the names of the loci: evaluate them
    # (with their arguments evaluated in the caller environment, as in
    # tidyselect) with the loci names as the current variables, which does
    # not require a data mask
    helper <- get(rlang::call_name(expr), envir = asNamespace("tidyselect"))
    helper_call <- rlang::call2(helper, !!!rlang::call_args(expr))
    pos <- tryCatch(
      tidyselect::with_vars(
        as.character(loci_names(.data)),
        eval(helper_call, rlang::quo_get_env(sel_quo))
      ),
      error = function(e) NULL
    )
    if (!is.numeric(pos) || anyNA(pos)) {
      return(NULL)
    }
    pos <- as.integer(pos)
  } else if (is_numeric_literal_selection(expr)) {
    # only numeric literals, so evaluating the expression has no side effects
    sel <- eval(expr, baseenv())
    if (!is_valid_positions(sel) || length(sel) == 0) {
      return(NULL)
    }
    sel_range <- range(sel)
    if (sel_range[1] > 0 && sel_range[2] <= n_loci) {
      pos <- as.integer(sel)
    } else if (sel_range[2] < 0 && sel_range[1] >= -n_loci && !negate &&
                 is.call(expr)) {
      # c(-1, -2) is equivalent to -c(1, 2) (but negative values in an
      # injected vector are an error in tidyselect)
      pos <- -as.integer(sel)
      negate <- TRUE
    } else {
      # zeros (-0 selects everything in tidyselect), positions out of range,
      # and mixed positive and negative positions (combined as set operations
      # by tidyselect) are left to tidyselect
      return(NULL)
    }
  } else {
    return(NULL)
  }

  # tidyselect returns each selected locus only once, in order of selection
  # (if positions are strictly increasing, they are already unique)
  if (is.unsorted(pos, strictly = TRUE)) {
    pos <- unique(pos)
  }
  if (negate) {
    # keep all other loci, in their original order
    if (length(pos) > 0) {
      pos <- seq_len(n_loci)[-pos]
    } else {
      pos <- seq_len(n_loci)
    }
  }
  pos
}

#' Evaluate an ordinary function call used as a selection
#'
#' `tidyselect` evaluates a call to any function other than its operators
#' (`c()`, `:`, `-`, `!`, `&`, `|`, arithmetic and `~`) as ordinary R code in
#' the caller environment (with its selection helpers available, but not the
#' loci names), and selects with its value. This function does the same for a
#' call (optionally negated with `-` or `!`), e.g. `which(...)`, `sample(...)`
#' or `sort(x)`, and replaces the call with its value, so that the selection
#' can be resolved by [select_loci_fast()]. The call is thus evaluated only
#' once, also when its value is then left to `tidyselect` (e.g. to give an
#' error). Calls that use selection helpers or the `.data` and `.env`
#' pronouns, which only exist within `tidyselect`, and calls that give an
#' error (so that `tidyselect` reports it) are left unchanged.
#' @param sel_quo the defused selection (a quosure)
#' @returns the quosure, with the call replaced by its value
#' @keywords internal
#' @noRd
select_loci_eval_call <- function(sel_quo) {
  if (rlang::quo_is_missing(sel_quo)) {
    return(sel_quo)
  }
  expr <- rlang::quo_get_expr(sel_quo)
  # a negated selection (-x or !x): evaluate the call that is negated
  negate <- rlang::is_call(expr, c("-", "!"), n = 1)
  sel_call <- if (negate) expr[[2]] else expr
  if (!is.call(sel_call)) {
    return(sel_quo)
  }
  # operators that have a meaning in tidyselect
  tidyselect_operators <- c(
    "(", "-", ":", "|", "&", "||", "&&", "!", "+", "*", "/", "^", "~", "c"
  )
  if (is.symbol(sel_call[[1]]) &&
        as.character(sel_call[[1]]) %in% tidyselect_operators) {
    return(sel_quo)
  }
  # selection helpers and pronouns that only exist within tidyselect
  tidyselect_only <- c(
    names(tidyselect::vars_select_helpers), ".data", ".env"
  )
  if (any(all.names(sel_call) %in% tidyselect_only)) {
    return(sel_quo)
  }
  value <- tryCatch(
    list(eval(sel_call, rlang::quo_get_env(sel_quo))),
    error = function(e) NULL
  )
  if (is.null(value)) {
    return(sel_quo)
  }
  # replace the call with its value (keeping the negation)
  value <- value[[1]]
  if (negate) {
    value <- rlang::call2(expr[[1]], value)
  }
  rlang::quo_set_expr(sel_quo, value)
}

#' Find the value of a bare variable used as a selection
#'
#' `tidyselect` resolves a bare symbol (e.g. `idx` in `select_loci(x, idx)`)
#' by looking it up first in a data mask that binds every locus name, and
#' then in the caller environment and its parents; a value that is not a
#' function is then used as the argument of `all_of()` would be (with a
#' deprecation warning, see [warn_external_vector()]). Building that data
#' mask is what makes selecting with a vector given by name quadratic with
#' millions of loci (it creates an R symbol for every locus name), so this
#' function does the same lookup directly, in the same order. It returns
#' `NULL`, so that [select_loci()] falls back on `tidyselect`, when the
#' symbol could be found in the data mask (a locus with the same name, which
#' `tidyselect` selects instead of the variable, or a name starting with a
#' dot), when the variable does not exist or can not be evaluated (so that
#' `tidyselect` gives the error), and when its value is a function (which
#' `tidyselect` treats as a locus name, or as a deprecated predicate).
#' @param sym the symbol
#' @param sel_quo the defused selection (a quosure)
#' @param .data the `gen_tibble`
#' @returns the value of the variable, wrapped in a list (as it may be
#'   `NULL`), or NULL if the selection should be left to `tidyselect`
#' @keywords internal
#' @noRd
select_loci_symbol_value <- function(sym, sel_quo, .data) {
  name <- as.character(sym)
  # besides the loci names, the tidyselect data mask binds the .data and
  # .env pronouns and internal objects (e.g. .top_env), all starting with a
  # dot, which hide variables of the same name
  if (startsWith(name, ".")) {
    return(NULL)
  }
  # the loci are in the data mask, so a locus named like the variable is
  # selected instead (this scans the loci names once, without building a
  # hash table of all of them, see match_loci_names())
  if (!is.na(match_loci_names(name, loci_names(.data)))) {
    return(NULL)
  }
  # look up the variable as tidyselect does, in the caller environment and
  # its parents (which also forces a lazy argument, as tidyselect would)
  value <- tryCatch(
    rlang::env_get(
      rlang::quo_get_env(sel_quo), name,
      default = rlang::missing_arg(), inherit = TRUE
    ),
    error = function(e) rlang::missing_arg()
  )
  if (rlang::is_missing(value) || is.function(value)) {
    return(NULL)
  }
  list(value)
}

#' Warn that a vector given by name was used as a selection
#'
#' `tidyselect` deprecates selecting with a bare variable that holds a vector
#' (e.g. `select_loci(x, idx)`, rather than `select_loci(x, all_of(idx))`,
#' which can not be confused with a locus named `idx`). When
#' [select_loci_fast()] resolves such a selection, this gives the warning
#' that `tidyselect` would give, by letting `tidyselect` itself resolve the
#' same variable (bound to a single position, in a child of the environment
#' of the selection) on a vector of one element: the warning is thus the one
#' of the installed version of `tidyselect` (whose message changes between
#' versions), with its condition class (`lifecycle_warning_deprecated`), and
#' it is shown, or not, as `tidyselect` would show it. As with any soft
#' deprecation by `lifecycle`, it is shown only when the selection is made in
#' the global environment (or in the tests of the package that makes it),
#' and at most once every 8 hours for each variable name. If `tidyselect`
#' ever turns this deprecation into an error, this gives that error.
#' @param name the name of the variable
#' @param user_env the environment in which the selection was made
#' @returns NULL, invisibly
#' @keywords internal
#' @noRd
warn_external_vector <- function(name, user_env) {
  env <- rlang::env(user_env)
  rlang::env_bind(env, !!name := 1L)
  # (the name of the element starts with a dot, unlike the variable, which
  # select_loci_symbol_value() would otherwise leave to tidyselect)
  tidyselect::eval_select(
    rlang::sym(name),
    data = c(.tidypopgen_locus = 1L),
    env = env
  )
  invisible(NULL)
}

#' Find the positions of the loci selected by a value
#'
#' The value of the argument of `all_of()` or `any_of()`, or a value used as
#' a selection (see [select_loci_fast()]), is turned into positions as
#' `tidyselect` does: numbers are positions (zeros are ignored), and strings
#' (or the labels of a factor) are loci names.
#' @param sel the value
#' @param .data the `gen_tibble`
#' @param strict whether all positions and names must exist (as in
#'   `all_of()`); if FALSE, those that do not exist are dropped (as in
#'   `any_of()`)
#' @returns an integer vector of positions, or NULL if `tidyselect` should
#'   handle the value (e.g. to give an error)
#' @keywords internal
#' @noRd
selection_value_positions <- function(sel, .data, strict) {
  n_loci <- count_loci(.data)
  if (is.null(sel)) {
    sel <- integer(0)
  }
  if (is.factor(sel)) {
    sel <- as.character(sel)
  }
  if (is.numeric(sel) && !is.object(sel)) {
    # positions must be whole numbers that can be stored as integers,
    # with no missing values and no negative values
    if (!is_valid_positions(sel) ||
          (length(sel) > 0 && min(sel) < 0)) {
      return(NULL)
    }
    # positions of 0 are ignored by tidyselect
    if (length(sel) > 0 && min(sel) == 0) {
      sel <- sel[sel != 0]
    }
    if (length(sel) > 0 && max(sel) > n_loci) {
      if (strict) {
        return(NULL)
      }
      sel <- sel[sel <= n_loci]
    }
    pos <- as.integer(sel)
  } else if (is.character(sel) && !is.object(sel)) {
    # empty strings are not valid names
    if (anyNA(sel) || any(sel == "")) {
      return(NULL)
    }
    pos <- match_loci_names(sel, loci_names(.data))
    if (anyNA(pos)) {
      if (strict) {
        return(NULL)
      }
      pos <- pos[!is.na(pos)]
    }
  } else {
    return(NULL)
  }
  pos
}

# tidyselect helpers that select based only on the names of the variables
tidyselect_name_helpers <- c(
  "starts_with", "ends_with", "contains", "matches", "num_range",
  "everything", "last_col"
)

#' Check that a numeric vector can be used as positions
#'
#' Positions must be whole numbers that can be stored as integers, with no
#' missing values (tidyselect errors otherwise).
#' @param sel a numeric vector
#' @returns TRUE or FALSE
#' @keywords internal
#' @noRd
is_valid_positions <- function(sel) {
  if (anyNA(sel)) {
    return(FALSE)
  }
  if (is.double(sel)) {
    if (length(sel) > 0 && max(abs(sel)) > .Machine$integer.max) {
      return(FALSE)
    }
    if (any(sel != trunc(sel))) {
      return(FALSE)
    }
  }
  TRUE
}

#' Check that an expression is made only of numeric literals
#'
#' Literals can be combined with `c()`, `:`, `(` and unary `-`; for such
#' expressions, evaluating them in base R gives the same positions as
#' `tidyselect`, as long as the result is all positive or all negative (and
#' contains no zeros, which `tidyselect` treats differently in `-0`).
#' @param expr an expression
#' @returns TRUE or FALSE
#' @keywords internal
#' @noRd
is_numeric_literal_selection <- function(expr) {
  if (is.numeric(expr)) {
    # a number, or a numeric vector injected with !!
    return(TRUE)
  }
  if (!is.call(expr) || !is.symbol(expr[[1]])) {
    return(FALSE)
  }
  fn <- as.character(expr[[1]])
  n_args <- length(expr) - 1
  # binary `-` is a set difference in tidyselect, so only unary `-` is allowed
  valid_call <- (fn == "c") ||
    (fn == ":" && n_args == 2) ||
    (fn %in% c("-", "(") && n_args == 1)
  if (!valid_call) {
    return(FALSE)
  }
  all(vapply(as.list(expr)[-1], is_numeric_literal_selection, logical(1)))
}

#' Match names against the loci names
#'
#' This gives the same result as `match(sel, loci_names)` (i.e. the first
#' occurrence of each name), but it builds the hash table on the (usually
#' short) selection rather than on the (possibly very long) loci names. If
#' loci names are stored as integers, `tidyselect` would compare `sel` to
#' their character representation; we avoid converting millions of integers
#' to strings by converting `sel` to integers instead.
#' @param sel a character vector of names to find
#' @param loci_names the loci names (character or integer)
#' @returns an integer vector of positions (NA for names not found)
#' @keywords internal
#' @noRd
match_loci_names <- function(sel, loci_names) {
  if (is.integer(loci_names)) {
    sel_int <- suppressWarnings(as.integer(sel))
    # only strings that are the canonical representation of an integer
    # (e.g. "30", but not "030" or "3e1") can match an integer name
    sel_int[is.na(sel_int) | as.character(sel_int) != sel] <- NA_integer_
    sel <- sel_int
  }
  # candidate loci, i.e. the loci whose name is in the selection (names
  # that are not valid, now NA, can not match a missing locus name)
  candidates <- which(loci_names %in% sel[!is.na(sel)])
  # and now match among the candidates (keeping the first occurrence)
  candidates[match(sel, loci_names[candidates], incomparables = NA)]
}
