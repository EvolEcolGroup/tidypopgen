#' Write site frequency spectra for dadi or fastsimcoal2
#'
#' Write site frequency spectra in the formats used by `dadi` (and `moments`)
#' and by `fastsimcoal2`, as written by `easySFS`:
#' - for an `sfs` object (as returned by [pop_sfs()], [pairwise_pop_sfs()] and
#'   [nwise_pop_sfs()] with `type = "array"`, or by the `bootstrap_*`
#'   functions, e.g. [bootstrap_pop_sfs()]), write that spectrum to `file`;
#' - for a grouped [`gen_tibble`], compute all the spectra that `easySFS`
#'   writes (the 1D spectrum of each population, the joint spectrum of each
#'   pair of populations, and the joint spectrum of all populations, in a
#'   single pass over the genotypes), and write them in `out_dir`, with the
#'   file names of `easySFS`, optionally with bootstrap replicates. Missing
#'   genotypes are dealt with by projecting each population down to `n_proj`
#'   allele copies (see [pop_sfs()]).
#'
#' @details
#' **dadi format** (`format = "dadi"`): a line with the dimensions of the
#' spectrum (the number of allele copies of each population plus 1),
#' `folded` or `unfolded` and the population labels in quotes; a line with
#' the cells of the spectrum; and a line with the mask (1 for masked cells,
#' see [pop_sfs()]). In dadi files, the cells are listed with the number of
#' copies in the last population changing fastest. With `fastest = "first"`,
#' the populations are written in reverse order, so that the number of
#' copies in the first population changes fastest and the file is still a
#' valid dadi file.
#'
#' **fastsimcoal2 format** (`format = "fastsimcoal2"`): for a 1D spectrum,
#' as a `_DAFpop0.obs` (unfolded) or `_MAFpop0.obs` (folded) file; for a 2D
#' spectrum, as a `_jointDAFpop1_0.obs` or `_jointMAFpop1_0.obs` file, a
#' matrix with the second population in the rows (labelled `d1_0`, `d1_1`,
#' ...) and the first population in the columns (labelled `d0_0`, `d0_1`,
#' ...); for a spectrum of three or more populations, as a multidimensional
#' spectrum (`_MSFS.obs` file). For an `sfs` object,
#' `format = "fastsimcoal2_multidimensional"` writes the multidimensional
#' format whatever the number of populations. In the multidimensional
#' format, the cells are listed in the order set by
#' `fastest`: "last" (the default) has the number of copies in the last
#' population changing fastest, as written by `easySFS`; "first" has the
#' first population changing fastest. Please check which order your version
#' of `fastsimcoal2` expects.
#'
#' **Monomorphic sites.** As in `easySFS`, `total_length` gives the total
#' number of sites sequenced (including monomorphic sites that are not in
#' the `gen_tibble`): the first cell of the spectrum (sites with no copies of
#' the alternate or minor allele) is then set to `total_length` minus the
#' number of sites in all the other cells.
#'
#' **The spectra of a `gen_tibble`** are written in `out_dir`, with the file
#' names of `easySFS` (`pop` stands for the population labels, `p` for their
#' projections, and `M` is replaced by `D` for unfolded spectra):
#' - dadi: `pop-p.sfs` for each population, `pop1-pop2.sfs` for each pair,
#'   and `pop1-pop2-...-popK.sfs` for all populations;
#' - fastsimcoal2: `pop_MAFpop0.obs` for each population,
#'   `prefix_jointMAFpopj_i.obs` for each pair (where `i` and `j` are the
#'   0-based indices of the two populations, so that the second population of
#'   the pair is in the rows), and `prefix_MSFS.obs` for all populations.
#'
#' `easySFS` writes these files in the sub-directories `dadi` and
#' `fastsimcoal2` of its output directory: use, for example,
#' `out_dir = file.path(dir, "dadi")` with `format = "dadi"` and
#' `out_dir = file.path(dir, "fastsimcoal2")` with `format = "fastsimcoal2"`
#' to get the same layout. Populations are in the order of the groups (see
#' [dplyr::group_keys()]): to change their order, which matters for `dadi`
#' and `fastsimcoal2` (and is set with `--order` in `easySFS`), change the
#' order of the levels of the grouping variable (e.g. with [factor()]).
#'
#' With `n_bootstrap` replicates of a block bootstrap (see
#' [bootstrap_pop_sfs()]), the spectra of each replicate are also written,
#' with the same file names, in the sub-directories `bootstrap_1`,
#' `bootstrap_2`, ... of `out_dir`, numbered with as many digits as
#' `n_bootstrap` (e.g. `bootstrap_001` to `bootstrap_100`), which can be
#' used to estimate confidence intervals with `fastsimcoal2`. Replicates use
#' the same projections as the spectra of all loci, and `total_length` and
#' `round` are applied to each of them. To draw the same replicates in both
#' formats, set the same seed (with [set.seed()]) before writing each of
#' them.
#'
#' @inheritParams nwise_pop_sfs
#' @param x an `sfs` object, or a grouped [`gen_tibble`] (as obtained by
#'   using [dplyr::group_by()]).
#' @param file for an `sfs` object, the name of the file to write.
#' @param out_dir for a `gen_tibble`, the directory in which to write the
#'   files, or `NULL` (the default) to use the location and name of the
#'   backing files, followed by `_dadi` or `_fastsimcoal2`.
#' @param format the format of the files: `"dadi"` or `"fastsimcoal2"`; for
#'   an `sfs` object, also `"fastsimcoal2_multidimensional"` (the
#'   multidimensional format of fastsimcoal2, whatever the number of
#'   populations). See details.
#' @param fastest the order of the cells in dadi files and in the
#'   multidimensional fastsimcoal2 format: "last" (the default) keeps the order
#'   of the populations, so that the number of allele copies in the last
#'   population changes fastest, as in `dadi` and `easySFS`; "first" makes the
#'   number of copies in the first population change fastest. See details.
#' @param total_length the total number of sites sequenced, including
#'   monomorphic ones, used to set the number of monomorphic sites (the first
#'   cell). If `NULL` (the default), the spectra are written as they are.
#' @param round boolean, whether to round the number of sites in each cell to
#'   the nearest integer (as `easySFS` with `--dtype int`). Defaults to
#'   `FALSE`.
#' @param overwrite boolean, whether to overwrite files that already exist.
#'   Defaults to `FALSE`.
#' @param prefix for a `gen_tibble` written in fastsimcoal2 format, the
#'   prefix of the names of the 2D and multidimensional files, or `NULL` (the
#'   default) to use the name of the backing files.
#' @param n_bootstrap for a `gen_tibble`, the number of bootstrap replicates
#'   to write (0, the default, for none).
#' @param bootstrap_block_size,bootstrap_size_unit for a `gen_tibble`, the
#'   size of the blocks of loci resampled for the bootstrap, and its unit
#'   ("snp", the default, "cM" or "bp"); see [bootstrap_pop_sfs()].
#' @param max_cells for a `gen_tibble`, the maximum number of cells of the
#'   joint spectrum of all populations; an error is raised for larger
#'   spectra.
#' @param ... further arguments passed to the methods.
#' @returns the path of the file (for an `sfs` object) or the paths of the
#'   files written (for a `gen_tibble`), invisibly.
#' @seealso [pop_sfs()], [pop_sfs_preview()], [bootstrap_pop_sfs()]
#' @export
#' @examplesIf all(rlang::is_installed(c("RhpcBLASctl", "data.table")))
#' \dontshow{
#' data.table::setDTthreads(2)
#' RhpcBLASctl::blas_set_num_threads(2)
#' RhpcBLASctl::omp_set_num_threads(2)
#' }
#' example_gt <- load_example_gt("grouped_gen_tbl")
#'
#' # write a single spectrum
#' sfs_2d <- example_gt %>%
#'   pairwise_pop_sfs(n_proj = 2, type = "array")
#' dadi_file <- tempfile(fileext = ".sfs")
#' write_sfs(sfs_2d$pop1.pop2, dadi_file)
#' readLines(dadi_file)
#' fastsimcoal2_file <- tempfile(fileext = "_jointMAFpop1_0.obs")
#' write_sfs(sfs_2d$pop1.pop2, fastsimcoal2_file, format = "fastsimcoal2")
#' readLines(fastsimcoal2_file)
#'
#' # write all the spectra written by easySFS, in the same directories
#' out_dir <- file.path(tempdir(), "example_sfs")
#' dadi_files <- example_gt %>%
#'   write_sfs(out_dir = file.path(out_dir, "dadi"), n_proj = 2)
#' basename(dadi_files)
#' fastsimcoal2_files <- example_gt %>%
#'   write_sfs(
#'     out_dir = file.path(out_dir, "fastsimcoal2"), format = "fastsimcoal2",
#'     n_proj = 2, prefix = "example"
#'   )
#' basename(fastsimcoal2_files)
#'
#' # with 3 bootstrap replicates, resampling blocks of 2 SNPs
#' set.seed(1)
#' fastsimcoal2_files <- example_gt %>%
#'   write_sfs(
#'     out_dir = file.path(out_dir, "fastsimcoal2"), format = "fastsimcoal2",
#'     n_proj = 2, prefix = "example", n_bootstrap = 3,
#'     bootstrap_block_size = 2, overwrite = TRUE
#'   )
#' list.files(file.path(out_dir, "fastsimcoal2"))
write_sfs <- function(x, ...) {
  UseMethod("write_sfs", x)
}

#' @rdname write_sfs
#' @export
write_sfs.sfs <- function(
  x,
  file,
  format = c("dadi", "fastsimcoal2", "fastsimcoal2_multidimensional"),
  fastest = c("last", "first"),
  total_length = NULL,
  round = FALSE,
  overwrite = FALSE,
  ...
) {
  rlang::check_dots_empty()
  format <- match.arg(format)
  fastest <- match.arg(fastest)
  .sfs_check_flag(round, "round")
  .sfs_check_flag(overwrite, "overwrite")
  if (!overwrite && file.exists(file)) {
    stop(
      "file ", file,
      " already exists; use 'overwrite = TRUE' to overwrite it"
    )
  }
  x <- .sfs_finalise(x, total_length = total_length, round = round)
  n_pops <- length(dim(x))
  if (format == "dadi") {
    lines <- .sfs_dadi_lines(x, fastest = fastest)
  } else if (format == "fastsimcoal2_multidimensional" || n_pops > 2) {
    lines <- .sfs_fastsimcoal2_msfs_lines(x, fastest = fastest)
  } else if (n_pops == 1) {
    lines <- .sfs_fastsimcoal2_1d_lines(x)
  } else {
    lines <- .sfs_fastsimcoal2_2d_lines(x, pop_index = c(0, 1))
  }
  writeLines(lines, file)
  invisible(file)
}

#' @rdname write_sfs
#' @export
write_sfs.gen_tbl <- function(
  x,
  out_dir = NULL,
  format = c("dadi", "fastsimcoal2"),
  n_proj = NULL,
  prefix = NULL,
  folded = TRUE,
  ancestral = NULL,
  total_length = NULL,
  round = FALSE,
  fastest = c("last", "first"),
  n_bootstrap = 0,
  bootstrap_block_size = NULL,
  bootstrap_size_unit = c("snp", "cM", "bp"),
  max_cells = 1e7,
  overwrite = FALSE,
  n_cores = bigstatsr::nb_cores(),
  block_size = bigstatsr::block_size(nrow(x), 1),
  ...
) {
  rlang::check_dots_empty()
  format <- match.arg(format)
  .sfs_write_set_files(
    x,
    format = format, out_dir = out_dir, n_proj = n_proj, prefix = prefix,
    folded = folded, ancestral = ancestral, total_length = total_length,
    round = round, fastest = fastest, n_bootstrap = n_bootstrap,
    bootstrap_block_size = bootstrap_block_size,
    bootstrap_size_unit = bootstrap_size_unit, max_cells = max_cells,
    overwrite = overwrite, n_cores = n_cores, block_size = block_size
  )
}

# a data frame that has lost the class of its tibble (e.g. a gen_tibble
# returned by pick() or reframe()): restore the class and dispatch again
#' @export
#' @noRd
write_sfs.data.frame <- function(x, ...) {
  if (!"genotypes" %in% names(x)) {
    # not a tibble of genotypes (e.g. the tidy output of pop_sfs()): explain
    # what x should be
    return(NextMethod())
  }
  write_sfs(restore_genotypes_tbl(x, "write_sfs"), ...)
}

#' @rdname write_sfs
#' @export
write_sfs.default <- function(x, ...) {
  stop(
    "`x` should be an sfs object (as returned by pop_sfs(), ",
    "pairwise_pop_sfs() or nwise_pop_sfs() with type = \"array\"), or a ",
    "grouped gen_tibble"
  )
}

# the work of write_sfs() for a gen_tibble: compute the spectra of each
# population, each pair and all populations (and their bootstrap
# replicates), and write them in the given format ("dadi" or
# "fastsimcoal2")
.sfs_write_set_files <- function(.x, format, out_dir, n_proj, prefix, folded,
                                 ancestral, total_length, round, fastest,
                                 n_bootstrap, bootstrap_block_size,
                                 bootstrap_size_unit, max_cells, overwrite,
                                 n_cores, block_size) {
  if (!inherits(.x, "gen_tbl") || !inherits(.x, "grouped_df")) {
    stop("`x` should be a grouped gen_tibble (see dplyr::group_by())")
  }
  fastest <- match.arg(fastest, c("last", "first"))
  bootstrap_size_unit <- match.arg(bootstrap_size_unit, c("snp", "cM", "bp"))
  .sfs_check_folded(folded)
  .sfs_check_flag(round, "round")
  .sfs_check_flag(overwrite, "overwrite")
  if (!identical(n_bootstrap, 0) && !identical(n_bootstrap, 0L)) {
    .sfs_check_bootstrap_args(n_bootstrap, bootstrap_block_size)
  }
  backing_name <- bigstatsr::sub_bk(attr(.x$genotypes, "fbm")$backingfile)
  if (is.null(out_dir)) {
    out_dir <- paste0(
      backing_name,
      if (format == "dadi") "_dadi" else "_fastsimcoal2"
    )
  }
  if (is.null(prefix)) {
    prefix <- basename(backing_name)
  }
  # the prefix is part of the file names (in fastsimcoal2 format), so a
  # path separator would write files outside out_dir
  is_name <- is.character(prefix) && length(prefix) == 1 && !is.na(prefix) &&
    nzchar(prefix) && !grepl("[/\\\\]", prefix)
  if (!is_name) {
    stop(
      "`prefix` should be a single string without '/' or '\\', as it is ",
      "used in file names"
    )
  }
  data <- .sfs_set_ancestral(.sfs_data(.x), ancestral, folded = folded)
  labels <- data$labels
  if (any(grepl("[/\\\\]", labels))) {
    stop(
      "group labels cannot contain '/' or '\\', as they are used in ",
      "file names"
    )
  }
  n_proj <- .sfs_n_proj(n_proj, data,
    n_cores = n_cores,
    block_size = block_size
  )
  n_cells <- prod(as.numeric(n_proj) + 1)
  if (n_cells > max_cells) {
    stop(
      "the spectrum of all populations would have ", format(n_cells),
      " cells (more than `max_cells`); use fewer populations or smaller ",
      "projections, or increase `max_cells`"
    )
  }

  # the spectra (1D, each pair, and all populations), as in easySFS
  n_groups <- data$n_groups
  if (n_groups >= 2) {
    pairs <- utils::combn(n_groups, 2)
  } else {
    pairs <- matrix(integer(0), nrow = 2)
  }
  spectra_def <- c(
    as.list(seq_len(n_groups)),
    lapply(seq_len(ncol(pairs)), function(i) pairs[, i]),
    list(seq_len(n_groups))
  )
  type <- c(
    rep("1d", n_groups), rep("2d", ncol(pairs)), "multi"
  )

  # the file names of easySFS (relative to out_dir, and to the directory of
  # each bootstrap replicate)
  if (format == "dadi") {
    pair_files <- character(0)
    if (ncol(pairs) > 0) {
      pair_files <- paste0(
        labels[pairs[1, ]], "-", labels[pairs[2, ]],
        ".sfs"
      )
    }
    set_files <- c(
      paste0(labels, "-", n_proj, ".sfs"),
      pair_files,
      paste0(paste(labels, collapse = "-"), ".sfs")
    )
  } else {
    m_or_d <- if (folded) "M" else "D"
    pair_files <- character(0)
    if (ncol(pairs) > 0) {
      pair_files <- paste0(
        prefix, "_joint", m_or_d, "AFpop", pairs[2, ] - 1, "_",
        pairs[1, ] - 1, ".obs"
      )
    }
    set_files <- c(
      paste0(labels, "_", m_or_d, "AFpop0.obs"),
      pair_files,
      paste0(prefix, "_MSFS.obs")
    )
  }
  # two different spectra cannot be written to the same file (e.g. in dadi
  # format, the spectrum of a group "pop1" projected to 2 copies and the
  # spectrum of the pair of "pop1" and a group "2" are both "pop1-2.sfs").
  # With two groups, the spectra of the pair and of all groups are the same,
  # and share their file, as in easySFS
  spectra_key <- vapply(spectra_def, paste, character(1), collapse = ",")
  is_shared <- tapply(spectra_key, set_files, function(k) {
    length(unique(k)) > 1
  })
  if (any(is_shared)) {
    stop(
      "several spectra would be written to the same file (",
      paste(names(is_shared)[is_shared], collapse = ", "),
      "); rename the groups"
    )
  }
  # the directory of each replicate (e.g. bootstrap_001)
  rep_dirs <- character(0)
  if (n_bootstrap > 0) {
    rep_dirs <- sprintf(
      "bootstrap_%0*d", nchar(n_bootstrap),
      seq_len(n_bootstrap)
    )
  }
  files <- unique(file.path(out_dir, c(
    set_files,
    as.vector(outer(set_files, rep_dirs, function(f, d) file.path(d, f)))
  )))
  if (!overwrite && any(file.exists(files))) {
    stop(
      "some files already exist in ", out_dir,
      "; use 'overwrite = TRUE' to overwrite them"
    )
  }

  # draw the blocks of the bootstrap replicates before writing any file, so
  # that we stop early if they cannot be drawn
  if (n_bootstrap > 0) {
    blocks <- bootstrap_blocks(data$loci,
      n_bootstrap = n_bootstrap,
      block_size = bootstrap_block_size,
      size_unit = bootstrap_size_unit
    )
  }

  configs <- .sfs_config_counts(data,
    spectra = spectra_def,
    min_n = n_proj,
    n_cores = n_cores,
    block_size = block_size
  )
  spectra <- .sfs_from_configs(configs, spectra_def, n_proj,
    folded = folded,
    n_cores = n_cores
  )
  # write the spectra of a set (all loci, or a bootstrap replicate)
  write_set <- function(spectra, set_dir) {
    spectra <- lapply(spectra, .sfs_finalise,
      total_length = total_length,
      round = round
    )
    dir.create(set_dir, recursive = TRUE, showWarnings = FALSE)
    for (i in seq_along(spectra)) {
      fs <- spectra[[i]]
      if (format == "dadi") {
        lines <- .sfs_dadi_lines(fs, fastest = fastest)
      } else {
        lines <- switch(type[i],
          "1d" = .sfs_fastsimcoal2_1d_lines(fs),
          "2d" = .sfs_fastsimcoal2_2d_lines(fs,
            pop_index = spectra_def[[i]] - 1
          ),
          "multi" = .sfs_fastsimcoal2_msfs_lines(fs, fastest = fastest)
        )
      }
      writeLines(lines, file.path(set_dir, set_files[i]))
    }
  }
  write_set(spectra, out_dir)

  # bootstrap replicates: the configurations are counted in each block, in a
  # second pass over the genotypes, and the replicates are projected in
  # chunks, to limit the memory used
  if (n_bootstrap > 0) {
    configs <- .sfs_config_counts(data,
      spectra = spectra_def,
      min_n = n_proj,
      n_cores = n_cores,
      block_size = block_size,
      block_ids = blocks$block_id
    )
    n_cells <- sum(vapply(spectra_def, function(x) {
      prod(as.numeric(n_proj[x]) + 1)
    }, numeric(1)))
    # (the maximum number of cells of a chunk can be set with an option, for
    # testing)
    max_chunk_cells <- getOption("tidypopgen.sfs_max_chunk_cells", 1e7)
    chunk_size <- max(1, floor(max_chunk_cells / n_cells))
    chunks <- split(
      seq_len(n_bootstrap),
      ceiling(seq_len(n_bootstrap) / chunk_size)
    )
    for (chunk in chunks) {
      projected <- .sfs_project(configs, spectra_def, n_proj,
        folded = folded,
        n_cores = n_cores,
        weights = blocks$weights[chunk, , drop = FALSE]
      )
      for (i in seq_along(chunk)) {
        write_set(
          lapply(projected, .sfs_replicate, i_rep = i),
          file.path(out_dir, rep_dirs[chunk[i]])
        )
      }
    }
  }
  invisible(files)
}

# check that an argument is TRUE or FALSE
.sfs_check_flag <- function(x, name) {
  if (!is.logical(x) || length(x) != 1 || is.na(x)) {
    stop("`", name, "` should be TRUE or FALSE")
  }
  invisible(TRUE)
}

# Set the number of monomorphic sites (the first cell) of a spectrum from the
# total number of sites sequenced, and round the cells, as easySFS does
# before writing (--total-length and --dtype int)
.sfs_finalise <- function(fs, total_length = NULL, round = FALSE) {
  if (!is.null(total_length)) {
    is_number <- is.numeric(total_length) && length(total_length) == 1 &&
      is.finite(total_length)
    if (!is_number) {
      stop("`total_length` should be a single finite number")
    }
    n_sites <- sum(unclass(fs))
    if (total_length <= n_sites) {
      stop(
        "`total_length` (", total_length, ") should be larger than the ",
        "number of sites in the spectrum (", base::round(n_sites), ")"
      )
    }
    fs[1] <- fs[1] + total_length - n_sites
  }
  if (round) {
    fs[] <- base::round(as.vector(unclass(fs)))
  }
  fs
}

# format numbers as dadi (numpy.savetxt with fmt = "%.16g")
.sfs_format_num <- function(x) {
  sprintf("%.16g", x)
}

# the cells (and mask) of a spectrum in the order of fastest
.sfs_ordered_values <- function(x, fastest) {
  dims <- dim(x)
  values <- array(as.vector(unclass(x)), dim = dims)
  mask <- attr(x, "mask")
  if (fastest == "last") {
    perm <- rev(seq_along(dims))
    values <- aperm(values, perm)
    mask <- aperm(mask, perm)
  }
  list(values = as.vector(values), mask = as.vector(mask))
}

# the lines of a dadi file (Spectrum.to_file()): the dimensions, folding and
# population labels; the cells (the last population changing fastest); and
# the mask. With fastest = "first", the populations are written in reverse
# order, so that the first population changes fastest
.sfs_dadi_lines <- function(x, fastest = "last") {
  dims <- dim(x)
  pops <- attr(x, "pop_ids")
  # the labels are written between double quotes on the first line, and
  # dadi splits that line on '"', so they cannot contain either
  if (any(grepl("[\"\r\n]", pops))) {
    stop(
      "population labels cannot contain '\"' or line breaks in dadi files"
    )
  }
  cells <- .sfs_ordered_values(x, fastest)
  if (fastest == "first") {
    dims <- rev(dims)
    pops <- rev(pops)
  }
  c(
    paste0(
      paste(dims, collapse = " "), " ",
      if (attr(x, "folded")) "folded" else "unfolded",
      paste0(" \"", pops, "\"", collapse = "")
    ),
    paste(.sfs_format_num(cells$values), collapse = " "),
    paste(as.integer(cells$mask), collapse = " ")
  )
}

# the lines of a 1D fastsimcoal2 file (as written by easySFS)
.sfs_fastsimcoal2_1d_lines <- function(x) {
  p <- dim(x)[1] - 1
  c(
    "1 observation",
    paste0("d0_", 0:p, collapse = "\t"),
    paste(.sfs_format_num(as.vector(unclass(x))), collapse = " "),
    ""
  )
}

# the lines of a 2D fastsimcoal2 file (as written by easySFS): the second
# population in the rows, the first in the columns, labelled with their
# (0-based) indices pop_index
.sfs_fastsimcoal2_2d_lines <- function(x, pop_index) {
  dims <- dim(x)
  values <- matrix(as.vector(unclass(x)), nrow = dims[1])
  rows <- vapply(seq_len(dims[2]), function(r) {
    paste0(
      "d", pop_index[2], "_", r - 1, "\t",
      paste(.sfs_format_num(values[, r]), collapse = " ")
    )
  }, character(1))
  c(
    "1 observation",
    paste0("\t", paste0("d", pop_index[1], "_", 0:(dims[1] - 1),
      collapse = "\t"
    )),
    rows
  )
}

# the lines of a multidimensional fastsimcoal2 file (as written by easySFS),
# with the cells in the order set by fastest
.sfs_fastsimcoal2_msfs_lines <- function(x, fastest = "last") {
  dims <- dim(x)
  cells <- .sfs_ordered_values(x, fastest)
  c(
    "1 observations. No. of demes and sample sizes are on next line.",
    paste0(length(dims), "\t", paste(dims - 1, collapse = " ")),
    paste(.sfs_format_num(cells$values), collapse = " "),
    ""
  )
}
